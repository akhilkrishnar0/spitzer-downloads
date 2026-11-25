"""Download original SEIP mosaics and display diagnostic figures."""
from pathlib import Path
import csv
import re
import warnings
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.nddata.utils import NoOverlapError
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
import astropy.units as u
import matplotlib.pyplot as plt
from spitzer_query import SpitzerQuery, classify, BANDS, PRODUCTS


class SpitzerPipeline:
    def __init__(self, base_dir='spitzer_data', galaxy_centers=None):
        self.base_dir = Path(base_dir)
        self.galaxy_centers = galaxy_centers or {}
        self.query = SpitzerQuery()

    def get_coord(self, target):
        if target in self.galaxy_centers:
            ra, dec = self.galaxy_centers[target]
            return SkyCoord(ra, dec, unit='deg', frame='icrs')
        try:
            return SkyCoord.from_name(target).icrs
        except Exception as exc:
            raise ValueError(f'Cannot resolve {target}; supply galaxy_centers={{"{target}": (RA, Dec)}} in degrees') from exc

    @staticmethod
    def clean_name(target):
        name = re.sub(r'[^A-Za-z0-9_.-]', '', target)
        if not name or name in ('.', '..'):
            raise ValueError('Invalid target name')
        return name

    def run(self, targets, bands=('3.6', '4.5', '5.8', '8.0', '24'),
            radius_arcmin=6, skip_download=True, do_plot=True,
            force_download=False, products=('science', 'uncertainty', 'coverage')):
        if isinstance(targets, str):
            targets = [targets]
        if isinstance(bands, (str, int, float)):
            bands = [bands]
        bands = tuple(str(b) for b in bands)
        if not set(bands) <= set(BANDS.values()):
            raise ValueError(f'Available bands: {tuple(BANDS.values())}')
        if isinstance(products, str):
            products = [products]
        if not set(products) <= set(PRODUCTS.values()):
            raise ValueError(f'Unknown products: {products}')
        return {target: self._run_one(target, bands, radius_arcmin, skip_download,
                                      do_plot, force_download, products) for target in targets}

    def _run_one(self, target, bands, radius, skip, do_plot, force, products):
        coord = self.get_coord(target)
        root = self.base_dir / self.clean_name(target)
        root.mkdir(parents=True, exist_ok=True)
        table = self.query.search(coord.ra.deg, coord.dec.deg, radius)
        table.write(root/f'archive_query_{self.clean_name(target)}.ecsv', format='ascii.ecsv', overwrite=True)
        groups = {}
        records = []
        seen = set()
        for row in table:
            url = row['access_url']
            if np.ma.is_masked(url):
                continue
            url = url.decode() if isinstance(url, bytes) else str(url)
            info = classify(url)
            if info is None or url in seen:
                continue
            seen.add(url)
            band, kind, prefix, name = info
            if band not in bands or kind not in products:
                continue
            # Preserve matching ancillary files by exact mosaic prefix.
            key = (band, prefix)
            path = root/'downloads'/f'{band}um'/prefix/name
            record = dict(band=band, mosaic=prefix, product=kind, url=url,
                          path=str(path), status='downloaded', center_pixel='')
            try:
                self.query.download(url, path, overwrite=force or not skip)
                data, wcs, header = self.read_image(path)
                x, y = skycoord_to_pixel(coord, wcs, origin=0)
                inside = bool(np.isfinite(x) and np.isfinite(y) and
                              -.5 <= x < data.shape[1]-.5 and -.5 <= y < data.shape[0]-.5)
                record['center_pixel'] = 'inside' if inside else 'outside'
                record['shape'] = f'{data.shape[0]}x{data.shape[1]}'
                record['bunit'] = header.get('BUNIT', '')
                groups.setdefault(key, {})[kind] = path
                print(f'{target}: {band} um {kind} {prefix} ({record["center_pixel"]})')
            except Exception as exc:
                record['status'] = f'failed: {exc}'
                warnings.warn(f'{target}: {name}: {exc}')
            records.append(record)
        fields = ('band', 'mosaic', 'product', 'url', 'path', 'status', 'center_pixel', 'shape', 'bunit')
        with (root/f'inventory_{self.clean_name(target)}.csv').open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(records)
        missing = [band for band in bands if not any(k[0] == band and 'science' in g for k, g in groups.items())]
        if 'science' in products and missing:
            print(f'{target}: no science mosaics downloaded for {missing}')
        for key, group in groups.items():
            absent = set(products) - set(group)
            if absent:
                print(f'{target}: {key} missing {sorted(absent)}')
        if do_plot:
            self.plot(target, coord, groups, products, radius)
        return groups

    @staticmethod
    def read_image(path):
        with fits.open(path, memmap=False) as hdul:
            hdu = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
            wcs = WCS(hdu.header).celestial
            if not wcs.has_celestial:
                raise ValueError(f'No celestial WCS in {path}')
            return hdu.data.copy(), wcs, hdu.header.copy()

    def crop(self, target, size_arcmin=3, center=None, bands=None, do_plot=True):
        """Crop existing matched science/uncertainty mosaics without resampling.

        size_arcmin is the full width of the square, not its radius.
        center is an optional (RA, Dec) pair in ICRS degrees.
        Output is a new directory for this centre and size; raw FITS are unchanged.
        """
        size = float(size_arcmin)
        if not np.isfinite(size) or size <= 0:
            raise ValueError("size_arcmin must be positive and finite")
        coord = self.get_coord(target) if center is None else SkyCoord(
            center[0], center[1], unit="deg", frame="icrs")
        if bands is None:
            bands = tuple(BANDS.values())
        elif isinstance(bands, (str, int, float)):
            bands = [bands]
        bands = tuple(str(b) for b in bands)
        if not set(bands) <= set(BANDS.values()):
            raise ValueError(f"Available bands: {tuple(BANDS.values())}")
        root = self.base_dir / self.clean_name(target)
        downloads = root / "downloads"
        if not downloads.is_dir():
            raise FileNotFoundError(f"{downloads}: run the download pipeline first")
        groups = {}
        for path in sorted(downloads.rglob("*.fits*")):
            if not path.is_file():
                continue
            info = classify(path.name)
            if info is None:
                continue
            band, kind, prefix, _ = info
            if band not in bands or kind not in ("science", "uncertainty"):
                continue
            group = groups.setdefault((band, prefix), {})
            if kind in group:
                raise ValueError(f"Duplicate {kind} files for {band} {prefix}")
            group[kind] = path
        if not groups:
            raise FileNotFoundError(f"No requested science/uncertainty FITS under {downloads}")
        out = root / "crops" / f"{size:g}arcmin_RA{coord.ra.deg:.6f}_Dec{coord.dec.deg:+.6f}"
        out.mkdir(parents=True, exist_ok=True)
        results = {}
        records = []
        for (band, prefix), group in groups.items():
            record = {"band": band, "mosaic": prefix, "status": "", "science": "", "uncertainty": ""}
            try:
                if not {"science", "uncertainty"} <= group.keys():
                    raise ValueError("Missing science or uncertainty partner")
                image, iw, ih = self.read_image(group["science"])
                unc, uw, uh = self.read_image(group["uncertainty"])
                if image.shape != unc.shape:
                    raise ValueError("Science and uncertainty shapes differ")
                ny, nx = image.shape
                xx = np.array([0, nx-1, 0, nx-1, (nx-1)/2])
                yy = np.array([0, 0, ny-1, ny-1, (ny-1)/2])
                ra, dec = iw.pixel_to_world_values(xx, yy)
                ux, uy = uw.world_to_pixel_values(ra, dec)
                if not (np.allclose(ux, xx, atol=1e-3, rtol=0) and
                        np.allclose(uy, yy, atol=1e-3, rtol=0)):
                    raise ValueError("Science and uncertainty WCS differ; no resampling performed")
                if iw.has_distortion or uw.has_distortion:
                    # Cutout2D does not handle FITS distortion lookup tables.
                    if any(getattr(w, k, None) is not None for w in (iw, uw)
                           for k in ("cpdis1", "cpdis2", "det2im1", "det2im2")):
                        raise ValueError("Unsupported WCS distortion lookup table")
                cut = Cutout2D(image.astype(float, copy=False), coord,
                               (size*u.arcmin, size*u.arcmin), wcs=iw,
                               mode="partial", fill_value=np.nan, copy=True)
                # Copy exactly the same native pixel rectangle for uncertainty.
                cropped_unc = np.full(cut.data.shape, np.nan, dtype=float)
                cropped_unc[cut.slices_cutout] = unc[cut.slices_original]
                folder = out / f"{band}um" / prefix
                folder.mkdir(parents=True, exist_ok=True)
                cropped = {}
                for kind, data, header in (("science", cut.data, ih),
                                            ("uncertainty", cropped_unc, uh)):
                    h = header.copy()
                    # A pure pixel translation: preserve the original calibration
                    # and WCS representation (including CD or PC and SIP terms).
                    dy = cut.slices_original[0].start - cut.slices_cutout[0].start
                    dx = cut.slices_original[1].start - cut.slices_cutout[1].start
                    h["CRPIX1"] = float(header.get("CRPIX1", 0.0)) - dx
                    h["CRPIX2"] = float(header.get("CRPIX2", 0.0)) - dy
                    if "IMAGEW" in h:
                        h["IMAGEW"] = data.shape[1]
                    if "IMAGEH" in h:
                        h["IMAGEH"] = data.shape[0]
                    h.add_history(f"Native pixel cutout: requested {size:g} arcmin square; no resampling.")
                    h.add_history("Outside original footprint is NaN.")
                    path = folder / f"{prefix}_{band}um_{kind}_{self.clean_name(target)}.fits"
                    fits.writeto(path, data, h, overwrite=True)
                    cropped[kind] = path
                    record[kind] = str(path)
                results[(band, prefix)] = cropped
                partial = any(s.start != 0 or s.stop != n
                              for s, n in zip(cut.slices_cutout, cut.data.shape))
                record["status"] = "partial (NaN padded)" if partial else "cropped"
                print(f"{target}: {band} um {prefix}: {record['status']} {cut.data.shape}")
            except (ValueError, OSError, NoOverlapError) as exc:
                record["status"] = f"skipped: {exc}"
                warnings.warn(f"{target}: {band} um {prefix}: {exc}")
            records.append(record)
        with (out / f"inventory_{self.clean_name(target)}.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=("band", "mosaic", "status", "science", "uncertainty"))
            writer.writeheader()
            writer.writerows(records)
        if do_plot and results:
            self.plot(target, coord, results, ("science", "uncertainty"), size,
                      title_suffix=f"{size:g} arcmin crop",
                      filename_suffix=f"_crop_{size:g}arcmin_RA{coord.ra.deg:.6f}_Dec{coord.dec.deg:+.6f}")
        print(f"Crop output: {out}")
        return results

    def plot(self, target, coord, groups, products, radius_arcmin, output_dir=None, title_suffix='native mosaics', filename_suffix=''):
        out = Path(output_dir) if output_dir is not None else self.base_dir/self.clean_name(target)/'plots'
        out.mkdir(parents=True, exist_ok=True)
        # One plot per mosaic prefix, not arbitrary pairing across observations.
        prefixes = sorted({key[1] for key in groups})
        for prefix in prefixes:
            available = [(band, groups[(band, prefix)]) for band in BANDS.values() if (band, prefix) in groups]
            rows = list(products)
            fig = plt.figure(figsize=(4.2*len(available), 3.7*len(rows)))
            try:
                for j, (band, group) in enumerate(available):
                    for i, kind in enumerate(rows):
                        idx = i*len(available)+j+1
                        if kind not in group:
                            ax = fig.add_subplot(len(rows), len(available), idx)
                            ax.text(.5,.5,'Not available',ha='center',va='center',transform=ax.transAxes)
                            ax.set_axis_off()
                        else:
                            data, wcs, header = self.read_image(group[kind])
                            ax = fig.add_subplot(len(rows), len(available), idx, projection=wcs)
                            vals = data[np.isfinite(data)]
                            lo, hi = np.percentile(vals, (1,99)) if vals.size else (0,1)
                            if hi <= lo:
                                hi = lo + max(abs(float(lo))*1e-6, 1e-12)
                            im = ax.imshow(data, origin='lower', cmap='viridis', vmin=lo, vmax=hi)
                            x,y = skycoord_to_pixel(coord, wcs, origin=0)
                            if np.isfinite(x) and np.isfinite(y) and -.5 <= x < data.shape[1]-.5 and -.5 <= y < data.shape[0]-.5:
                                ax.plot(x,y,'+',color='white',ms=10)
                            ax.set_xlabel('RA'); ax.set_ylabel('Dec')
                            ax.coords.grid(color='white',ls=':',alpha=.4)
                            fig.colorbar(im, ax=ax, fraction=.046, pad=.04, label=header.get('BUNIT',''))
                        ax.set_title(f'{band} μm {kind}')
                fig.suptitle(f'{target} — {prefix} ({title_suffix})')
                fig.tight_layout()
                fig.savefig(out/f'{prefix}{filename_suffix}_{self.clean_name(target)}.png',dpi=150,bbox_inches='tight')
            finally:
                plt.close(fig)
