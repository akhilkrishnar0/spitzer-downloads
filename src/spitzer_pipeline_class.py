"""Download original SEIP mosaics and display diagnostic figures."""
from pathlib import Path
import csv
import re
import warnings
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
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
        table.write(root/'archive_query.ecsv', format='ascii.ecsv', overwrite=True)
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
        with (root/'inventory.csv').open('w', newline='') as f:
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

    def plot(self, target, coord, groups, products, radius_arcmin):
        out = self.base_dir/self.clean_name(target)/'plots'
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
                fig.suptitle(f'{target} — {prefix} (native mosaics)')
                fig.tight_layout()
                fig.savefig(out/f'{prefix}.png',dpi=150,bbox_inches='tight')
            finally:
                plt.close(fig)
