# Spitzer Downloads

Download original Spitzer Enhanced Imaging Products (SEIP) mosaics from NASA/IPAC IRSA, together with their available uncertainty and coverage FITS files. Data acquisition and visual diagnostics only; no photometry, reprojection, PSF matching or scientific analysis.

## Install

```bash
pip install --upgrade git+https://github.com/akhilkrishnar0/spitzer-downloads.git
```

## Example: NGC1566

```python
from spitzer_pipeline_class import SpitzerPipeline

pipe = SpitzerPipeline(galaxy_centers={'NGC1566': (65.0042, -54.9378)})
pipe.run(['NGC1566'], bands=['3.6','4.5','5.8','8.0','24'],
         radius_arcmin=6, skip_download=True, do_plot=True)
```

The radius is an archive search radius, not a requested cutout size. All matching archive mosaics are preserved with native WCS and units. The archive search can include mosaics that intersect the search circle but do not contain the target centre. The inventory records whether the centre lies inside each FITS grid; this is not a test for valid exposure or full galaxy coverage. Multiple mosaics are not merged. Missing products are reported, not manufactured. `skip_download=True` reuses validated FITS files; `force_download=True` refreshes them. The search is repeated to discover the current archive inventory.

## Archive query

The SIA VOTable results are read using FIELD names rather than XML IDs. The parser identifies the data-access URL from the results resource and reports archive errors and returned columns if the field is missing.

## Products

Default: `science` (`mosaic.fits`), `uncertainty` (`unc.fits`), `coverage` (`cov.fits`). Optional products include `std`, `mask`, `median_science`, `median_uncertainty`; availability is archive-dependent. Uncertainty and coverage retain their archive definitions, not converted to inverse variance.

## Outputs

```
spitzer_data/NGC1566/
    archive_query.ecsv
    inventory.csv
    downloads/<band>um/<mosaic-prefix>/<original-filename>
    plots/<mosaic-prefix>.png
```

Each diagnostic shows available products by band for the same mosaic prefix, with a target-centre marker when inside the grid. Native shapes may differ. Full mosaics can be large. No actual NGC1566 availability is assumed.

## References

- [IRSA SEIP image search](https://irsa.ipac.caltech.edu/data/SPITZER/Enhanced/SEIP/)
- [IRSA SEIP SIA v2 example](https://irsa.ipac.caltech.edu/docs/program_interface/siav2_seip.html)
- [IRSA SIA v2 API](https://irsa.ipac.caltech.edu/ibe/sia.html)
- [SEIP overview](https://irsa.ipac.caltech.edu/data/SPITZER/Enhanced/SEIP/overview.html)
