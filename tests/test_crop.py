import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from spitzer_pipeline_class import SpitzerPipeline


def test_native_crop(tmp_path):
    pipe = SpitzerPipeline(base_dir=tmp_path, galaxy_centers={"test": (65., -55.)})
    prefix = "example"
    folder = tmp_path/"test"/"downloads"/"3.6um"/prefix
    folder.mkdir(parents=True)
    w = WCS(naxis=2)
    w.wcs.crpix = [101., 101.]
    w.wcs.crval = [65., -55.]
    w.wcs.cdelt = [-1./3600, 1./3600]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    h = w.to_header()
    h["BUNIT"] = "MJy/sr"
    image = np.arange(201*201, dtype=float).reshape(201,201)
    fits.writeto(folder/"example.IRAC.1.mosaic.fits", image, h)
    fits.writeto(folder/"example.IRAC.1.unc.fits", image+1, h)
    result = pipe.crop("test", size_arcmin=3, do_plot=False)
    pair = result[("3.6", prefix)]
    sci, sh = fits.getdata(pair["science"], header=True)
    unc = fits.getdata(pair["uncertainty"])
    assert sci.shape == (180,180)
    assert np.array_equal(unc, sci+1)
    assert sh["BUNIT"] == "MJy/sr"
    x,y = WCS(sh).world_to_pixel_values(65., -55.)
    assert abs(x-90)<1 and abs(y-90)<1
    assert fits.getdata(folder/"example.IRAC.1.mosaic.fits").shape == (201,201)
