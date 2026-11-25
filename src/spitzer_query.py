"""NASA/IPAC IRSA SEIP archive queries and FITS downloads."""
from io import BytesIO
from pathlib import Path
from urllib.parse import urlsplit, unquote
import re
import time
import requests
from astropy.table import Table
from astropy.io import fits, votable

SIA_URL = 'https://irsa.ipac.caltech.edu/SIA'
BANDS = {'IRAC.1': '3.6', 'IRAC.2': '4.5', 'IRAC.3': '5.8', 'IRAC.4': '8.0', 'MIPS.1': '24'}
PRODUCTS = {'mosaic': 'science', 'unc': 'uncertainty', 'cov': 'coverage', 'std': 'std', 'extmask': 'mask', 'median_mosaic': 'median_science', 'median_mosaic_unc': 'median_uncertainty'}


def classify(url):
    name = unquote(Path(urlsplit(url).path).name)
    match = re.fullmatch(r'(.+)\.(IRAC\.[1-4]|MIPS\.1)\.([a-z_]+)\.fits(?:\.gz)?', name, re.I)
    if not match:
        return None
    prefix, channel, suffix = match.groups()
    kind = PRODUCTS.get(suffix.lower())
    if kind is None:
        return None
    return BANDS[channel.upper()], kind, prefix, name


class SpitzerQuery:
    def __init__(self, timeout=(30, 180)):
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers['User-Agent'] = 'spitzer-downloads/0.1.0'

    def search(self, ra, dec, radius_arcmin=6):
        if not 0 <= float(ra) < 360 or not -90 <= float(dec) <= 90:
            raise ValueError('RA and Dec must be in valid ICRS degrees')
        if not 0 < float(radius_arcmin) <= 180:
            raise ValueError('radius_arcmin must be in (0, 180]')
        params = {'COLLECTION': 'spitzer_seip', 'POS': f'CIRCLE {ra} {dec} {radius_arcmin/60}', 'RESPONSEFORMAT': 'VOTABLE'}
        response = self.session.get(SIA_URL, params=params, timeout=self.timeout)
        response.raise_for_status()
        doc = votable.parse(BytesIO(response.content))
        for resource in doc.resources:
            for info in resource.infos:
                if info.name == 'QUERY_STATUS' and info.value != 'OK':
                    raise RuntimeError(f'IRSA query status: {info.value} {info.content}')
        table = doc.get_first_table().to_table()
        if 'access_url' not in table.colnames:
            raise RuntimeError('IRSA response has no access_url column')
        return table

    @staticmethod
    def valid_fits(path):
        try:
            with fits.open(path, memmap=False) as hdul:
                hdul.verify('exception')
                return any(h.data is not None and h.data.ndim == 2 for h in hdul)
        except Exception:
            return False

    def download(self, url, path, overwrite=False):
        path = Path(path)
        if path.exists() and not overwrite and self.valid_fits(path):
            return path
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name('part-' + path.name)
        for attempt in range(3):
            try:
                with self.session.get(url, stream=True, timeout=self.timeout) as response:
                    response.raise_for_status()
                    with tmp.open('wb') as f:
                        for chunk in response.iter_content(1024 * 1024):
                            if chunk:
                                f.write(chunk)
                if not self.valid_fits(tmp):
                    raise ValueError(f'Invalid FITS: {url}')
                tmp.replace(path)
                return path
            except (requests.RequestException, OSError, ValueError):
                tmp.unlink(missing_ok=True)
                if attempt == 2:
                    raise
                time.sleep(2 ** attempt)
