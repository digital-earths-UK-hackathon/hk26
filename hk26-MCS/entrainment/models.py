"""
Central configuration for models, analysis regions, S3 URLs, and filesystem paths.

All pipeline scripts import from this module so that model/region choices propagate
consistently without duplication.

Pipeline usage
--------------
    import models
    models.add_model_arg(parser)
    models.add_region_arg(parser)
    ...
    zarr = models.data_dir(args.model) / f'entrainment_{args.region}.zarr'
"""
from pathlib import Path

S3_BASE     = 'https://hackathon-o.s3-ext.jc.rl.ac.uk/sim-data/analysis/PyFLEXTRKR'
CATALOG_URL = 'https://digital-earths-global-hackathon.github.io/catalog/catalog.yaml'
DATE_RANGE  = '20200201.0000_20210301.0000'

# Four global UM models with PyFLEXTRKR tracking data at HEALPix zoom=9.
MODELS = {
    'um_glm_n2560_RAL3p3_tuned_hk26': {'zoom': 9, 'display': 'UM N2560 RAL3.3'},
    'um_glm_n2560_CoMA9_hk26':         {'zoom': 9, 'display': 'UM N2560 CoMA9'},
    'um_glm_n1280_GAL9_v2_hk26':       {'zoom': 9, 'display': 'UM N1280 GAL9'},
    'um_glm_n1280_CoMA9_hk26':         {'zoom': 9, 'display': 'UM N1280 CoMA9'},
}

# Analysis regions.  lon_min/lon_max are in [0, 360]; if lon_min > lon_max the
# region wraps across 0°.  buf_* bounds are in [-180, 180] and include a buffer
# used for pre-filtering MCS track centroids.
REGIONS = {
    'wam': {
        'display':     'WAM',
        'lon_min':     340, 'lon_max': 20,    # wraps: 340–360° + 0–20°
        'lat_min':     2,   'lat_max': 15,
        'buf_lon_min': -25, 'buf_lon_max': 25,
        'buf_lat_min': -3,  'buf_lat_max': 20,
    },
}


# ---------------------------------------------------------------------------
# URL constructors
# ---------------------------------------------------------------------------

def _pyflex_key(model):
    """PyFLEXTRKR S3 directory name, e.g. um_glm_n2560_RAL3p3_tuned_z9."""
    return model.replace('_hk26', f'_z{MODELS[model]["zoom"]}')


def mask_url(model):
    zoom = MODELS[model]['zoom']
    return (f'{S3_BASE}/{_pyflex_key(model)}/mcstracking/'
            f'mcs_mask_hp{zoom}_{DATE_RANGE}.zarr')


def stats_url(model):
    return (f'{S3_BASE}/{_pyflex_key(model)}/stats/'
            f'mcs_tracks_final_{DATE_RANGE}.nc')


# ---------------------------------------------------------------------------
# Filesystem path helpers  (all relative to the entrainment/ working dir)
# ---------------------------------------------------------------------------

def data_dir(model):
    return Path('/work/scratch-nopw2/mmuetz/hk26/hk26-MCS/entrainment/data') / model


def figs_dir(model):
    return Path('figs') / model


def done_dir(model):
    return Path('donefiles') / model


def init_donefile(model, region):
    return done_dir(model) / f'init_{region}.done'


def chunk_donefile(model, chunk_idx):
    return done_dir(model) / f'chunk_{chunk_idx:03d}.done'


# ---------------------------------------------------------------------------
# Argparse helpers
# ---------------------------------------------------------------------------

def add_model_arg(parser):
    parser.add_argument('--model', required=True, choices=list(MODELS),
                        metavar='MODEL', help='Model key (see models.MODELS)')


def add_region_arg(parser, default='wam'):
    parser.add_argument('--region', default=default, choices=list(REGIONS),
                        help='Analysis region (default: %(default)s)')
