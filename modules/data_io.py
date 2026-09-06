"""Repository paths and explicit data-file selection.

Synthetic data must never silently replace the proprietary dataset used
for manuscript estimates.
"""

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = PROJECT_ROOT / "data"
FIGURE_DIR = PROJECT_ROOT / "figure"
PDATA_DIR = PROJECT_ROOT / "pdata"

REAL_DATA_PATH = DATA_DIR / "M_1920_2023.csv"
SYNTHETIC_DATA_PATH = DATA_DIR / "M_1920_2023_synthetic.csv"


def get_real_data_path() -> Path:
    """Return the proprietary-data path or fail explicitly."""
    if not REAL_DATA_PATH.exists():
        raise FileNotFoundError(
            "The proprietary manuscript dataset is not available. "
            "Synthetic data are not substituted automatically."
        )
    return REAL_DATA_PATH


def get_synthetic_data_path() -> Path:
    """Return the public synthetic-data path."""
    if not SYNTHETIC_DATA_PATH.exists():
        raise FileNotFoundError(
            f"Synthetic dataset not found: {SYNTHETIC_DATA_PATH}"
        )
    return SYNTHETIC_DATA_PATH
