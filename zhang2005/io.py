"""I/O helpers for MATLAB files produced by the original replication code."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import loadmat, savemat


def load_mat(path: str | Path) -> dict[str, Any]:
    """Load a MATLAB ``.mat`` file and drop scipy metadata keys."""
    data = loadmat(path, squeeze_me=True, struct_as_record=False)
    return {k: v for k, v in data.items() if not k.startswith("__")}


def save_npz(path: str | Path, **arrays: np.ndarray) -> None:
    """Save Python replication arrays in NumPy's compressed format."""
    np.savez_compressed(path, **arrays)


def save_mat(path: str | Path, **arrays: np.ndarray) -> None:
    """Save arrays for MATLAB interoperability."""
    savemat(path, arrays)
