from importlib.resources import files as irfiles
from pathlib import Path
import numpy as np

DEFAULT_APP_DIR = Path.home() / ".local" / "share" / "preph"


def load_np_asset(name: str) -> np.ndarray:
    primitive_dir = irfiles().joinpath("resources")
    return np.load(str(primitive_dir.joinpath(name)))


def load_static_data():
    stacking_matrix = load_np_asset("stacking_matrix.npy")
    bulge_list = load_np_asset("bulge_list.npy")
    intl11_matrix = load_np_asset("intl11_matrix.npy")
    intl12_matrix = load_np_asset("intl12_matrix.npy")
    intl22_matrix = load_np_asset("intl22_matrix.npy")

    return stacking_matrix, bulge_list, intl11_matrix, intl12_matrix, intl22_matrix


def get_kmer_table_path(k: int, g: int, folder: Path = DEFAULT_APP_DIR, suffix: str = "mers_stacking_energy_binary.npy"):
    return folder / (str(k) + str(g) + suffix)


def load_kmer_table(k: int, g: int, folder: Path = DEFAULT_APP_DIR) -> np.ndarray:
    table_path = get_kmer_table_path(k, g, folder)
    if not table_path.exists():
        raise FileNotFoundError(f"k-mer stacking energy table not found at {table_path}. Build it with preph-precalculate-energies.")
    return np.load(table_path)