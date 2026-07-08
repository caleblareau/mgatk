import shutil
import subprocess
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parent


def _r_packages_available(packages):
    if shutil.which("Rscript") is None:
        return False
    expr = "sapply(c({}), requireNamespace, quietly=TRUE)".format(
        ", ".join('"{}"'.format(p) for p in packages)
    )
    try:
        out = subprocess.run(
            ["Rscript", "-e", expr], capture_output=True, text=True, timeout=60
        )
    except (subprocess.TimeoutExpired, OSError):
        return False
    return out.returncode == 0 and "FALSE" not in out.stdout


requires_snakemake = pytest.mark.skipif(
    shutil.which("snakemake") is None, reason="snakemake is required to run the pipeline"
)

requires_mgatk_r_packages = pytest.mark.skipif(
    not _r_packages_available(
        ["data.table", "SummarizedExperiment", "GenomicRanges", "Matrix"]
    ),
    reason="R with data.table/SummarizedExperiment/GenomicRanges/Matrix is required",
)

requires_dplyr = pytest.mark.skipif(
    not _r_packages_available(["dplyr"]), reason="R with dplyr is required"
)


@pytest.fixture(scope="session")
def repo_root():
    return REPO_ROOT


@pytest.fixture(scope="session")
def humanbam_dir():
    return TESTS_DIR / "humanbam"


@pytest.fixture(scope="session")
def barcode_dir():
    return TESTS_DIR / "barcode"


@pytest.fixture(scope="session")
def pearsonbam_dir():
    return TESTS_DIR / "pearsonbam"


@pytest.fixture(scope="session")
def rcrs_fasta(repo_root):
    return repo_root / "mgatk" / "bin" / "anno" / "fasta" / "rCRS.fasta"
