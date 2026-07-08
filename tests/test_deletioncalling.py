import csv

import pytest
from click.testing import CliRunner

from mgatk.deletioncalling import clidel, clifind
from conftest import requires_dplyr, requires_snakemake


@pytest.fixture
def runner():
    return CliRunner()


def test_del_find_produces_clip_and_sa_tables(runner, tmp_path, pearsonbam_dir):
    bam = pearsonbam_dir / "CACCACTAGGAGGCGA-1.qc.bam"
    out_prefix = tmp_path / "out"

    result = runner.invoke(
        clifind.main,
        ["-i", str(bam), "-mc", "chrM", "-o", str(out_prefix), "-ml", "16569"],
    )
    assert result.exit_code == 0, result.output

    clip_file = tmp_path / "out.clip.tsv"
    sa_file = tmp_path / "out.SA.tsv"
    assert clip_file.exists()
    assert sa_file.exists()

    with open(clip_file) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    assert rows[0] == ["position", "coverage", "clip_count", "SA"]
    # one row per mitochondrial base, plus the header
    assert len(rows) == 16569 + 1

    with open(sa_file) as fh:
        sa_rows = list(csv.reader(fh, delimiter="\t"))
    assert sa_rows[0] == ["out1", "out2"]


@requires_snakemake
@requires_dplyr
def test_del_end_to_end_on_pearsonbam(runner, tmp_path, pearsonbam_dir):
    out_dir = tmp_path / "out"
    result = runner.invoke(
        clidel.main,
        [
            "-i", str(pearsonbam_dir), "-o", str(out_dir), "-z",
            "-lc", "6073,5000", "-rc", "13095,5000", "--snake-stdout",
        ],
    )
    assert result.exit_code == 0, result.output

    heteroplasmy_file = out_dir / "final" / "mgatk_del.deletion_heteroplasmy.tsv"
    assert heteroplasmy_file.exists()
    assert heteroplasmy_file.stat().st_size > 0
