import shutil

import pysam
import pytest

from mgatk import mgatkHelp as helpmod


# ---------------------------------------------------------------------------
# Pure string / list helpers
# ---------------------------------------------------------------------------

def test_string_hamming_distance_counts_mismatches():
    assert helpmod.string_hamming_distance("karolin", "kathrin") == 3
    assert helpmod.string_hamming_distance("AAAA", "AAAA") == 0


def test_string_hamming_distance_empty_strings():
    assert helpmod.string_hamming_distance("", "") == 0


def test_rev_comp_reverses_and_complements():
    assert helpmod.rev_comp("ACGT") == "ACGT"
    assert helpmod.rev_comp("AATTCCGG") == "CCGGAATT"
    assert helpmod.rev_comp("N") == "N"


def test_gettime_matches_expected_format():
    # e.g. "Sun Jul 05 21:28:24 EDT 2026: "
    stamp = helpmod.gettime()
    parts = stamp.strip(": ").split(" ")
    assert len(parts) == 6
    assert stamp.endswith(": ")


def test_findIdx_returns_matching_positions():
    assert helpmod.findIdx(["a", "b", "c", "b"], ["b"]) == [1, 3]
    assert helpmod.findIdx(["a", "b"], ["z"]) == []


# ---------------------------------------------------------------------------
# Filesystem helpers
# ---------------------------------------------------------------------------

def test_make_folder_creates_and_is_idempotent(tmp_path):
    target = tmp_path / "nested" / "dir"
    helpmod.make_folder(str(target))
    assert target.is_dir()
    # Calling again on an already-existing folder must not raise
    helpmod.make_folder(str(target))
    assert target.is_dir()


def test_file_len_counts_lines(tmp_path):
    f = tmp_path / "lines.txt"
    f.write_text("one\ntwo\nthree\n")
    assert helpmod.file_len(str(f)) == 3


def test_file_len_empty_file_is_zero(tmp_path):
    f = tmp_path / "empty.txt"
    f.write_text("")
    assert helpmod.file_len(str(f)) == 0


def test_split_barcodes_file_returns_original_when_under_threshold(tmp_path):
    barcode_file = tmp_path / "barcodes.txt"
    barcode_file.write_text("\n".join(f"BC{i}" for i in range(5)) + "\n")

    result = helpmod.split_barcodes_file(str(barcode_file), nsamples=10, output=str(tmp_path))
    assert result == [str(barcode_file)]


def test_split_barcodes_file_returns_original_when_nsamples_zero(tmp_path):
    barcode_file = tmp_path / "barcodes.txt"
    barcode_file.write_text("\n".join(f"BC{i}" for i in range(5)) + "\n")

    result = helpmod.split_barcodes_file(str(barcode_file), nsamples=0, output=str(tmp_path))
    assert result == [str(barcode_file)]


def test_split_barcodes_file_splits_when_over_threshold(tmp_path):
    barcode_file = tmp_path / "barcodes.txt"
    barcodes = [f"BC{i}" for i in range(10)]
    barcode_file.write_text("\n".join(barcodes) + "\n")

    result = helpmod.split_barcodes_file(str(barcode_file), nsamples=4, output=str(tmp_path))

    assert len(result) == 3  # ceil(10 / 4)
    all_lines = []
    for chunk in result:
        with open(chunk) as fh:
            lines = [line.strip() for line in fh if line.strip()]
        assert len(lines) <= 4
        all_lines.extend(lines)
    assert sorted(all_lines) == sorted(barcodes)


def test_parse_fasta_reads_single_contig(rcrs_fasta):
    sequences = helpmod.parse_fasta(str(rcrs_fasta))
    assert list(sequences.keys()) == ["chrM"]
    assert len(sequences["chrM"]) == 16569
    assert set(sequences["chrM"].upper()) <= set("ACGTN")


# ---------------------------------------------------------------------------
# Software / package availability checks
# ---------------------------------------------------------------------------

def test_check_software_exists_passes_for_real_binary():
    real_tool = shutil.which("python3") or shutil.which("ls")
    assert real_tool is not None
    helpmod.check_software_exists(real_tool.rsplit("/", 1)[-1])


def test_check_software_exists_exits_for_missing_binary():
    with pytest.raises(SystemExit):
        helpmod.check_software_exists("definitely_not_a_real_tool_xyz")


# ---------------------------------------------------------------------------
# BAM-aware helpers
# ---------------------------------------------------------------------------

def test_verify_bai_creates_missing_index(tmp_path, humanbam_dir):
    src_bam = humanbam_dir / "MGH60-P6-A11.mito.bam"
    dest_bam = tmp_path / "sample.bam"
    shutil.copyfile(src_bam, dest_bam)
    assert not (tmp_path / "sample.bam.bai").exists()

    helpmod.verify_bai(str(dest_bam))

    assert (tmp_path / "sample.bam.bai").exists()


def test_verify_sample_mitobam_true_for_matching_length(humanbam_dir):
    bam = humanbam_dir / "MGH60-P6-A11.mito.bam"
    assert helpmod.verify_sample_mitobam(str(bam), "chrM", 16571) is True


def test_verify_sample_mitobam_false_for_mismatched_length(humanbam_dir):
    bam = humanbam_dir / "MGH60-P6-A11.mito.bam"
    assert helpmod.verify_sample_mitobam(str(bam), "chrM", 16569) is False


def test_verify_sample_mitobam_false_for_missing_chromosome(humanbam_dir):
    bam = humanbam_dir / "MGH60-P6-A11.mito.bam"
    assert helpmod.verify_sample_mitobam(str(bam), "chrDOESNOTEXIST", 16571) is False


# ---------------------------------------------------------------------------
# Reference genome inference
# ---------------------------------------------------------------------------

@pytest.fixture
def supported_genomes(repo_root):
    fasta_dir = repo_root / "mgatk" / "bin" / "anno" / "fasta"
    return [f.stem for f in fasta_dir.glob("*.fasta")]


def test_handle_fasta_inference_known_genome_no_write(repo_root, supported_genomes, tmp_path):
    fastaf, mito_chr, mito_length = helpmod.handle_fasta_inference(
        "rCRS", supported_genomes, str(repo_root / "mgatk"), "call", str(tmp_path),
        write_files=False,
    )
    assert mito_chr == "chrM"
    assert mito_length == 16569
    assert fastaf.endswith("rCRS.fasta")
    # write_files=False must not create any output folders
    assert not (tmp_path / "fasta").exists()


def test_handle_fasta_inference_writes_expected_files(repo_root, supported_genomes, tmp_path):
    fastaf, mito_chr, mito_length = helpmod.handle_fasta_inference(
        "rCRS", supported_genomes, str(repo_root / "mgatk"), "call", str(tmp_path),
        write_files=True,
    )
    assert mito_length == 16569
    assert (tmp_path / "fasta" / "chrM.fasta").exists()
    assert (tmp_path / "final" / "chrM_refAllele.txt").exists()

    ref_lines = (tmp_path / "final" / "chrM_refAllele.txt").read_text().splitlines()
    assert len(ref_lines) == 16569
    assert ref_lines[0].split("\t")[0] == "1"


def test_handle_fasta_inference_unknown_genome_exits(repo_root, supported_genomes, tmp_path):
    with pytest.raises(SystemExit):
        helpmod.handle_fasta_inference(
            "not_a_real_genome", supported_genomes, str(repo_root / "mgatk"), "call",
            str(tmp_path), write_files=False,
        )


# ---------------------------------------------------------------------------
# Misc
# ---------------------------------------------------------------------------

def test_available_cpu_count_returns_positive_int():
    count = helpmod.available_cpu_count()
    assert isinstance(count, int)
    assert count >= 1
