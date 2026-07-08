import pandas as pd

from mgatk.processing.variant_calling import _gather_possible_variants


def _make_base_df(mito_length, cell_name, nonzero):
    """One-cell coverage DataFrame; `nonzero` maps 1-based position -> count."""
    row = {pos: nonzero.get(pos, 0) for pos in range(1, mito_length + 1)}
    return pd.DataFrame([row], index=[cell_name])


def test_gather_possible_variants_breaks_position_ties_deterministically(tmp_path):
    mito_length = 5
    cell = "cell1"

    # Two distinct alt alleles (C and G) both observed at position 3 on both
    # strands: a genuine tie on position that the sort must still resolve
    # the same way regardless of input order.
    base_coverage_dict = {
        "A": (_make_base_df(mito_length, cell, {}), _make_base_df(mito_length, cell, {})),
        "T": (_make_base_df(mito_length, cell, {}), _make_base_df(mito_length, cell, {})),
        "C": (_make_base_df(mito_length, cell, {3: 5}), _make_base_df(mito_length, cell, {3: 5})),
        "G": (_make_base_df(mito_length, cell, {3: 7}), _make_base_df(mito_length, cell, {3: 7})),
    }

    reference_file = tmp_path / "ref.txt"
    reference_file.write_text("\n".join(f"{pos}\tA" for pos in range(1, mito_length + 1)) + "\n")

    variants = _gather_possible_variants(base_coverage_dict, str(reference_file), mito_length)

    assert (3, "A", "C") in variants
    assert (3, "A", "G") in variants
    # The whole list must already be in (position, observed base) order: ties on
    # position alone (e.g. 3C vs 3G) must not depend on set-iteration order, which
    # is hash-randomized per Python process and therefore not reproducible.
    assert variants == sorted(variants, key=lambda v: (v[0], v[2]))
