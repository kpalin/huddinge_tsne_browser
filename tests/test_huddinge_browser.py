import pytest
from hypothesis import given
import hypothesis.strategies as st


@pytest.fixture
def tsne_obj():
    from huddinge_tsne_browser.tsne_mapper import TsneMapper
    import os.path
    dir = os.path.dirname(__file__)
    fname = os.path.join(dir, "../examples/enriched_kmers_z99_HOXB13.tsv")

    d = TsneMapper(os.path.join(dir, "TFAP2A-head-1000.dist"))
    return d


@pytest.fixture
def tsne_laidout():
    from huddinge_tsne_browser.tsne_mapper import TsneMapper
    import os.path
    dir = os.path.dirname(__file__)
    laidout_file = os.path.join(dir, "TFAP2A-head-1000.laidout")
    if os.path.exists(laidout_file):
        d = TsneMapper(laidout_file)
    else:
        d = TsneMapper(os.path.join(dir, "TFAP2A-head-1000.dist"))
        d.compute_mds()
        #d.write_data(laidout_file)

    return d


def test_tsne_read(tsne_obj):
    assert tsne_obj is not None


def test_matrix_formatting(tsne_obj):
    import numpy as np
    assert np.diag(tsne_obj.matrix).sum() == 0
    np.testing.assert_allclose(tsne_obj.matrix, tsne_obj.matrix.T)


def matrix_dtype(tsne_obj):
    assert tsne_obj.distance.dtype == int
    assert tsne_obj.matrix.shape == (tsne_obj.N, tsne_obj.N)
    assert tsne_obj.matrix.dtype == float


@pytest.mark.skip(reason="Slow test")
def test_layout_mds(tsne_obj):
    tsne_obj.compute_mds()
    assert tsne_obj.laidout()


@pytest.mark.skip(reason="Slow test")
def test_layout_tsne(tsne_obj):
    tsne_obj.compute_tsne()
    assert tsne_obj.laidout()


@pytest.mark.skip(reason="Slow test")
def test_layout_spectral(tsne_obj):
    tsne_obj.compute_spectral()
    assert tsne_obj.laidout()


def test_layout_clear(tsne_laidout):
    assert tsne_laidout.laidout()
    tsne_laidout.clear_layout()
    assert not tsne_laidout.laidout()


@given(seqs=st.lists(
    st.sampled_from(list(tsne_obj().sequences[0])),
    min_size=1,
    max_size=10,
    unique=True))
def test_subsetting_count(tsne_obj, seqs):
    from copy import deepcopy
    tsne_obj = deepcopy(tsne_obj)

    tsne_obj.subset_sequences(seqs)
    assert len(tsne_obj) == len(seqs)
    assert set(seqs) == set(tsne_obj.sequences[0])


@given(seqs=st.lists(
    st.sampled_from(list(tsne_obj().sequences[0])),
    min_size=1,
    max_size=10,
    unique=True))
def test_subsetting_distances(tsne_obj, seqs):
    from copy import deepcopy
    tsne_old = deepcopy(tsne_obj)

    import numpy as np
    tsne_new = deepcopy(tsne_old)

    tsne_new.subset_sequences(seqs)
    assert len(tsne_new) == len(seqs)

    seqs_order = list(tsne_new.sequences[0])
    old_i = tsne_old.sequences.set_index(0).index

    old_idx = [old_i.get_loc(x) for x in seqs_order]

    
    old_dist = tsne_old.matrix.iloc[old_idx, old_idx]
    new_dist = tsne_new.matrix
    np.testing.assert_allclose(old_dist, new_dist)


@pytest.mark.skip(reason="Slow test")
@given(seqs=st.lists(
    st.sampled_from(list(tsne_obj().sequences[0])),
    min_size=3,
    max_size=10,
    unique=True))
def test_layout_write_read(tsne_obj, tmpdir, seqs):
    from copy import deepcopy
    import os.path
    import numpy as np
    tsne_obj = deepcopy(tsne_obj)

    tsne_obj.subset_sequences(seqs)
    tsne_obj.compute_mds()

    tf = tmpdir.join("new_layout.mds")
    tsne_obj.write_data(str(tf))
    from huddinge_tsne_browser.tsne_mapper import TsneMapper
    loaded_obj = TsneMapper(str(tf))

    assert (loaded_obj.sequences[0].reset_index(
        drop=True) == tsne_obj.sequences[0].reset_index(drop=True)).all()

    np.testing.assert_allclose(
        loaded_obj.embedding, tsne_obj.embedding, atol=1e-7, rtol=1e-3)

    np.testing.assert_allclose(
        loaded_obj.fit_measure_, tsne_obj.fit_measure_, atol=1e-7, rtol=1e-3)
