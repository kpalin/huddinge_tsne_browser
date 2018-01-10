import pytest


@pytest.fixture
def tsne_obj():
    from huddinge_tsne_browser.tsne_mapper import TsneMapper
    import os.path
    dir = os.path.dirname(__file__)
    d = TsneMapper(os.path.join(dir, "TFAP2A-head-1000.dist"))
    return d


def test_tsne_read(tsne_obj):
    assert tsne_obj is not None


def test_matrix_formatting(tsne_obj):
    import numpy as np
    assert np.diag(tsne_obj.matrix).sum() == 0
    assert np.allclose(tsne_obj.matrix, tsne_obj.matrix.T)
