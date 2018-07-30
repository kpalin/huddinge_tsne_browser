import pytest
from hypothesis import given
import hypothesis.strategies as st


@pytest.fixture
def polarmap_obj():
    from huddinge_tsne_browser.polarmapper import PolarMapper
    import os.path
    dir = os.path.dirname(__file__)
    fname = os.path.join(dir, "../examples/enriched_kmers_z99_HOXB13.tsv")
    d = PolarMapper(fname)
    return d


def test_polarmap_init(polarmap_obj):
    print(str(polarmap_obj))

def test_get_local_maxima(polarmap_obj):
    n=10
    local_maximas,representatives  = polarmap_obj.get_local_maxima(polarmap_obj.sole_binder,n=n)
    print(local_maximas)
    print(representatives)
    print("Number unique: %d"%(representatives.nunique()))
    print("Number unique(index): %d"%(representatives.index.to_series().nunique()))
    print(representatives.value_counts())

    assert (representatives.isin(local_maximas) | (representatives.values==representatives.index.values)).all()
    vals = (polarmap_obj.selected_kmers.loc[polarmap_obj.sole_binder,polarmap_obj.enrichment_column])
    assert (vals.loc[representatives.index].values <= vals.loc[representatives].values).all()
    #assert len(local_maximas)<n






def test_plot_polar(polarmap_obj):
    import holoviews as hv
    import os.path
    plot = polarmap_obj.plot_polar()
    
    hv.renderer('bokeh').save(plot,"test_enrichment_kmers_HOXB13")

    assert os.path.exists("test_enrichment_kmers_HOXB13.html")



def test_plot_points(polarmap_obj):
    import holoviews as hv
    import os.path
    points, _,_ = polarmap_obj.plot_points("Default")
    
    hv.renderer('bokeh').save(points,"test_enrichment_kmers_HOXB13_points_bare")

    assert os.path.exists("test_enrichment_kmers_HOXB13_points_bare.html")


def test_bokeh_points(polarmap_obj):
    import os.path
    polarmap_obj.save_bokeh_points("test_enrichment_kmers_HOXB13_points")


    assert os.path.exists("test_enrichment_kmers_HOXB13_points.html")


