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



def test_plot_polar(polarmap_obj):
    import holoviews as hv
    import os.path
    plot = polarmap_obj.plot_polar()
    
    hv.renderer('bokeh').save(plot,"test_enrichment_kmers_HOXB13")

    assert os.path.exists("test_enrichment_kmers_HOXB13.html")



def test_plot_points(polarmap_obj):
    import holoviews as hv
    import os.path
    points, enrichment_r = polarmap_obj.plot_points("Default")
    
    hv.renderer('bokeh').save(points,"test_enrichment_kmers_HOXB13_points")

    assert os.path.exists("test_enrichment_kmers_HOXB13_points.html")


def test_bokeh_points(polarmap_obj):
    import holoviews as hv
    import os.path
    polarmap_obj.save_bokeh_points("Default","test_enrichment_kmers_HOXB13_points")


    assert os.path.exists("test_enrichment_kmers_HOXB13_points.html")


