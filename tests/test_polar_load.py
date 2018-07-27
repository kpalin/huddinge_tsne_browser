from hypothesis import given
from hypothesis.strategies import text,floats,integers
import pytest
from hypothesis.extra.pandas import column, data_frames, indexes

ALPHABET=list("ACGT")


tsv_strategy =data_frames(columns=[column("mean_ln_fold",elements=floats()),
                                    column("some_other_value",elements=floats())],
                        index=indexes(integers(min_value=2,max_value=15).flatmap(lambda n:text(ALPHABET,min_size=n,max_size=n)),
                                    min_size=10)
                    )

@pytest.mark.skip(reason="Fails for many things")
@given(tsv_strategy)
def test_fileread(tmpdir,dataf):
    from huddinge_tsne_browser.polarmapper import PolarMapper
    import os.path
    fname=os.path.join(tmpdir,"sample.tsv")
    dataf.to_csv(os.path.join(tmpdir,"sample.tsv"),sep="\t",header=True,index=True)
    if len(dataf)>1:
        d = PolarMapper(fname)

        d.save_bokeh_points(os.path.join(tmpdir,"silly"))
        assert os.path.exists(os.path.join(tmpdir,"silly.html"))