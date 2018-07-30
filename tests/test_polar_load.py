from hypothesis import given,assume,reproduce_failure
from hypothesis.strategies import text,floats,integers
import pytest
from hypothesis.extra.pandas import column, data_frames, indexes
import numpy as np


ALPHABET=list("ACGT")




tsv_strategy =integers(min_value=2,max_value=15).flatmap(lambda n:data_frames(columns=[column("mean_ln_fold",elements=floats()),
                                    column("some_other_value",elements=floats())],
                        index=indexes(text(ALPHABET,min_size=n,max_size=n),min_size=2,max_size=200)
                    ))


tsv_strategy_finite =integers(min_value=2,max_value=15).flatmap(lambda n:data_frames(columns=[
                                    column("mean_ln_fold",elements=floats(allow_nan=False,allow_infinity=False,
                                        min_value=np.finfo(float).min/2,
                                        max_value=np.finfo(float).max/2)),
                                    column("some_other_value",elements=floats())],
                        index=indexes(text(ALPHABET,min_size=n,max_size=n),min_size=2,max_size=200)
                    ))


#@pytest.mark.skip(reason="Fails for many things")
@given(tsv_strategy_finite)
def test_fileread(tmpdir,dataf):
    from huddinge_tsne_browser.polarmapper import PolarMapper
    import os.path
    fname=os.path.join(tmpdir,"sample.tsv")
    dataf.to_csv(os.path.join(tmpdir,"sample.tsv"),sep="\t",header=True,index=True)
    
    d = PolarMapper(fname)
    d.tsne_obj[d.sole_binder].matrix.head()


@given(tsv_strategy)
def test_fileread_fail(tmpdir,dataf):
    assume(dataf.mean_ln_fold.isnull().any())

    from huddinge_tsne_browser.polarmapper import PolarMapper
    import os.path
    fname=os.path.join(tmpdir,"sample.tsv")
    dataf.to_csv(os.path.join(tmpdir,"sample.tsv"),sep="\t",header=True,index=True)
    
    with pytest.raises(ValueError):
        _ = PolarMapper(fname)
