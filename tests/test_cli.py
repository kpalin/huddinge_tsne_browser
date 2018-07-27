import pytest
from hypothesis import given
import hypothesis.strategies as st
from huddinge_tsne_browser import cli


def test_cli():
    assert cli.cli(["--input", "TFAP2A-head-10000.dist"]) is not None


def test_cli_no_input():
    with pytest.raises(SystemExit):
        assert cli.cli([]) is None

@given(st.sampled_from(["HOXB13","HNF4A"]),st.sampled_from(["mean_ln_fold","fold_z"]))
def test_main_make_html(binder,enrich_column):
    import os
    import os.path
    from huddinge_tsne_browser.main import main
    fout="examples/test_enriched_kmers_z99_{}_{}".format(binder,enrich_column)
    try:
        os.remove(fout+".html")
    except FileNotFoundError:
        pass

    try:
        os.remove(fout+".js")
    except FileNotFoundError:
        pass

    args = ["-i","examples/enriched_kmers_z99_{}.tsv".format(binder),
        "--html",fout,
        "--enrichment_column",enrich_column,
        "--verbose"]
    print(args)
    main(args)


    assert os.path.exists(fout+".html")
    assert os.path.exists(fout+".js")
