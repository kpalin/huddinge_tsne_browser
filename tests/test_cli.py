from huddinge_tsne_browser import cli


def test_cli():
    assert cli.cli([]) is not None
