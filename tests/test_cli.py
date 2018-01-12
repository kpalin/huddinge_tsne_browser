import pytest

from huddinge_tsne_browser import cli


def test_cli():
    assert cli.cli(["--input", "TFAP2A-head-10000.dist"]) is not None


def test_cli_no_input():
    with pytest.raises(SystemExit):
        assert cli.cli([]) is None
