from huddinge_tsne_browser import util


def test_memory_usage():
    assert util.memory_usage() > 0


def test_memory_usage_resource():
    assert util.memory_usage_resource() > 0


def test_memory_usage_psutil():
    assert util.memory_usage_psutil() > 0
