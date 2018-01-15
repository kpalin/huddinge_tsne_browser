huddinge\_tsne\_browser
=======================

[![image](https://img.shields.io/travis/kpalin/huddinge_tsne_browser.svg)](https://travis-ci.org/kpalin/huddinge_tsne_browser)

Tool to browse sequence kmers laid out in 2D with TSNE approximating
Huddinge distance.

To use the package efficiently, you should download, install and run the
jupyter notebook examples/Huddinge TSNE browser 2.ipynb from the github
repository. A static sample can be [browsed in nbviewer](http://nbviewer.jupyter.org/github/kpalin/huddinge_tsne_browser/blob/master/examples/Huddinge%20TSNE%20browser%202.ipynb)


Installation
============

The system has bit of dependencies with holoviews \>=1.9.2, jellyfish and jupyter
covering most of them. Those are easiest to install from conda-forge and bioconda.

If you have [bioconda](https://bioconda.github.io/) installed (you should), then it's enough to say

	conda install -c kpalin huddinge_tsne_browser
