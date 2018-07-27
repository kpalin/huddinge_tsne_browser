huddinge\_tsne\_browser
=======================

[![image](https://img.shields.io/travis/kpalin/huddinge_tsne_browser.svg)](https://travis-ci.org/kpalin/huddinge_tsne_browser)

Tool to browse kmer sequences laid out in 2D approximating Huddinge
distance and displaying selex enrichment.  Current visualisation shows
10 local maxima laid out in polar coordinates such that the distance
from origin is the mean log fold enrichment in the selex
experiments. The other kmers that are enriched with z-score>2.58 are
displayed in the (jittered) direction of their representative local
maxima. The representative is the most enriched reached with steps of
huddinge size one.  The points are colored according to their distance
to the representative. The circle is at the highest enrichment peak.




For running the software you need a tab separated file with header line, kmers on first column and annotation and enrichment number values on other columns. Example file is at `examples/enriched_kmers_z99_HOXB13.tsv`




To use the package efficiently, you should download, install and run the
jupyter notebook examples/Huddinge TSNE browser 2.ipynb from the github
repository. A static sample can be [browsed in nbviewer](http://nbviewer.jupyter.org/github/kpalin/huddinge_tsne_browser/blob/master/examples/Huddinge%20TSNE%20browser%202.ipynb)


Installation
============

The system has bit of dependencies with holoviews \>=1.9.2, jellyfish and jupyter
covering most of them. Those are easiest to install from conda-forge and bioconda.

If you have [bioconda](https://bioconda.github.io/) installed (you should), then it's enough to say

	conda install -c kpalin huddinge_tsne_browser pyhuddinge
