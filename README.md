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




For running the software you need a tab separated file with header line, kmers on first column and annotation and enrichment number values on other columns. Example file is at `tests/enriched_kmers_z99_HOXB13.tsv`. You should get a session:

```
$ huddinge_tsne_browser -i tests/enriched_kmers_z99_HNF4A.tsv -t  enriched_kmers_HNF4A 2>/dev/null 
<script src="enriched_kmers_HNF4A.js" id="bd73b7e5-bd13-42c5-8a5c-c0bd4e3a1914"></script>
$ ls enriched_kmers_HNF4A*
enriched_kmers_HNF4A.html  enriched_kmers_HNF4A.js
$
```

The file `enriched_kmers_HNF4A.html` contains a sample HTML file you can open in a web browser.  To embed the plot in your own file, copy-paste the script tag from standard output where you want to have the plot.






Installation
============

The system has bit of dependencies with holoviews \>=1.9.2, jellyfish and jupyter
covering most of them. Those are easiest to install from conda-forge and bioconda.

If you have [bioconda](https://bioconda.github.io/) installed (you should), then it's enough to say

	conda install -c kpalin huddinge_tsne_browser pyhuddinge
