import logging as log


class TsneMapper(object):
    """Reader and tsne transformer for huddinge distance files
    """

    def __init__(self, input_file):
        """
        
        Arguments:
        - `input_file`:
        """
        self._input_file = input_file
        self.read_data()

    def read_data(
            self, ):
        """Read data from 'moder'
        """
        import pandas as pd
        from cStringIO import StringIO

        fin = open(self._input_file)
        l = fin.readline()
        bytes_read = len(l)
        p = l.split()
        self.N = int(p[0])
        if len(p) > 1:
            self.KLdivergence_ = float(p[1])

        l = [fin.readline() for _ in range(self.N)]
        bytes_read += sum(len(x) for x in l)

        self.sequences = pd.read_table(
            StringIO("".join(l)),
            nrows=self.N,
            sep="\t",
            header=None,
            engine="python")

        assert self.N == len(self.sequences)

        log.info("Read %d sequnces.", self.N)

        log.info("Reading distances.")
        self.read_distances(fin)

        if self.sequences.shape[1] > 1:
            log.info("Setting embedding from input data")
            self.embedding = self.sequences.set_index(0)
            self.embedding.columns = ["tsne0", "tsne1"]
            self.embedding.index.name = "Sequence"

    def read_distances(self, fin):
        """Read distance data from open file fin
        
        Arguments:
        - `fin`:
        """
        import pandas as pd
        import numpy as np

        d_pos = fin.tell()
        d = np.fromfile(fin, dtype="uint8", count=-1)
        assert len(d) == self.N * (self.N - 1) / 2
        log.info("Reshaping..")
        i, j = np.tril_indices(self.N, k=-1)
        self.distances = pd.DataFrame({0: i, 1: j, 2: d})

        #self.distances = pd.read_table(
        #    fin,
        #    sep="\t",
        #    header=None,
        #    engine="python",
        #    dtype={0: int,
        #           1: int,
        #           2: float})

        assert len(self.distances) == self.N * (self.N - 1) / 2

    def laidout(
            self, ):
        """Are the sequences laid out properly
        """
        return self.sequences.shape[1] > 1

    def _get_matrix(self):
        if not hasattr(self, "_matrix"):
            import numpy as np
            tot = self.distances[2].sum()

            M = self.distances.set_index([0, 1])[2].unstack().fillna(0)
            M[M.index[-1]] = 0.0
            M.loc[M.columns[0]] = 0.0

            self._matrix = M + M.T
            assert np.abs((tot * 2) - self._matrix.sum().sum()) < 1e-5, (
                2 * tot - self._matrix.sum().sum())

        return self._matrix

    matrix = property(_get_matrix)

    def compute_tsne(self, perplexity=30.0, fake=False):
        """
        """
        import pandas as pd
        if fake:
            import numpy as np
            self.embedding = np.random.normal(size=(len(self.sequences), 2))
        else:
            from sklearn.manifold import TSNE
            self.seq_tsne = TSNE(
                perplexity=perplexity,
                metric="precomputed",
                verbose=9,
                init="random")
            self.embedding = self.seq_tsne.fit_transform(self.matrix)

        self.embedding = pd.DataFrame(
            self.embedding,
            columns=["tsne0", "tsne1"],
            index=self.sequences[0])

        self.embedding.index.name = "Sequence"
        self.KLdivergence_ = self.seq_tsne.kl_divergence_

    def write_data(self, outfile):
        """
        
        Arguments:
        - `outfile`:
        """

        log.info("Writing %s" % (outfile))
        with open(outfile, "w") as outf:
            outf.write("%d\t%g\n" % (len(self.sequences), self.KLdivergence_))
            self.embedding.to_csv(outf, sep="\t", header=False, index=True)
            self.distances[2].copy().astype("uint8").values.tofile(outf)
            #self.distances.to_csv(outf, sep="\t", header=False, index=False)

    def holoview_plot(
            self, ):
        """
        """
        import holoviews as hv
        from bokeh.models import HoverTool
        hover = HoverTool(tooltips=[("index", "$index"), ("seq", "@Sequence")])

        p = hv.Scatter(
            self.embedding.reset_index(),
            kdims=["tsne0", "tsne1"],
            vdims=["Sequence"],
            label="Sequences",
        ).opts(plot=dict(tools=[hover, "box_select"])).opts(plot=dict(
            width=800, height=800))
        #p.pprint()

        return None
        #+hv.Histogram(self.distances[0], label="Distances", bins=20)

    def html(
            self, ):
        """
        """
        import holoviews as hv
        from bokeh.io import curdoc
        renderer = hv.renderer('bokeh')
        p = self.holoview_plot()
        #hvplot_p = renderer.get_plot(self.holoview_plot(), curdoc())
        doc = renderer.server_doc(p)
        doc.title = "TSNE browser"

        return doc
