import logging as log

from . import util


def read_jf(fname):
    "Read jellyfish output file (kmer counts)"
    import pandas as pd
    import json
    fin = open(fname)
    b = int(fin.read(9))

    j = fin.read(b)
    j = json.loads(j.strip("\x00"))
    d = pd.read_table(fin, sep=" ", header=None, names=["Sequence", "Count"])
    return j, d.set_index("Sequence").Count


class TsneMapper(object):
    """Reader and tsne transformer for huddinge distance files
    """

    def __init__(self, input_file=None,force_distances=False):
        """
        
        Arguments:
        - `input_file`:
        """

        self._input_file = input_file
        self.data_dims = []
        self.coord_dims = ["tsne0", "tsne1"]

        self.read_data(force_distances)

    def __len__(
            self, ):
        """Return number of kmers
        """
        return len(self.sequences)

    def read_data(
            self, force_distances=False):
        """Read data from 'moder'
        """
        import pandas as pd
        from io import StringIO

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

        log.info("Read %d sequences.", self.N)

        if self.sequences.shape[1] > 1:
            log.info("Setting embedding from input data")
            self.embedding = self.sequences.set_index(0)
            self.embedding.columns = self.coord_dims
            self.embedding.index.name = "Sequence"
            
        if (not hasattr(self,"embedding")) or force_distances:
            log.info("Memory usage %gMB" % (util.memory_usage()))
            log.info("Reading distances.")
            self.read_distances(fin)
            log.info("Memory usage %gMB" % (util.memory_usage()))

    def set_kmer_values(self, annot, append=False):
        """Set kmer annotation values
        
        Arguments:
        - `annot`:
        - `append`:
        """
        new_emb = self.embedding
        if not append:
            new_emb = new_emb.drop(self.data_dims, axis=1)
            self.data_dims = []

        annot = annot.drop(
            [x for x in annot.columns if x in self.coord_dims], axis=1)
        new_emb = new_emb.join(annot)

        self.data_dims.extend(annot)
        self.embedding = new_emb

    def add_kmercounts(self, name, filename):
        """
        
        Arguments:
        - `name`:
        - `filename`:
        """
        _, counts = read_jf(filename)

        try:
            self.embedding[name] = counts.reindex(index=self.embedding.index,fill_value=0)
            if name not in self.data_dims:
                self.data_dims.append(name)
        except KeyError as e:
            kmer_lens = sorted(set(len(x) for x in counts.index))
            log.warning(
                "Coudn't add %s kmers. They are of length %s while the embedding is for kmers of length %d"
                % (name, ",".join(str(x) for x in kmer_lens), self.kmer_size))

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

        self.distances = d

        #self.distances = pd.DataFrame(M,index=np.arange(self.N),columns=np.arange(self.N))

        #self.distances = pd.read_table(
        #    fin,
        #    sep="\t",
        #    header=None,
        #    engine="python",
        #    dtype={0: int,
        #           1: int,
        #           2: float})

        #assert len(self.distances) == self.N * (self.N - 1) / 2

    def _get_kmer_size(self):
        return len(self.sequences[0])

    kmer_size = property(_get_kmer_size)

    def laidout(
            self, ):
        """Are the sequences laid out properly
        """
        return self.sequences.shape[1] > 1

    def _get_matrix(self):
        if not hasattr(self, "_matrix"):
            import numpy as np
            log.info("Memory usage before matrix formatting %gMB" %
                     (util.memory_usage()))

            log.info("Reshaping..")
            log.info("sizeof(distances) = %gMB" % (self.distances.nbytes /
                                                   (2.0**20)))
            M = np.zeros((self.N, self.N), dtype="float32")
            log.info("sizeof(M) = %gMB" % (M.nbytes / (2.0**20)))

            idx = np.tril_indices(self.N, k=-1)
            log.info("sizeof(idx) = %gMB" % (sum(x.nbytes
                                                 for x in idx) / (2.0**20)))

            M[idx] = self.distances

            M += M.T
            self._matrix = M

            # Due to output format of all_pairs_huddinge, the following does *not* work.
            #from scipy.spatial.distance import squareform
            #self._matrix = squareform(self.distances.astype("float32"))
            #assert np.allclose(M, self._matrix)
            log.info("Memory usage after matrix formatting %gMB" %
                     (util.memory_usage()))
            del (M)
            del (idx)
            log.info("Memory usage after cleaning %gMB" %
                     (util.memory_usage()))

        return self._matrix

    matrix = property(_get_matrix)

    def compute_tsne(self, perplexity=30.0, fake=False):
        """
        """
        import pandas as pd
        if fake:
            import numpy as np
            self.embedding = np.random.normal(size=(self.kmer_size, 2))
        else:
            from sklearn.manifold import TSNE
            self.seq_tsne = TSNE(
                perplexity=perplexity,
                metric="precomputed",
                verbose=9,
                n_iter=4000,
                init="random")
            self.embedding = self.seq_tsne.fit_transform(self.matrix)
            log.info("Memory usage after embedding fit %gMB" %
                     (util.memory_usage()))

        self.embedding = pd.DataFrame(
            self.embedding, columns=self.coord_dims, index=self.sequences[0])

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
            assert self.distances.dtype == "uint8"
            assert self.distances.ndim == 1
            self.distances.tofile(outf)
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
            kdims=self.coord_dims,
            vdims=["Sequence"],
            label="Sequences",
        ).opts(plot=dict(tools=[hover, "box_select"])).opts(plot=dict(
            width=800, height=800))
        #p.pprint()

        return p
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
