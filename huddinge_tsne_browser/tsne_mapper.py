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

    def __init__(self, input_file=None, force_distances=False):
        """
        
        Arguments:
        - `input_file`:
        """

        self._input_file = input_file
        self.data_dims = []
        self.coord_dims = ["tsne0", "tsne1"]

        if input_file is not None:
            self.read_data(force_distances)

    def __len__(
            self, ):
        """Return number of kmers
        """
        return len(self.sequences)

    def read_data(self, force_distances=False):
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
            self.fit_measure_ = float(p[1])

        l = [fin.readline() for _ in range(self.N)]
        bytes_read += sum(len(x) for x in l)

        self.sequences = pd.read_table(
            StringIO(u"".join(l)),
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

        if (not hasattr(self, "embedding")) or force_distances:
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

    def focus_sequences(self, prop_signal=0.1, prop_background=0.01):
        "Return list of sequences including top prop_signal of highest count sequences for each data value and  prop_background of random sequences"
        seqs = set(
            self.sequences.sample(int(prop_background * len(self.sequences)))[
                0])
        lims = self.embedding[self.data_dims].quantile(1 - prop_signal)
        for d, l in lims.iteritems():
            s = self.embedding.loc[self.embedding[d] >= l].index
            seqs.update(s)

        return list(sorted(seqs))

    def subset_sequences(self, seqs):
        "Drop all information about other sequences but seqs"
        keepers = self.sequences[0].isin(seqs).nonzero()[0]

        log.info("Keeping a subset of {} sequences.".format(len(keepers)))
 
        if hasattr(self, "_matrix"):
            del (self._matrix)

        self.sequences = self.sequences.iloc[keepers]
        self.N = len(self.sequences)

        if hasattr(self, "embedding"):
            self.embedding = self.embedding.reindex(index=seqs)

    def add_kmercounts(self, name, filename):
        """
        
        Arguments:
        - `name`:
        - `filename`:
        """
        _, counts = read_jf(filename)

        try:
            self.embedding[name] = counts.reindex(
                index=self.embedding.index, fill_value=0)
            if name not in self.data_dims:
                self.data_dims.append(name)
        except KeyError as e:
            kmer_lens = sorted(set(len(x) for x in counts.index))
            log.warning(
                "Coudn't add %s kmers. They are of length %s while the embedding is for kmers of length %d"
                % (name, ",".join(str(x) for x in kmer_lens), self.kmer_size))

    def read_distances(self, fin):
        """Read distance data from open file fin

        Doesn't do anything since the distances are computed online.

        Arguments:
        - `fin`:
        """
        pass

    def _get_kmer_size(self):
        return len(self.sequences[0])

    kmer_size = property(_get_kmer_size)

    def clear_layout(self):
        "Remove layout information"
        self.sequences = self.sequences[[0]]
        self.embedding = self.embedding.drop(self.coord_dims, axis=1)
        del (self.fit_measure_)

    def laidout(
            self, ):
        """Are the sequences laid out properly
        """
        return hasattr(self, "embedding") and set(self.coord_dims).issubset(
            self.embedding.columns)

    def _get_matrix(self):
        if not hasattr(self, "_matrix"):
            import pandas as pd
            import numpy as np
            import pyhuddinge as ph

            log.info("Memory usage before matrix formatting %gMB" %
                     (util.memory_usage()))

            distances = ph.all_pairs_huddinge_distance(self.sequences[0],reverse_complements=True)

            M = ph.squareform(distances,self.sequences[0])
            self._matrix = M

            log.info("Memory usage after matrix formatting %gMB" %(util.memory_usage()))

        return self._matrix

    matrix = property(_get_matrix)

    def compute_tsne(self, perplexity=30.0, fake=False):
        """
        """
        if self.laidout():
            self.clear_layout()
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

        self.coord_dims = ["tsne0", "tsne1"]

        self.embedding = pd.DataFrame(
            self.embedding, columns=self.coord_dims, index=self.sequences[0])

        self.embedding.index.name = "Sequence"
        self.fit_measure_ = self.seq_tsne.kl_divergence_

    def compute_mds(self):
        """
        """
        if self.laidout():
            self.clear_layout()

        import pandas as pd
        from sklearn.manifold import MDS
        self.seq_mds = MDS(dissimilarity="precomputed", verbose=1, n_jobs=-2)
        self.embedding = self.seq_mds.fit_transform(self.matrix)
        log.info("Memory usage after embedding fit %gMB" %
                 (util.memory_usage()))

        self.coord_dims = ["mds0", "mds1"]

        self.embedding = pd.DataFrame(
            self.embedding, columns=self.coord_dims, index=self.sequences[0])

        self.embedding.index.name = "Sequence"
        self.fit_measure_ = self.seq_mds.stress_

    def _get_adjacency_matrix(self, random_reconnect=True):
        """Return adjacency matrix for the kmers. Possibly connect singleton kmers to random closeby kmer 
        
        Arguments:
        - `random_reconnect`:
        """
        # Compute adjacency matrix
        import logging
        import numpy as np
        import pandas as pd

        # Kmers on huddinge distance 1 are adjacent and distance two are half way.
        adjacency = np.zeros(self.matrix.shape)
        is_smallish = self.matrix.values < 2.5
        adjacency[is_smallish] = 1.0 / self.matrix.values[is_smallish]
        np.fill_diagonal(adjacency, 0.0)

        disconnected = np.nonzero(~adjacency.any(axis=0))[0]

        if len(disconnected) > 0:
            if random_reconnect:
                # Kmers further away are adjacent to a random kmer at huddinge distance 2
                dis_idx, dis_adj = np.nonzero(
                    (self.matrix.iloc[disconnected] == 2).values)

                import numpy as np
                rand_neighbor = pd.Series(dis_adj).groupby(dis_idx).apply(
                    lambda x: np.random.choice(x, 1)[0])

                x, y = disconnected[rand_neighbor.index].astype(
                    int), rand_neighbor.values.astype(int)
                assert (self.matrix.iloc[x, y] == 2).all().all()

                adjacency[x, y] = True
                adjacency[y, x] = True

            else:
                logging.info("Dropping %d disconnected kmers" %
                             (len(disconnected)))
                x = adjacency.any(axis=0)
                adjacency = adjacency[x][:, x]

        assert (adjacency == adjacency.T).all()

        disconnected = np.nonzero(~adjacency.any(axis=0))[0]

        if len(disconnected) > 0:
            import warnings
            warnings.warn("Disconnected kmers in adjacency graph")
        return adjacency

    def compute_spectral(self):
        """
        """
        if self.laidout():
            self.clear_layout()

        import pandas as pd
        from sklearn.manifold import SpectralEmbedding
        self.seq_spectral = SpectralEmbedding(affinity="precomputed")

        embedding = self.seq_spectral.fit_transform(
            self._get_adjacency_matrix().astype(float))

        log.info("Memory usage after embedding fit %gMB" %
                 (util.memory_usage()))

        self.coord_dims = ["spectral0", "spectral1"]
        embedding = pd.DataFrame(
            embedding, columns=self.coord_dims, index=self.sequences[0])

        embedding.index.name = "Sequence"

        if len(self.data_dims) > 0:
            self.embedding = embedding.join(self.embedding[self.data_dims])
        else:
            self.embedding = embedding

        self.fit_measure_ = -1.0
        self.data_dims = []

    def write_data(self, outfile):
        """
        
        Arguments:
        - `outfile`:
        """
        raise NotImplementedError("Don't try to do this.")

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
