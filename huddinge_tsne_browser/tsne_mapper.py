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


def polar2cartesian(r, thetas):
    import numpy as np
    import pandas as pd
    import scipy.spatial.distance as ssd

    cartes = np.array([r * np.cos(thetas), r * np.sin(thetas)]).T
    if hasattr(thetas, "index"):
        cartes = pd.DataFrame(cartes, columns=["x", "y"], index=thetas.index)
    return cartes


def normalize_selex(kmers, pseudocount=1.0):
    """Normalize the kmer counts and return

    1. Normalized counts (count of a kmer divided by median count on that cycle)
    2. Natural logarithm of fold change of normalized kmer count between each cycle
    3. Mean of the log fold changes
    4. z-score of the log fold change (mean divided by standard deviation)


    kmer.shape == (kmers,cycles).

    Return value is a tuple of
    (median normalized counts,
       log fold change for each cycle,
       mean fold change per cycle,
       fold change z-score) """

    import pandas as pd
    import numpy as np

    assert kmers.shape[0] > kmers.shape[1]
    kmers = kmers.sort_index(axis=1) + pseudocount

    # Median normalisation
    norm_cnt = kmers / kmers.median(axis=0)

    # log fold change per cycle
    ln_fold_change = np.log(norm_cnt).diff(axis=1)
    ln_fold_change = ln_fold_change.drop(ln_fold_change.columns[0], axis=1)

    # Mean fold change per cycle
    mean_fold = ln_fold_change.mean(axis=1)

    fold_z = mean_fold / ln_fold_change.std(axis=1)

    d = {}
    for i, c in enumerate(norm_cnt.columns):
        d[(c, "normalized")] = norm_cnt[c]
        if i > 0:
            d[(c, "ln_fold")] = ln_fold_change[c]

    d["mean_ln_fold"] = mean_fold
    d["fold_z"] = fold_z
    return pd.DataFrame(d)


class PolarMapper(object):
    "Lay out kmers according to polar coordinate setup"

    def __init__(self, distance_file, config_file):
        import json
        import pandas as pd
        import os.path

        if not os.path.exists(distance_file):
            raise ValueError(
                "Can't find distance file {}! You might want to download http://www.cs.helsinki.fi/u/kpalin/all8mers_min_rev_complement.dists.gz".
                format(distance_file))

        self.config = json.load(open(config_file))
        log.info(str(self.config))

        self.sole_binder = [x for x in self.config.keys() if x != "_config"]
        if len(self.sole_binder) != 1:
            import sys
            log.critical(
                "Currently config file can only take one binder at a time. Found {}".
                format(str(self.sole_binder)))
            sys.exit(1)
        else:
            self.sole_binder = self.sole_binder[0]

        self.kmer_size = self.config["_config"].setdefault("kmer_size", 8)

        # Load
        kmers = dict()
        for binder, files in self.config.items():
            if binder == "_config": continue
            for f in files:
                log.info("Loading cycle {} of {} from {}".format(f[
                    "cycle"], binder, f["filename"]))
                _, kmers[(binder, f["cycle"])] = read_jf(f["filename"])

        self.kmers = pd.DataFrame(kmers)
        del (kmers)

        # Normalise
        pseudocount = self.config["_config"].setdefault("pseudocount", 1.0)

        self.data = dict()
        for binder, files in self.config.items():
            if binder == "_config": continue
            self.data[binder] = normalize_selex(self.kmers[binder],
                                                pseudocount)
        self.data = pd.concat(self.data)

        # Select enriched kmers
        zscore_limit = self.config["_config"].setdefault("zscore", 2.58)
        self.selected_kmers = self.data.loc[self.data["fold_z"] > zscore_limit,
                                            "mean_ln_fold"]

        # Read distances etc.
        self.local_maxima = dict()
        self.tsne_obj = dict()

        for binder, files in self.config.items():
            if binder == "_config": continue
            tsne_obj = TsneMapper(distance_file, force_distances=True)
            tsne_obj.subset_sequences(
                list(self.selected_kmers.loc[binder].index))
            self.tsne_obj[binder] = tsne_obj

            self.local_maxima[binder] = self.get_local_maxima(binder)

    def get_local_maxima(self, binder, n=10, max_local_dist=1.0):
        "Find the local (hudding space) maxima (enrichment) kmers"

        import pandas as pd

        huddinge_mat = self.tsne_obj[binder].matrix

        skmers = self.selected_kmers.loc[binder]

        candidates = skmers.sort_values(ascending=False).iteritems()
        #Find local maxima
        local_maxima = []
        rep_maxima = pd.Series(None, index=skmers.index)
        for kmer, v in candidates:
            if rep_maxima.loc[[kmer]].notnull()[0]: continue
            neighbours = skmers[huddinge_mat.loc[kmer] <= max_local_dist]
            neighbours = neighbours.drop(kmer)
            if len(neighbours) == 0:
                log.info("No {:.0f} hudding distance neighbours for kmer {}".
                         format(max_local_dist, kmer))
                rep_maxima.loc[kmer] = kmer
                continue
            if neighbours.max() < v:
                # All neighbours are lower, hence we have local maxima
                local_maxima.append(kmer)
                new_neighbours = set(neighbours.index) - set(
                    rep_maxima.loc[rep_maxima.notnull()].index)
                rep_maxima.loc[list(new_neighbours) + [kmer]] = kmer

            else:
                # Not local maxima. Allocate self to highest rep neighbour.
                allocated_neighbours = set(
                    neighbours[neighbours >= v].index) & set(
                        rep_maxima.loc[rep_maxima.notnull()].index)
                c_reps = rep_maxima.loc[list(allocated_neighbours)]

                if c_reps.nunique() != 1:
                    log.info(str(skmers.loc[c_reps.drop_duplicates()]))

                #new_neighbours = set(neighbours.index) - set(
                #    rep_maxima.loc[rep_maxima.notnull()].index)
                rep_maxima.loc[  #list(new_neighbours) +
                    [kmer]] = skmers.loc[c_reps].argmax()

                #print(rep_maxima.loc[list(allocated_neighbours)])
                #log.info("Skipping non maximal locus {}".format(kmer))
        rep_maxima.index.name = "kmer"
        rep_maxima.name = "representative"
        assert (skmers.loc[rep_maxima.index].values <=
                skmers.loc[rep_maxima].values).all()

        log.info(str(skmers.loc[local_maxima]))
        return local_maxima, rep_maxima

    def circle_map_anchors(self, binder, anchors):
        "Find radial (r=enrichment, theta = x) placing for the given anchors. Trying to match euclidean and hudding distance"
        import logging as log
        import scipy.spatial.distance as ssd
        import numpy as np
        import pandas as pd

        jitter = np.random.uniform(
            -0.1, 0.1, size=(len(anchors), len(anchors)))
        np.fill_diagonal(jitter, 0)

        jittered_dist = self.tsne_obj[binder].matrix.loc[anchors, anchors] + (
            jitter + jitter.T) / 2.0
        jittered_D = ssd.squareform(jittered_dist)

        def circle_distances(r_thetas):

            import numpy as np
            import scipy.spatial.distance as ssd

            r, thetas = r_thetas[0], r_thetas[1:]
            cartes = np.array([r * np.cos(thetas), r * np.sin(thetas)]).T

            D = ssd.pdist(cartes, "euclidean")
            return D

        def ssq_distance_diff(r_thetas):
            return ((circle_distances(r_thetas) - jittered_D)**2).sum()

        import scipy.optimize as so
        init_params = np.append(
            np.array([np.max(jittered_D) / 2]),
            np.random.uniform(0, 2 * np.pi, size=len(anchors)))

        bounds = [(0.1, None)] + [(0.0, 2.0 * np.pi)] * len(anchors)

        #opt = so.minimize(ssq_distance_diff,x0=init_params,bounds = bounds)
        opt = so.basinhopping(
            ssq_distance_diff,
            x0=init_params,
            minimizer_kwargs=dict(bounds=bounds))

        log.info(opt.message)
        r = opt.x[0]
        thetas = pd.Series(opt.x[1:], index=anchors)
        D = pd.DataFrame(
            ssd.squareform(circle_distances(opt.x)),
            index=anchors,
            columns=anchors)
        _n = len(anchors) * (len(anchors) - 1) / 2
        return r, thetas, D, opt.fun / _n, ((D - jittered_dist)
                                            **2).sum().sum() / 2

    def plot_polar(self, binder=None, theta_angle=None):
        "Arguments: Name of the binder tf and list of representative kmers (with angle positions)"
        import numpy as np
        import holoviews as hv

        if binder is None:
            binder = self.sole_binder
        else:
            assert binder == self.sole_binder

        tf = binder

        loc_max, representatives = self.local_maxima[binder]
        loc_max = loc_max[:10]
        log.info("local maxima: {}".format(str(loc_max)))

        _, theta_angle, _, _, _ = self.circle_map_anchors(binder, loc_max)

        from holoviews import streams

        huddinge_mat = self.tsne_obj[binder].matrix

        # Position anchors
        enrichment_r = self.data.loc[tf].loc[theta_angle.index, "mean_ln_fold"]
        anchors_x = polar2cartesian(enrichment_r, theta_angle)
        anchors_x["enrichment"] = enrichment_r

        # Select the closest (in hudding distance) local maxima as representative. 
        # Break ties according to enrichment
        import pandas as pd
        single_rep = pd.DataFrame(representatives)
        single_rep["distance"] = [
            huddinge_mat.at[x, y]
            for x, y in single_rep.representative.iteritems()
        ]

        #single_rep["enrichment"] = self.data.loc[tf].loc[single_rep.index,
        #                                                 "mean_ln_fold"].values
        ENRICHMENT = "mean_ln_fold"
        single_rep = single_rep.join(self.data.loc[tf])
        single_rep.columns = map(str, single_rep.columns)
        single_rep["theta"] = theta_angle[single_rep.representative].values

        # Add angular jitter
        jitter_span = 2 * np.pi / 200.0 * single_rep.distance

        single_rep[["x", "y"]] = polar2cartesian(
            single_rep[ENRICHMENT], single_rep.theta + np.random.uniform(
                low=-jitter_span, high=jitter_span, size=len(jitter_span)))

        # Declare some points
        x = single_rep.reset_index()
        max_distance = x.distance.max()
        print(max_distance)
        #x["distance"] = x.distance.apply("{:.0f}".format)
        points = hv.Points(
            x,
            kdims=["x", "y"],
            vdims=[_x for _x in x.columns if _x not in ["x", "y"]],
            extents=(-1, -1, 1, 1)).opts(
                plot=dict(
                    tools=["hover", 'box_select', 'lasso_select'],
                    width=600,
                    height=500,
                    scaling_factor=13,
                    size_index=ENRICHMENT,
                    bgcolor="lightgray",
                    show_grid=True,
                    color_index="distance",
                    colorbar=True,
                    colorbar_opts=dict(title="Distance to rep"),
                    #  cmap='PiYG', color_levels=int(max_distance+1.1)
                ),
                style=dict(cmap="inferno_r"),
                norm=dict(axiswise=True, framewise=False))
        #points = points.options(background_fill_color="ligthgray")
        #points.opts(plot=dict(aspect_weight=1.0,weight=1.0))

        #points = hv.Points(tsneD[tf].embedding.join(data.loc[tf]).reset_index(),
        #                   kdims=kdims,
        #                   vdims=['mean_ln_fold',"Sequence"]).relabel(tf)

        # Declare points as source of selection stream
        selection = streams.Selection1D(source=points)

        # Write function that uses the selection indices to slice points and compute stats
        def selected_histogram(index):
            selected = points.iloc[index]
            if index:
                #label = str(selected.dframe().mean(axis=0))[:15]
                label = "Mean {} {}: {:.3g}".format(
                    tf, ENRICHMENT, selected.dframe()[ENRICHMENT].mean())
                #label = 'Mean x, y: %.3f, %.3f' % tuple(selected.array().mean(axis=0))
            else:

                selected = points
                label = 'No selection'
            from holoviews.operation import histogram

            h = histogram(
                selected, dimension=ENRICHMENT,
                dynamic=False).relabel(label)  #.opts(style=dict(color='red'))
            return h

        def selected_table(index):
            selected = points.iloc[index]
            if index:
                label = "Mean {} {}: {:.3g}".format(
                    tf, ENRICHMENT, selected.dframe()[ENRICHMENT].mean())
            else:

                selected = points
                label = 'No selection'

            t = selected.table()
            html = t.dframe().sort_values(
                ENRICHMENT, ascending=False).head(50).drop(
                    ["x", "y"], axis=1).set_index("kmer").to_html()
            return hv.Div("<div>{}</div>".format(html)).opts(plot=dict(
                width=200))

        def selected_matrix(index):
            import pandas as pd
            import holoviews as hv
            selected = points.iloc[index]
            if not index:
                selected = points
            d = selected.data.kmer

            if len(d) < 50:
                counts = align_logo(list(d)).T
            else:
                counts = pd.DataFrame(tuple(x) for x in d).apply(
                    lambda x: x.value_counts(), axis=0)
                #counts = counts.fillna(0)        

            counts = counts.reindex(
                columns=pd.RangeIndex(-counts.shape[1] / 2,
                                      int(counts.shape[1] * 1.5)),
                fill_value=0)
            counts.columns = map(str, counts.columns)

            ##return hv.Div(counts.astype(int).to_html(notebook=True)).opts(plot=dict(width=400,height=700))
            return hv.Table(counts)

        def align_logo(seqs):
            "Align fixed length kmers to the first one and generate a count matrix"
            import pandas as pd
            seqs = [
                pd.get_dummies(
                    pd.Categorical(list(x), categories=["A", "C", "G", "T"]))
                for x in seqs
            ]

            m_len = len(seqs[0])
            motif = seqs[0].reindex(
                index=pd.RangeIndex(-m_len, m_len * 2 - 1), fill_value=0)
            for s in seqs[1:]:
                score, shift = max(((motif.shift(shift) * s).sum().sum(),
                                    shift)
                                   for shift in range(-m_len + 1, m_len - 1))

                rc_s = s.copy()
                rc_s.index = list(s.index)[::-1]
                rc_s.columns = list(s.columns)[::-1]

                rc_score, rc_shift = max(
                    ((motif.shift(shift) * rc_s).sum().sum(), shift)
                    for shift in range(-m_len + 1, m_len - 1))
                if rc_score > score:
                    print("rc better")
                    shift, score, s = rc_shift, rc_score, rc_s.sort_index()
                motif.loc[-shift:-shift + m_len - 1] += s.set_index(
                    pd.RangeIndex(-shift, -shift + m_len))

            return motif.loc[motif.sum(axis=1) > 0].copy()

        def show_logo(mat, height=100, width=400):
            import svgwrite as sw
            #assert len(mat)==4
            mat = (mat + 1.0 / (mat.shape[0] * mat.shape[1]))
            mat = mat / mat.sum().max()
            w = mat.shape[1] * 15
            #        dr=sw.Drawing(size=("%dpt"%(w),"30pt"))
            dr = sw.Drawing(size=("{:d}pt".format(width),
                                  "{:d}pt".format(height)))
            #dr.add(dr.rect((0,0),(4100,1130),fill="pink"))
            g = sw.container.Group()
            cmap = dict(zip("ACGT", ["green", "blue", "orange", "red"]))

            for c, clab in enumerate(mat.columns):

                xpos = 10 * int(c)
                ypos = 10
                #g.add(dr.text(str(c),fill="black" )).translate(xpos,ypos)

                for i, base in enumerate(mat.index):
                    yscale = mat.loc[base, clab]

                    t = g.add(dr.text(base, fill=cmap.get(base)))
                    t.translate(xpos, ypos)
                    t.scale(sx=1, sy=yscale)
                    ypos -= (10 * yscale)
            dr.add(g)
            g.scale(sx=width / (len(mat.columns) * 10.0), sy=height / (10.0))
            return hv.Div(dr.tostring()).opts(plot=dict(
                width=width * 2, height=height))

        def selected_heatmap(index):
            import pandas as pd
            import holoviews as hv
            selected = points.iloc[index]
            if not index:
                selected = points
            d = selected.data.sort_values(ENRICHMENT, ascending=False).kmer
            if len(d) < 50:
                counts = align_logo(list(d)).T
            else:
                counts = pd.DataFrame(tuple(x) for x in d).apply(
                    lambda x: x.value_counts(), axis=0)
                #counts = counts.fillna(0)        
            #        counts = pd.DataFrame(tuple(x) for x in d).apply(lambda x:x.value_counts(),axis=0)
            counts = counts.fillna(0)
            try:
                return show_logo(counts, height=100)
            except (ImportError, AttributeError):
                counts = counts.stack().reset_index()
                counts.columns = map(str, counts.columns)

                t = hv.Table(counts, kdims=["level_1", "level_0"])
                t = t.redim(level_1="Position", level_0="Base")

                return t.to.heatmap().opts(plot={
                    "colorbar": True,
                    "tools": ["hover"],
                    "invert_yaxis": True,
                    "sizing_mode": "fixed",
                    "width": 400,
                    "height": 200
                })

        # Combine points and DynamicMap

        ellipse_diameter = enrichment_r.max() * 2
        r = points.hist(dimension=ENRICHMENT) * hv.Ellipse(
            0, 0, ellipse_diameter) << hv.DynamicMap(
                selected_histogram, streams=[selection])
        r = r.relabel(tf)
        r= r+hv.DynamicMap(selected_table, streams=[selection]) + \
            hv.DynamicMap(selected_matrix, streams=[selection]) + \
            hv.DynamicMap(selected_heatmap, streams=[selection])

        return r.cols(2)


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
        if hasattr(self, "distances"):
            import itertools as it
            import numpy as np
            ij = np.fromiter(
                it.combinations(keepers, 2), dtype=[('j', int), ('i', int)])
            ij.sort(order=["i", "j"])
            idx = ij["i"] * ((ij["i"] - 1.0) / 2.0) + ij["j"]

            idx = idx.astype(int)
            self.distances = self.distances[idx]

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
        
        Arguments:
        - `fin`:
        """
        import pandas as pd
        import numpy as np

        d_pos = fin.tell()
        d = np.fromfile(fin, dtype="uint8", count=-1)
        assert len(d) == self.N * (self.N - 1) / 2

        self.distances = d

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

            self._matrix = pd.DataFrame(
                self._matrix,
                index=self.sequences[0].squeeze(),
                columns=self.sequences[0].squeeze())

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
        is_smallish = self.matrix < 2.5
        adjacency[is_smallish] = 1.0 / self.matrix[is_smallish]
        np.fill_diagonal(adjacency, 0.0)

        disconnected = np.nonzero(~adjacency.any(axis=0))[0]

        if len(disconnected) > 0:
            if random_reconnect:
                # Kmers further away are adjacent to a random kmer at huddinge distance 2
                dis_idx, dis_adj = np.nonzero(self.matrix[disconnected] == 2)

                import numpy as np
                rand_neighbor = pd.Series(dis_adj).groupby(dis_idx).apply(
                    lambda x: np.random.choice(x, 1)[0])

                x, y = disconnected[rand_neighbor.index].astype(
                    int), rand_neighbor.values.astype(int)
                assert (self.matrix[x, y] == 2).all()

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

        log.info("Writing %s" % (outfile))
        with open(outfile, "w") as outf:
            outf.write("%d\t%g\n" % (len(self.sequences), self.fit_measure_))
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
