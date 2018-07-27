
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

    def __init__(self, input_file, enrichment_column="mean_ln_fold"):
        import holoviews as hv
        hv.extension('bokeh')
        self.enrichment_column = enrichment_column
        self.read_kmer_data(input_file)


        
        
    def read_config(self, config_file, distance_file=None):
        "Read configuration file for 'old style' input that computes the enrichments etc. by itself."
        import json
        import pandas as pd
        import os.path

        if not os.path.exists(distance_file):
            raise ValueError(
                "Can't find distance file {}! You might want to download http://www.cs.helsinki.fi/u/kpalin/all8mers_min_rev_complement.dists.gz".
                format(distance_file))

        self.config = json.load(open(config_file))
        log.info(str(self.config))

        self._binders = [x for x in self.config.keys() if x != "_config"]
        if len(self._binders) != 1:
            self.sole_binder = self._binders[0]
            
            log.warning(
                "Using only {} as binder. Can't currently handle more (given:  {})".
                format(self.sole_binder, str(self._binders)))

        else:
            self.sole_binder = self._binders[0]

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

    def read_kmer_data(self,input_file,binder="Default"):
        "Read a tab separated input file"
        import pandas as pd
        import logging as log
        from .tsne_mapper import TsneMapper

        self.data = dict()
        data = pd.read_table(input_file,sep="\t")
        kmer_name = data.columns[0]
        if not data[kmer_name].str.contains("^[ACGT]+$").all():
            idx = ~data[kmer_name].str.contains("^[ACGT]+$")
            false_kmer_eg = data.loc[idx,kmer_name].head()
            log.critical("First column ({}) contains non ACGT kmers (e.g. {})!".format(kmer_name," ".join(false_kmer_eg)))
        else:
            data = data.set_index(kmer_name)

        self.data[binder] = data
        self.data = pd.concat(self.data)
        self.selected_kmers = self.data   # Tautology
        
        #self.sequences = pd.DataFrame(self.selected_kmers.index.get_level_values(1).to_series())
        self.tsne_obj = dict()

        tsne_obj = TsneMapper()
        tsne_obj.sequences=pd.DataFrame({0:data.index.to_series()})
        log.info(tsne_obj.sequences)
        self.tsne_obj[binder] = tsne_obj

        log.info("Sequences: "+str(self.tsne_obj[binder].sequences.head()))
        self.N = len(self.tsne_obj[binder].sequences)
        if self.N<=1:
            log.critical("At most one sequence found! Will not work.")
            import sys
            sys.exit(1)
        if self.N<10:
            log.warning("Found very few sequences: "+str(self.tsne_obj[binder].sequences))

        self.sole_binder = binder
        self._binders = [ binder ]


        self.kmer_size = max(len(x) for x in data.index)


        self.local_maxima = dict()

        self.local_maxima[binder] = self.get_local_maxima(binder)


    def get_local_maxima(self, binder, n=10, max_local_dist=1.0):
        "Find the local (hudding space) maxima (enrichment) kmers"
        import logging as log

        import pandas as pd

        huddinge_mat = self.tsne_obj[binder].matrix

        skmers = self.selected_kmers.loc[binder,self.enrichment_column]
        print(skmers.head())
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


                rep_maxima.loc[[kmer]] = skmers.loc[c_reps].idxmax()

                #print(rep_maxima.loc[list(allocated_neighbours)])
                #log.info("Skipping non maximal locus {}".format(kmer))
        rep_maxima.index.name = "kmer"
        rep_maxima.name = "representative"
        assert (skmers.loc[rep_maxima.index].values <=
                skmers.loc[rep_maxima].values).all()

        log.info(str(skmers.loc[local_maxima]))
        return local_maxima, rep_maxima

    def circle_map_anchors(self, binder, anchors, seed=0):
        "Find radial (r=enrichment, theta = x) placing for the given anchors. Trying to match euclidean and hudding distance"
        import logging as log
        import scipy.spatial.distance as ssd
        import numpy as np
        import pandas as pd

        np.random.seed(seed)

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
            cartes = np.vstack([[0.0,r],cartes]) # Top anchor always faces north.
            D = ssd.pdist(cartes, "euclidean")
            return D

        def ssq_distance_diff(r_thetas):
            return ((circle_distances(r_thetas) - jittered_D)**2).sum()

        import scipy.optimize as so
        init_params = np.append(
            np.array([np.max(jittered_D) / 2]),
            np.random.uniform(0, 2 * np.pi, size=len(anchors)-1))

        bounds = [(0.1, None)] + [(0.0, 2.0 * np.pi)] * (len(anchors)-1)

        #opt = so.minimize(ssq_distance_diff,x0=init_params,bounds = bounds)
        opt = so.basinhopping(
            ssq_distance_diff,
            x0=init_params,
            minimizer_kwargs=dict(bounds=bounds))

        log.info(opt.message)
        r = opt.x[0]
        thetas = pd.Series(np.concatenate(([np.pi/2.0],opt.x[1:])), index=anchors)
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
        import logging as log

        from holoviews import streams

        if binder is None:
            binder = self.sole_binder
        else:
            assert binder in self._binders

        
        points, enrichment_r = self.plot_points(binder)





        # Declare points as source of selection stream
        selection = streams.Selection1D(source=points)

        # Write function that uses the selection indices to slice points and compute stats
        def selected_histogram(index):
            selected = points.iloc[index]
            if index:
                #label = str(selected.dframe().mean(axis=0))[:15]
                label = "Mean {} {}: {:.3g}".format(
                    binder, self.enrichment_column, selected.dframe()[self.enrichment_column].mean())
                #label = 'Mean x, y: %.3f, %.3f' % tuple(selected.array().mean(axis=0))
            else:

                selected = points
                label = 'No selection'
            from holoviews.operation import histogram

            h = histogram(
                selected, dimension=self.enrichment_column,
                dynamic=False).relabel(label)  #.opts(style=dict(color='red'))
            return h

        def selected_table(index):
            selected = points.iloc[index]
            if index:
                label = "Mean {} {}: {:.3g}".format(
                    binder, self.enrichment_column, selected.dframe()[self.enrichment_column].mean())
            else:

                selected = points
                label = 'No selection'

            t = selected.table()
            html = t.dframe().sort_values(
                self.enrichment_column, ascending=False).head(50).drop(
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
            d = selected.data.sort_values(self.enrichment_column, ascending=False).kmer
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
        r = points.hist(dimension=self.enrichment_column) * hv.Ellipse(
            0, 0, ellipse_diameter) << hv.DynamicMap(
                selected_histogram, streams=[selection])
        #r = r.relabel(binder)
        r= r+hv.DynamicMap(selected_table, streams=[selection]) + \
            hv.DynamicMap(selected_matrix, streams=[selection]) + \
            hv.DynamicMap(selected_heatmap, streams=[selection])

        return r.cols(2)

    def save_bokeh_points(self,output_file_name,binder=None):
        """Write the points plot to output_file_name.html. This plot
        can be included in other HTML pages."""
        import holoviews as hv
        from bokeh.io import save
        from bokeh.resources import CDN
        import logging as log

        from bokeh.embed import autoload_static

        if binder is None:
            binder = self.sole_binder

        renderer = hv.renderer('bokeh')

        points, enrichment_r = self.plot_points(binder, extra_tools=[])
        ellipse_diameter = enrichment_r.max() * 2

        points = points* hv.Ellipse(0, 0, ellipse_diameter) 
        #points = points.relabel(binder)
        # Convert to bokeh figure then save using bokeh
        plot = renderer.get_plot(points).state
        plot.legend.visible = False


        t="""<html><h1>Enrichment plot for {binder}</h1>
        The plot displays 10 locally maximal kmers and their neighbourhoods. 
        The locality is defined by huddinge distance and the maximality by {enrichment_column}.
        In the plot, distance from the origin is the {enrichment_column} value for each point. 
        The 10 local maxima are set to locations minimizing the sum of squared error 
        between their Euclidean distance on 2D plot and their Huddinge distances.
        The non-maximal points are allocated to the representative local maxima reached by up-hill climbing 
        from the point in question. The point color represents its distance to 
        its representative local maxima. There is random i.i.d jitter on the angle of the 
        non-maximal points. The circle is at the globally maximal enrichment<p>
        {tag}
        This file is generated by {pkg.key} {pkg.version} at {time}
        </html>"""
        
        from datetime import datetime
        import os.path
        import pkg_resources
        js_name   = "./{}.js".format(output_file_name)
        html_name = "./{}.html".format(output_file_name)
        js, tag = autoload_static(plot, CDN, os.path.basename(js_name))

        with open(js_name,"w") as js_f, open(html_name,"w") as html_f:
            js_f.write(js)

            pkg_resources.get_distribution("huddinge_tsne_browser")



            html_f.write(t.format(binder=binder,
                enrichment_column=self.enrichment_column,
                tag=tag,
                pkg = pkg_resources.get_distribution("huddinge-tsne-browser"),
                time = datetime.now().isoformat()
                ))

        print(tag)
        log.info("Wrote {} and {}.".format(js_name, html_name))
        #print("#"*80)



    def plot_points(self, binder, extra_tools =  ['box_select', 'lasso_select']):
        import pandas as pd
        import logging as log
        import numpy as np
        import holoviews as hv
        from bokeh.models import HoverTool


        loc_max, representatives = self.local_maxima[binder]
        loc_max = loc_max[:10]
        log.info("local maxima: {}".format(str(loc_max)))

        _, theta_angle, _, _, _ = self.circle_map_anchors(binder, loc_max)


        huddinge_mat = self.tsne_obj[binder].matrix

        # Position anchors
        enrichment_r = self.data.loc[binder].loc[theta_angle.index, self.enrichment_column]
        anchors_x = polar2cartesian(enrichment_r, theta_angle)
        anchors_x["enrichment"] = enrichment_r

        # Select the closest (in hudding distance) local maxima as representative. 
        # Break ties according to enrichment

        single_rep = pd.DataFrame(representatives)
        single_rep["distance"] = [
            huddinge_mat.at[x, y]
            for x, y in single_rep.representative.iteritems()
        ]


        single_rep = single_rep.join(self.data.loc[binder])
        single_rep.columns = map(str, single_rep.columns)
        single_rep["theta"] = theta_angle[single_rep.representative].values

        # Add angular jitter
        jitter_span = 2 * np.pi / 200.0 * single_rep.distance

        # TODO: instead of jitter, bias the kmers to "closest" maxima direction
        np.random.seed(0)
        single_rep[["x", "y"]] = polar2cartesian(
            single_rep[self.enrichment_column], single_rep.theta + np.random.uniform(
                low=-jitter_span, high=jitter_span, size=len(jitter_span)))


        p_data = single_rep.reset_index()


        from .util import isidentifier
        col_renamer = {x:"A"+x.translate(str.maketrans(",() '`'","_______")).strip("_")+"B" for x in p_data.columns if not isidentifier(x)}


        col_renamer_r = {y:x for x,y in col_renamer.items()}

        p_data.columns = [col_renamer.get(x,x) for x in p_data.columns]

        log.info(p_data.columns)
        hover = HoverTool(tooltips=[(col_renamer_r.get(c,c), "@{"+c+"}") for c in p_data.columns \
                                         if c not in ("x","y","distance","theta")])

        points = hv.Points(
            p_data,
            kdims=["x", "y"],
            vdims=[_x for _x in p_data.columns if _x not in ["x", "y"]],
            extents=(-1, -1, 1, 1)).opts(
                plot=dict(
                    tools=[hover ] + extra_tools,
                    width=600,
                    height=500,
                    scaling_factor=13,
                    size_index=self.enrichment_column,
                    bgcolor="lightgray",
                    show_grid=True,
                    color_index="distance",
                    colorbar=True,
                    colorbar_opts=dict(title="Distance to rep"),
                    #  cmap='PiYG', color_levels=int(max_distance+1.1)
                ),
                style=dict(cmap="inferno"),
                norm=dict(axiswise=True, framewise=False))
        return points.relabel(binder), enrichment_r

