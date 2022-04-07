# -*- coding: utf-8 -*-
import numpy as np
import holoviews as hv

from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application
from bokeh.io import curdoc
from bokeh.layouts import layout
from bokeh.models import Slider, Button

from .datashaderselect import DataShaderSelect


def eq_hist_percentiles_symmetric(data, mask=None, nbins=256):
    """Return a numpy array after histogram equalization.
    For use in `shade`.
    Parameters
    ----------
    data : ndarray
    mask : ndarray, optional
       Boolean array of missing points. Where True, the output will be `NaN`.
    nbins : int, optional
        Number of bins to use. Note that this argument is ignored for integer
        arrays, which bin by the integer values directly.
    Notes
    -----
    This function is adapted from the implementation in scikit-image [1]_.
    References
    ----------
    .. [1] http://scikit-image.org/docs/stable/api/skimage.exposure.html#equalize-hist
    """
    print(np.histogram(data[~mask]))
    if not isinstance(data, np.ndarray):

        #raise TypeError("data must be np.ndarray")
        data2 = np.array(data)
        data = data2
    else:
        data2 = data if mask is None else data[~mask]
    if np.issubdtype(data2.dtype, np.integer):
        assert ValueError("Something funny with input")
    else:
        # # Split by percentile
        # n = 100
        # xpercent = 1 - np.logspace(0, -2, num=n)

        # bin_edges = np.unique(np.percentile(data2, np.linspace(0, 100)))
        # hist, _ = np.histogram(data2, bins=bin_edges)
        # cdf = hist.cumsum()
        # cdf = cdf / float(cdf[-1] + 1.0)
        # bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        hist, bin_edges = np.histogram(data2[data2 <= 0], bins=256)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        cdf = hist.cumsum()
        cdf = cdf / (2.0 * float(cdf[-1]))

        hist2, bin_edges2 = np.histogram(data2[data2 > 0], bins=256)
        bin_centers2 = (bin_edges2[:-1] + bin_edges2[1:]) / 2

        cdf2 = hist2.cumsum()
        cdf2 = cdf2 / (2.0 * float(cdf2[-1])) + 0.5

        bin_centers = np.concatenate([bin_centers, bin_centers2])
        cdf = np.concatenate([cdf, cdf2])
        #print(hist.sum(), hist2.sum(), data2.shape, cdf)
        #print(np.histogram(data2))

    #    print("Range:", data2.min(), data2.max())
    #    print(
    #        "bin_centers,value:",  #map("{0[0]:.3g}:{0[1]:.2g}".format,
    #        zip(bin_centers, cdf))
    out = np.interp(data.flat, bin_centers, cdf).reshape(data.shape)
    return out if mask is None else np.where(mask, np.nan, out)


class HuddingBrowser(object):
    """
    """

    def __init__(self, data):
        """
        
        Arguments:
        - `data`: an object with attribute "embedding" e.g. TsneMapper()
        """
        self._data = data

        if len(self._data.data_dims) == 0:
            self._points = hv.Points(
                self._data.embedding.reset_index(),
                #tsne.embedding.reset_index(),
                kdims=self._data.coord_dims,
                vdims=["Sequence"] + self._data.data_dims,
                label="Sequences", )

        else:
            overlay_dict = {}
            for d in self._data.data_dims:
                p = hv.Points(
                    self._data.embedding[self._data.coord_dims + [d]].rename(
                        columns={d: "Value"}).reset_index(),
                    kdims=self._data.coord_dims,
                    vdims=["Sequence", "Value"],
                    label="Sequences")
                overlay_dict[d] = p  #.opts(plot={"color_index"=d})
            #self._points = hv.NdOverlay(overlay_dict, kdims=["Counts"])
            self._points = hv.HoloMap(overlay_dict, kdims=["Counts"])

    def datashade(self, column=None, reduction=None):
        """Return datashader shaded with reduction of given column value
        
        Arguments:
        - `column`:
        - `reduction`:
        """
        import numpy as np
        import datashader as ds
        from bokeh import palettes
        from holoviews.operation.datashader import aggregate, shade, datashade, dynspread
        if reduction is None:
            reduction = ds.mean

        #print "datashade", column, reduction
        shade_opts = {"cmap": palettes.inferno(64)}
        shade_opts = {"cmap": palettes.Viridis256[::-1]}

        if column is not None:
            d = self._points.dframe()[column]
            d_min, d_max = np.percentile(d, [1.0, 99.0])  #.min(), d.max()
            #print "Value Range:", d_min, d_max
            #shade_opts["clims"] = (d_min, d_max)

            if d_min * d_max < 0.0:
                d_extreme = max(-d_min, d_max)
                print("Diverging palette", (-d_extreme, d_extreme),
                      (d_min, d_max))
                shade_opts = {
                    "cmap": palettes.RdYlBu11,
                    #"clims": (-d_extreme, d_extreme),
                    "normalization": "eq_hist"
                    #"normalization": eq_hist_percentiles_symmetric
                }

        def _linear_norm_debug(v, masked):
            #import pandas as pd
            import numpy as np
            #return np.log1p(np.where(masked, np.nan, v))
            #            print args
            #            print kwargs.keys()
            #v1d = v[~masked]
            try:
                min_v, max_v = np.nanmin(v.flat), np.nanmax(v.flat)
            except Exception:
                print("Poikkeus:", v)
                v = np.array(v)
                min_v, max_v = np.nanmin(v.flat), np.nanmax(v.flat)

            o = (v - (min_v - 1.0))
            print("o range after sub:", np.nanmin(o.flat), np.nanmax(o.flat))
            o /= (max_v - min_v)
            print(min_v, max_v, np.nanmin(o.flat), np.nanmax(o.flat))
            r = o
            #r = 50.0**(o - 1.0)
            print("Out", np.nanmin(r.flat), np.nanmax(r.flat))
            return np.where(masked, np.nan, r)

        #shade_opts["normalization"] = "linear"  #_percentiles
        #shade_opts["span"] = None
        #del (shade_opts["clims"])
        plot = dynspread(
            datashade(
                self._points,
                aggregator=ds.count() if column is None else reduction(column),
                **shade_opts))

        return plot

    def holoview_plot(
            self, ):
        """
        """
        import datashader as ds
        from holoviews.operation.datashader import aggregate, shade, datashade, dynspread
        from holoviews.streams import RangeXY

        self.ds_points = self.datashade("Value" if len(self._data.data_dims) >
                                        0 else None)
        self.ds_points = self.ds_points.opts(plot=dict(width=600, height=600))

        # Hover and zoom grid tool.
        self._hover_grid = hv.util.Dynamic(
            aggregate(
                self._points,
                aggregator=ds.mean("Value"),
                width=15,
                height=15,
                streams=[RangeXY(source=self.ds_points)]),
            operation=hv.QuadMesh).opts(
                plot=dict(tools=["hover"]),
                style=dict(alpha=0, hover_alpha=0.2))

        # Get the points in tapped rectangle
        self._posxy = DataShaderSelect(
            source=self._hover_grid,
            dataset=self._data.embedding,
            x=self._data.coord_dims[0],
            y=self._data.coord_dims[1])

        #self._posxy = hv.streams.Tap(source=self._hover_grid)

        def _dss_logger(**kwargs):
            import logging as log
            log.info("Handling event from datashader select: %s", str(kwargs))

        self._posxy.add_subscriber(_dss_logger)

        # Make layout
        self.tap_indicators = hv.DynamicMap(
            self.tap_points, kdims=[], streams=[self._posxy])
        self.selected_table = hv.DynamicMap(
            self.tap_table, streams=[self._posxy])

        self.tap_zoom = hv.DynamicMap(
            self.focus_plot, streams=[self._posxy],
            kdims=["Counts"]).opts(norm=dict(framewise=True)).redim.values(
                Counts=self._data.data_dims)

        def _h(Counts, index, fine_index, **kwargs):
            #print index, kwargs
            from holoviews.operation import histogram

            m = {
                Counts: "Value",
                "{}_frequency": "frequency",
            }
            if len(index) > 0:
                d = self._data.embedding.iloc[index]
                if len(fine_index) > 0:
                    try:
                        d = d.iloc[fine_index]
                    except IndexError:
                        pass

                label = "{} {} points".format(Counts, len(d))
                r = histogram(
                    hv.Points(d),
                    #self.selected_table, 
                    dimension=Counts,
                    dynamic=False).redim(**m)
            else:
                label = "{} {} points".format(Counts,
                                              len(self._data.embedding))
                #print "Alt", Counts
                r = histogram(
                    self._points[Counts], dimension="Value",
                    dynamic=False).redim(Value_frequency="frequency")

            #print(r)
            return r.relabel(label)

        from holoviews import streams
        self.zoom_selection = streams.Selection1D(source=self.tap_zoom)

        self.p = hv.DynamicMap(
            _h,
            kdims=["Counts"],
            streams=[
                self.zoom_selection.rename(index="fine_index"), self._posxy
            ]).redim.values(
                Counts=self._data.data_dims).opts(norm=dict(framewise=True))


        self._layout = self.ds_points * self._hover_grid * self.tap_indicators \
                       + self.selected_table + self.tap_zoom + self.p

        self._layout = self._layout.cols(2).opts(plot={"shared_axes": False})

        return self._layout

    def tap_table(self, left_x, bottom_y, right_x, top_y, index):
        e = self._posxy._dataset[self._data.data_dims]
        #print left_x, bottom_y, right_x, top_y, index
        if len(index) == 0 or left_x is None or bottom_y is None:
            #return hv.Points(([0],[0]))
            return hv.Table(e.head(0).reset_index())

        d = e.iloc[index]
        if len(index) < 1000:
            return hv.Table(d.reset_index(), label="%d sequences" % (len(d)))
        else:
            p = hv.Table(
                e.head(0).reset_index(),
                label="Not showing %d sequences" % (len(d)))
        return p

    def tap_points(self, left_x, bottom_y, right_x, top_y, index):
        #def tap_points(self, **kwargs):
        if left_x is None:
            return hv.Points(([], []))
            #return hv.Table(e.head(0).reset_index())

        p = hv.Points(
            ([left_x, left_x, right_x, right_x],
             [bottom_y, top_y, top_y, bottom_y])).opts(style=dict(color="red"))

        return p

    def focus_plot(self, left_x, bottom_y, right_x, top_y, index, Counts):
        from holoviews import streams
        import holoviews as hv
        e = self._posxy._dataset
        d = []
        if len(index) == 0 or left_x is None or bottom_y is None or len(
                index) > 100000:
            p = hv.Points(
                e.reset_index().head(0),
                kdims=self._data.coord_dims,
                vdims=["Sequence"] + self._data.data_dims,
                label="%d sequences" % (len(e)))

        else:
            d = e.iloc[index]
            p = hv.Points(
                d.reset_index(),
                kdims=self._data.coord_dims,
                vdims=["Sequence"] + self._data.data_dims,
                label="%d sequences" % (len(d)))

        self.zoom_selection.reset()
        p = p.opts(plot=dict(
            color_index=Counts,
            tools=["hover", "lasso_select", "box_select"],
            width=400,
            height=400))
        #return p.hist(dimension=Counts, adjoin=False)
        #print(p.dframe().head())
        return p
        p = p.hist(dimension="MeanFold", adjoin=False)
        print(p)
        return p

    def _get_selected(self):
        self._selected = self._data.embedding.head(0)
        if self._posxy.selected is None or len(self._posxy.selected) > 0:
            self._selected = self._posxy._dataset.iloc[self._posxy.selected]
            if len(self.zoom_selection.index) > 0:
                self._selected = self._selected.iloc[self.zoom_selection.index]
        else:
            self._selected = self._data.embedding

        return self._selected

    selected = property(_get_selected)


# Define valid function for FunctionHandler
# when deploying as script, simply attach to curdoc
def modify_doc2(doc, data):
    renderer = hv.renderer('bokeh').instance(mode='server')
    stream = hv.streams.Stream.define('Phase', phase=0.)()
    dmap = hv.DynamicMap(sine, streams=[stream])

    # Create HoloViews plot and attach the document
    hvplot = renderer.get_plot(dmap, doc)

    # Create a slider and play buttons
    def animate_update():
        year = slider.value + 0.2
        if year > end:
            year = start
        slider.value = year

    def slider_update(attrname, old, new):
        # Notify the HoloViews stream of the slider update 
        stream.event(phase=new)

    start, end = 0, np.pi * 2
    slider = Slider(start=start, end=end, value=start, step=0.2, title="Phase")
    slider.on_change('value', slider_update)

    def animate():
        if button.label == '► Play':
            button.label = '❚❚ Pause'
            doc.add_periodic_callback(animate_update, 50)
        else:
            button.label = '► Play'
            doc.remove_periodic_callback(animate_update)

    button = Button(label='► Play', width=60)
    button.on_click(animate)

    p = data.holoview_plot()
    hvplot_p = renderer.get_plot(p, doc)

    # Combine the holoviews plot and widgets in a layout
    plot = layout(
        [[hvplot.state], [slider, button], hvplot_p.state],
        sizing_mode='fixed')
    plot = plot.opts(plot={"shared_axes": False})
    #plot = layout([p.state], sizing_mode='fixed')

    doc.add_root(plot)
    return doc


def modify_doc(doc, data):
    from holoviews.operation.datashader import aggregate, shade, datashade, dynspread
    import holoviews as hv

    renderer = hv.renderer('bokeh').instance(mode='server')

    p = data.holoview_plot()

    p = p.source.table()
    hvplot_p = renderer.get_plot(p, doc)

    # Combine the holoviews plot and widgets in a layout
    plot = layout([hvplot_p.state], sizing_mode='fixed')
    #plot = layout([p.state], sizing_mode='fixed')

    doc.add_root(plot)
    return doc


def modify_doc(doc, data):
    from holoviews.operation.datashader import aggregate, shade, datashade, dynspread
    import holoviews as hv

    renderer = hv.renderer('bokeh').instance(mode='server')

    #p = points + hv.DynamicMap(selected_info, streams=[selection])
    import pdb
    pdb.set_trace()
    if isinstance(data, TsneMapper):

        hb = HuddingBrowser(data)
        p = hb.holoview_plot()
    else:
        p = data.plot_polar("HNF4A")

    hvplot_p = renderer.get_plot(p, doc)

    # Combine the holoviews plot and widgets in a layout
    plot = layout([hvplot_p.state], sizing_mode='fixed')
    #plot = layout([p.state], sizing_mode='fixed')

    doc.add_root(plot)
    return doc


def huddinge_app(data):

    # To display in the notebook
    #handler = FunctionHandler(modify_doc)
    #app = Application(handler)
    #show(app, notebook_url='localhost:8888')

    # To display in a script
    doc = modify_doc(curdoc(), data)
    return doc


def huddinge_app(data):
    import logging as log
    import holoviews as hv

    #p = points + hv.DynamicMap(selected_info, streams=[selection])
    hb = HuddingBrowser(data)
    p = hb.holoview_plot()

    doc = hv.renderer('bokeh').server_doc(p)
    doc.title = "Huddinge Browser"

    log.info("Returning doc")

    return doc
