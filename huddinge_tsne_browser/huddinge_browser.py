# -*- coding: utf-8 -*-
import numpy as np
import holoviews as hv

from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application
from bokeh.io import curdoc
from bokeh.layouts import layout
from bokeh.models import Slider, Button

from datashaderselect import DataShaderSelect


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
                kdims=["tsne0", "tsne1"],
                vdims=["Sequence"] + self._data.data_dims,
                label="Sequences", )

        else:
            overlay_dict = {}
            for d in self._data.data_dims:
                p = hv.Points(
                    self._data.embedding[["tsne0", "tsne1", d]].rename(
                        columns={d: "Value"}).reset_index(),
                    kdims=["tsne0", "tsne1"],
                    vdims=["Sequence", "Value"],
                    label="Sequences", )
                overlay_dict[d] = p
            #self._points = hv.NdOverlay(overlay_dict, kdims=["Counts"])
            self._points = hv.HoloMap(overlay_dict, kdims=["Counts"])

    def datashade(self, column=None, reduction=None):
        """Return datashader shaded with reduction of given column value
        
        Arguments:
        - `column`:
        - `reduction`:
        """
        import datashader as ds
        from holoviews.operation.datashader import aggregate, shade, datashade, dynspread
        if reduction is None:
            reduction = ds.mean

        print "datashade", column, reduction, self._points
        plot = dynspread(
            datashade(
                self._points,
                aggregator=ds.count()
                if column is None else reduction(column)))
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
        #        else:
        #            agg_dict = {k: ds.mean(k) for k in self._data.data_dims}
        #            agg = ds.summary(**agg_dict)
        #            print agg
        #            print agg_dict
        #agg = ds.mean(self._data.data_dims[0])

        #        self.ds_points = dynspread(
        #            datashade(self._points, aggregator=agg)).opts(
        #                plot=dict(width=600, height=600))
        self.ds_points = self.ds_points.opts(plot=dict(width=600, height=600))

        # Hover and zoom grid tool.
        self._hover_grid = hv.util.Dynamic(
            aggregate(
                self._points,
                aggregator=ds.count(),
                width=15,
                height=15,
                streams=[RangeXY(source=self.ds_points)]),
            operation=hv.QuadMesh).opts(
                plot=dict(tools=["hover"]),
                style=dict(alpha=0, hover_alpha=0.2))

        # Get the points in tapped rectangle
        self._posxy = DataShaderSelect(
            source=self._hover_grid, dataset=self._data.embedding)

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
            self.focus_plot,
            streams=[self._posxy]).opts(norm=dict(framewise=True))

        self._layout = self.ds_points * self._hover_grid * self.tap_indicators + self.selected_table + self.tap_zoom
        #self._layout = ds_points * self._hover_grid
        self._layout = self._layout.cols(2).opts(plot={"shared_axes": False})

        return self._layout

    def tap_table(self, left_x, bottom_y, right_x, top_y, index):
        e = self._posxy._dataset
        print index
        if len(index) == 0 or left_x is None or bottom_y is None:
            #return hv.Points(([0],[0]))
            return hv.Table(e.head(0).reset_index())

        d = e.iloc[index]
        if len(index) < 1000:
            return hv.Table(d.reset_index(), label="%d sequences" % (len(d)))
        else:
            p = hv.Table(
                e.head(0).reset_index()[["Sequence"] + self._data.data_dims],
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

    def focus_plot(self, left_x, bottom_y, right_x, top_y, index):
        e = self._posxy._dataset
        if len(index) == 0 or left_x is None or bottom_y is None or len(
                index) > 1000:
            p = hv.Points(
                e.reset_index().head(0),
                kdims=["tsne0", "tsne1"],
                vdims=["Sequence"] + self._data.data_dims,
                label="%d sequences" % (len(e)))
        else:
            d = e.iloc[index]
            p = hv.Points(
                d.reset_index(),
                kdims=["tsne0", "tsne1"],
                vdims=["Sequence"] + self._data.data_dims,
                label="%d sequences" % (len(d)))

        return p.opts(plot=dict(
            tools=["hover", "lasso_select", "box_select"],
            width=400,
            height=400))

    def _get_selected(self):
        self._selected = self._data.embedding.head(0)
        if len(self._posxy.selected) > 0:
            self._selected = self._posxy._dataset.iloc[self._posxy.selected]

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
    hb = HuddingBrowser(data)
    p = hb.holoview_plot()

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
    import holoviews as hv

    #p = points + hv.DynamicMap(selected_info, streams=[selection])
    hb = HuddingBrowser(data)
    p = hb.holoview_plot()

    doc = hv.renderer('bokeh').server_doc(p)

    return doc
