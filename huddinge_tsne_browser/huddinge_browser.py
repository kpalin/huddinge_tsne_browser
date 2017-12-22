# -*- coding: utf-8 -*-
import numpy as np
import holoviews as hv

from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application
from bokeh.io import curdoc
from bokeh.layouts import layout
from bokeh.models import Slider, Button


# Create the holoviews app again
def sine(phase):
    xs = np.linspace(0, np.pi * 4)
    return hv.Curve((xs, np.sin(xs + phase))).opts(plot=dict(width=800))


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
    #plot = layout([p.state], sizing_mode='fixed')

    doc.add_root(plot)
    return doc


def modify_doc(doc, data):
    from holoviews.operation.datashader import aggregate, shade, datashade, dynspread
    import holoviews as hv
    from holoviews.streams import RangeXY

    renderer = hv.renderer('bokeh').instance(mode='server')

    p = data.holoview_plot()
    #mesh = hv.util.Dynamic(aggregate(points, width=10, height=10, streams=[RangeXY]), operation=hv.QuadMesh).opts(plot=dict(tools=["hover"]),style=dict(alpha=0,hover_alpha=0.2))
    #p = dynspread(datashade(p)) + p.table()

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
    from holoviews.streams import RangeXY

    renderer = hv.renderer('bokeh').instance(mode='server')

    points = hv.Points(
        data.embedding.reset_index(),
        kdims=["tsne0", "tsne1"],
        vdims=["Sequence"],
        label="Sequences", ).opts(plot=dict(
            tools=["lasso_select", "box_select"], width=800, height=800))
    #

    #points = hv.Points(np.random.randn(1000,2 )).opts(plot=dict(tools=['box_select', 'lasso_select']))
    # Declare points as source of selection stream
    selection = hv.streams.Selection1D(source=points)

    # Write function that uses the selection indices to slice points and compute stats
    def selected_info(index):
        #print index
        selected = points.iloc[index]
        if index:
            #raise str(index)
            label = 'Selected:'
        else:
            label = 'No selection'
        return hv.Table(selected).relabel(
            label)  #.opts(style=dict(color='red'))

    p = points + hv.DynamicMap(selected_info, streams=[selection])

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
