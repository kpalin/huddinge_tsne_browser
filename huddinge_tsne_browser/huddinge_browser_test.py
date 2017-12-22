import numpy as np
import holoviews as hv

from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application
from bokeh.io import curdoc
from bokeh.layouts import layout
from bokeh.models import Slider, Button
from holoviews.operation.datashader import aggregate, shade, datashade, dynspread

renderer = hv.renderer('bokeh').instance(mode='server')


# Create the holoviews app again
def sine(phase):
    xs = np.linspace(0, np.pi * 4)
    return hv.Curve((xs, np.sin(xs + phase))).opts(plot=dict(width=800))


stream = hv.streams.Stream.define('Phase', phase=0.)()
dmap = hv.DynamicMap(sine, streams=[stream])


# Define valid function for FunctionHandler
# when deploying as script, simply attach to curdoc
def modify_doc(doc):
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

    # Combine the holoviews plot and widgets in a layout
    plot = layout([[hvplot.state], [slider, button]], sizing_mode='fixed')

    doc.add_root(plot)
    return doc


# To display in the notebook
#handler = FunctionHandler(modify_doc)
#app = Application(handler)
#show(app, notebook_url='localhost:8888')

# To display in a script
doc = modify_doc(curdoc())
