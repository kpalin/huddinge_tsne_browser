# -*- coding: utf-8 -*-
import logging as log
from huddinge_tsne_browser import cli
from .tsne_mapper import TsneMapper, PolarMapper
from .huddinge_browser import huddinge_app

import holoviews as hv
hv.extension('bokeh')

args = cli.cli()

data = PolarMapper(args.input, args.json_config)

log.info("Executing module %s", __name__)


def serve_embedding(data):
    from .huddinge_browser import HuddingBrowser
    import holoviews as hv
    import traceback

    renderer = hv.renderer('bokeh')

    #import pdb
    #pdb.set_trace()
    p = data.plot_polar("HNF4A")

    app = renderer.app(p)

    from tornado.ioloop import IOLoop
    from bokeh.server.server import Server

    loop = IOLoop.current()
    server = Server({'/': app}, port=0, io_loop=loop)

    def show_callback():
        server.show('/')

    loop.add_callback(show_callback)
    server.start()
    try:
        log.info("Starting loop")
        loop.start()
    except Exception as e:
        traceback.print_exc()

        server.stop()


if __name__ == "__main__":
    if args.html is not None:
        log.info("Outputting %s", args.html)
        from bokeh.io import output_file, show, save

        save(data.html(), args.html)
    else:
        serve_embedding(data)
else:
    log.info("serving plots")
    #huddinge_app(data)
    serve_embedding(data)
