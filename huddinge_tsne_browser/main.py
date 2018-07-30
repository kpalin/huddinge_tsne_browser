# -*- coding: utf-8 -*-
import logging as log
from huddinge_tsne_browser import cli
from .tsne_mapper import TsneMapper
from .polarmapper import PolarMapper
from .huddinge_browser import huddinge_app

import holoviews as hv

hv.extension('bokeh')


def main(args=None):
    
    args = cli.cli(args)
    log.info("Running in main()")

    log.info(str(args))
    data = PolarMapper(args.input,args.enrichment_column)
    
    if args.json_config is not None:
        data.read_config(args.json_config)

    log.info("Executing module %s", __name__)
    if args.html is not None:
        log.info("Outputting %s", args.html)
        data.save_bokeh_points(args.html)
    else:
        log.info("Serving plots from main")
        serve_embedding(data)


def serve_embedding(data):
    
    import holoviews as hv
    import traceback

    renderer = hv.renderer('bokeh')

    #import pdb
    #pdb.set_trace()
    p = data.plot_polar()

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
    main()

