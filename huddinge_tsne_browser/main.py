# -*- coding: utf-8 -*-
import logging as log
from huddinge_tsne_browser import cli
from tsne_mapper import TsneMapper
from huddinge_browser import huddinge_app

import holoviews as hv
hv.extension('bokeh')

args = cli.cli()

data = TsneMapper(args.input)

if not data.laidout():
    data.compute_tsne()

try:
    data.write_data(args.output)
except TypeError:
    pass

log.info("Executing module %s", __name__)

for n, f in args.kmers:
    log.info("Adding %s as %s", f, n)
    data.add_kmercounts(n, f)

data.add_kmercounts(
    "HNF4A",
    "/home/kpalin/software/MODER/data/HNF4A_TGACAG20NGA_AF_1.mer_counts.jf")


def serve_embedding(data):
    from huddinge_browser import HuddingBrowser
    import holoviews as hv
    import traceback

    renderer = hv.renderer('bokeh')

    hb = HuddingBrowser(data)
    p = hb.holoview_plot()

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
    except Exception, e:
        traceback.print_exc()

        server.stop()


if __name__ == "__main__":
    if args.html is not None:
        log.info("Outputting %s", args.html)
        from bokeh.io import output_file, show, save

        #save(data.html(), args.html)
        serve_embedding(data)
else:
    log.info("serving plots")
    huddinge_app(data)
