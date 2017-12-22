# -*- coding: utf-8 -*-
import logging as log
from huddinge_tsne_browser import cli
from tsne_mapper import TsneMapper
from huddinge_browser import huddinge_app

args = cli.cli()

data = TsneMapper(args.input)

if not data.laidout():
    data.compute_tsne()

try:
    data.write_data(args.output)
except TypeError:
    pass

log.info("Executing module %s", __name__)

if __name__ == "__main__":
    if args.html is not None:
        log.info("Outputting %s", args.html)
        from bokeh.io import output_file, show, save

        save(data.html(), args.html)
else:
    log.info("serving plots")
    huddinge_app(data)
