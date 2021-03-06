#!/usr/bin/env python
#-*- coding: utf-8 -*-
""" 

Created at Mon Dec 18 15:15:47 2017 by Kimmo Palin <kpalin@merit.ltdk.helsinki.fi>
"""


def cli(args=None):
    import argparse
    parser = argparse.ArgumentParser(description="Launch tsne browser")

    parser.add_argument(
        "-i",
        "--input",
        help="Input file. Either with or without TSNE layout. Without version can be generated with all_pairs_huddinge  [default:%(default)s]",
        required=True)

    parser.add_argument(
        "-o",
        "--output",
        help="Store the input in this file along with the sequence placements")
    parser.add_argument(
        "-t",
        "--html",
        help="Output as html file  [default:%(default)s]",
        default="huddinge_tsne.html")

    parser.add_argument(
        "-k",
        "--kmers",
        help="Additional annotation kmer data. Currently only jellyfish output. Format name:/path/to/file.jf  [default:%(default)s]",
        default=[],
        action="append")

    parser.add_argument(
        "-V",
        "--verbose",
        default=False,
        const=True,
        nargs="?",
        help="Be more verbose with output [and log to a file] [default:%(default)s]"
    )

    import sys
    args = parser.parse_args(sys.argv[1:] if args is None else args)

    import logging
    if args.verbose:
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s:%(funcName)s:%(levelname)s:%(message)s')
        if args.verbose != True:
            log_file_handler = logging.FileHandler(args.verbose)
            log_file_handler.setFormatter(
                logging.getLogger().handlers[0].formatter)
            logging.getLogger().addHandler(log_file_handler)

    kmers_args = []
    from os.path import basename, exists
    for x in args.kmers:
        p = x.split(":")
        if len(p) == 1:
            p = (basename(p[0]).split(".")[0], p[0])
        elif len(p) == 2:
            p = tuple(p)
        else:
            raise ValueError(str(p))
        if not exists(p[1]):
            raise ValueError(p[1])
        kmers_args.append(p)
    args.kmers = kmers_args

    return args
