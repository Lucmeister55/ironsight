#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import argparse
import sys

# Local imports
import ironsight as pkg
from ironsight.common import *
from ironsight.Predict_Ferro import Predict_Ferro

# ~~~~~~~~~~~~~~TOP LEVEL ENTRY POINT~~~~~~~~~~~~~~#
def main(args=None):
    """ Main entry point for ironsight command line interface"""
    # Parser and subparsers for command
    parser = argparse.ArgumentParser(description=pkg.__description__)
    parser.add_argument("--version", action="version", version="{} v{}".format(pkg.__name__, pkg.__version__))
    subparsers = parser.add_subparsers(description="%(prog)s implements the following subcommands", dest="subcommands")
    subparsers.required = True

    # Predict_Ferro subparser
    f = Predict_Ferro
    sp_met = subparsers.add_parser("Predict_Ferro", description=doc_func(f))
    sp_met.set_defaults(func=f)
    sp_met_io = sp_met.add_argument_group("Input/Output options")
    arg_from_docstr(sp_met_io, f, "modbam", "i")
    arg_from_docstr(sp_met_io, f, "outdir", "o")
    sp_met_ms = sp_met.add_argument_group("Misc options")
    arg_from_docstr(sp_met_ms, f, "max_missing", "m")

    # Add common group parsers
    for sp in [sp_met, sp_cr, sp_cgi]:
        sp_vb = sp.add_argument_group("Verbosity options")
        sp_vb.add_argument("-v", "--verbose", action="store_true", default=False, help="Increase verbosity")
        sp_vb.add_argument("-q", "--quiet", action="store_true", default=False, help="Reduce verbosity")
        sp_vb.add_argument("-p", "--progress", action="store_true", default=False, help="Display a progress bar")
    
    # Parse args and call subfunction
    args = parser.parse_args()
    
    args.func(**vars(args))


if __name__ == "__main__":
    main()