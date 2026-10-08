"""python -m spheral_frag {prepare,analyse} ... (spec §4.2; M1 plan, Task 9).

  prepare   finite-element run directory -> Spheral input files    (spheral_frag/prepare.py)
  analyse   Spheral run outputs -> fragment record, debris log, mass accounts    (spheral_frag/analyse.py)

Exit codes as in the finite-element model: 0 ok, 1 a model failure (a contract or check failure, accounts that do not
close; message printed), 2 bad arguments, a missing input or a missing optional library."""
import argparse
import sys

USAGE = ("usage: python -m spheral_frag {prepare,analyse} [options]\n"
         "  prepare   finite-element run directory -> Spheral input files\n"
         "  analyse   Spheral run outputs -> fragment record, debris log, mass accounts\n")


def build_parser():
    from . import analyse, prepare
    p = argparse.ArgumentParser(prog="python -m spheral_frag", description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", metavar="{prepare,analyse}")
    prepare.add_arguments(sub.add_parser("prepare", help="finite-element run directory -> Spheral input files"))
    analyse.add_arguments(sub.add_parser("analyse", help="Spheral run outputs -> fragment record, debris log, "
                                                         "mass accounts"))
    return p


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] not in ("prepare", "analyse", "-h", "--help"):
        sys.stderr.write(USAGE)
        return 2
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as e:                          # argparse: 2 on bad arguments, 0 on --help
        return int(e.code or 0)
    args.argv = argv
    if args.command == "prepare":
        from .prepare import prepare
        return prepare(args)
    from .analyse import analyse
    return analyse(args)


if __name__ == "__main__":
    sys.exit(main())
