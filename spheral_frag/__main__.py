"""python -m spheral_frag {prepare,analyse} ... (the subcommands are wired in M1 Task 9; this is the stub).

Exit codes as in the finite-element model: 0 ok, 1 a model failure, 2 bad arguments or a missing optional library."""
import sys

USAGE = ("usage: python -m spheral_frag {prepare,analyse} [options]\n"
         "  prepare   finite-element run directory -> Spheral input files (not implemented yet)\n"
         "  analyse   Spheral run outputs -> fragment record, debris log, mass accounts (not implemented yet)\n")


def main(argv=None):
    sys.stderr.write(USAGE)
    return 2


if __name__ == "__main__":
    sys.exit(main())
