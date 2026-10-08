"""`python -m reentry_model` runs the `run`/`compare` command line (see reentry_model.cli)."""
import sys

from .cli import main

sys.exit(main())
