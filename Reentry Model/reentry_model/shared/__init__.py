"""Code shared by the model's two architectures, skyscrapers and the level set (design of record: the level-set
design's section 4.2 under docs/superpowers/specs/levelset/).

Modules here import nothing from reentry_model outside this package except the package-root constants __version__,
PROJECT_ROOT, REPO_ROOT, DATA_DIR, OUTPUT_DIR, DRAMA_DIR and DRAMA_DATA_DIR (tests/test_reentry_model_shared.py checks
it). The old module paths (reentry_model.aero, ...) are aliases of these modules. A change here must keep both
architectures' tests passing."""
