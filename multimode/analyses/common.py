"""
Run context and console helpers shared by the analysis tasks.
"""

from types import SimpleNamespace

from ..config import PATHS, RUN_OPTIONS


def make_context(save_csvs=None):
    """Run context passed to every task; save_csvs=None takes RUN_OPTIONS."""
    return SimpleNamespace(
        save_csvs=RUN_OPTIONS["save_csvs"] if save_csvs is None else save_csvs,
        paths=PATHS,
    )


def csv_dir(ctx):
    d = ctx.paths["output"]
    d.mkdir(parents=True, exist_ok=True)
    return d


def section(title, char="="):
    bar = char * 78
    print(f"\n{bar}\n{title}\n{bar}")
