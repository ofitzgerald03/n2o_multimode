#!/usr/bin/env python3
"""Run the analysis tasks from the project root:"""

#     python run.py --list                            # list the tasks
#     python run.py all                               # every result task
#     python run.py thruster_figures                  # one task
#     python run.py check_thruster check_spacecraft   # the two check suites

import argparse
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run the thesis analysis tasks.",
        epilog="Use --list to see the available tasks.",
    )
    parser.add_argument("tasks", nargs="*", help="task name(s), or 'all'")
    parser.add_argument(
        "--list", action="store_true", help="list the available tasks and exit"
    )
    parser.add_argument(
        "--no-save", action="store_true", help="do not write CSV and JSON outputs"
    )
    args = parser.parse_args(argv)

    # Backend before pyplot is imported anywhere
    import matplotlib

    matplotlib.use("Agg")

    from multimode import analyses  # noqa: F401
    from multimode.analyses.common import make_context
    from multimode.analyses.registry import TASKS, get, all_task_names

    if args.list or not args.tasks:
        width = max(len(n) for n in TASKS) + 2
        print("Available tasks (in the order of `all`):\n")
        for name, t in TASKS.items():
            flag_str = "" if t["in_all"] else "  [run by name only]"
            print(f"  {name:<{width}} {t['help']}{flag_str}")
        print("\nUsage: python run.py <task> [<task> ...] | all [--no-save]")
        return 0

    if args.tasks == ["all"]:
        names = all_task_names()
    else:
        names = args.tasks
        for n in names:
            get(n)

    ctx = make_context(save_csvs=not args.no_save)

    n_total = len(names)
    for k, name in enumerate(names, start=1):
        print(f"\n[{k}/{n_total}] {name}: {TASKS[name]['help']}\n")
        TASKS[name]["func"](ctx)

    print(f"\nDone: {n_total} task(s) complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
