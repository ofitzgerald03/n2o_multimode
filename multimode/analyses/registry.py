"""
Task registry. Each analysis registers itself with the @task decorator;
run.py looks tasks up here.
"""

from collections import OrderedDict

# name -> {'func', 'help', 'in_all'}
TASKS = OrderedDict()


def task(name, help_text, *, in_all=True):
    """Register the decorated function as task `name`. Tasks with
    in_all=False run only when named."""

    def _register(func):
        if name in TASKS:
            raise ValueError(f"Duplicate task name: {name!r}")
        TASKS[name] = {
            "func": func,
            "help": help_text,
            "in_all": in_all,
        }
        return func

    return _register


def get(name):
    if name not in TASKS:
        known = ", ".join(TASKS)
        raise KeyError(f"Unknown task {name!r}. Known tasks: {known}")
    return TASKS[name]


def all_task_names():
    """Names of the tasks run by `run.py all`, in registration order."""
    return [n for n, t in TASKS.items() if t["in_all"]]
