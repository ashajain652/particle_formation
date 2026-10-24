"""Assemble docs/superpowers/plans/2026-09-20-melt-spraying.md from the verified prototype files (proto3)."""
import os, sys

PROTO = os.environ.get("PROTO", os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "proto3"))
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "plan.md")


def read(path):
    full = os.path.join(PROTO, path)
    if not os.path.exists(full):
        return "@@@MISSING::%s@@@\n" % path
    with open(full) as fh:
        return fh.read().rstrip("\n") + "\n"


def block(text, lang="python"):
    return "```" + lang + "\n" + text + "```\n"


def tail_from(path, marker):
    if not os.path.exists(os.path.join(PROTO, path)):
        return "@@@MISSINGTAIL::%s@@@\n" % path
    text = read(path)
    if marker not in text:
        return "@@@MISSINGTAIL::%s@@@\n" % path
    i = text.index(marker)
    return text[i:]


def file_block(path, lang="python"):
    return block(read(path), lang)


parts = []
def emit(s):
    parts.append(s.rstrip("\n") + "\n")

exec(open(os.path.join(HERE, "part_header.py")).read())
exec(open(os.path.join(HERE, "part_tasks_a.py")).read())
exec(open(os.path.join(HERE, "part_tasks_b.py")).read())
exec(open(os.path.join(HERE, "part_tasks_c.py")).read())
exec(open(os.path.join(HERE, "part_tasks_c2.py")).read())
exec(open(os.path.join(HERE, "part_tasks_c3.py")).read())
exec(open(os.path.join(HERE, "part_tasks_d.py")).read())

with open(OUT, "w") as fh:
    fh.write("\n".join(parts))
print("wrote", OUT, sum(p.count("\n") for p in parts) + len(parts), "lines")
