"""
Appendix B tables (the verification checks, grouped by test level) and the
chapter 5 summary table, from the logs of the two check suites:
"""

#     python run.py check_thruster --no-save > thruster.log
#     python run.py check_spacecraft --no-save > spacecraft.log
#     python scripts/make_appendix_b.py --thruster-log thruster.log \
#         --spacecraft-log spacecraft.log --out-dir output
#
# Statements, tolerances, layers and levels are in appendix_b_spec.py. IDs
# are T.n for the thruster suite and Sg.n for the spacecraft suite, g
# numbering its sections in order.

import argparse
import re
import sys
from collections import Counter, OrderedDict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from appendix_b_spec import SL, CHK, GROUPS, LEVEL, LEVELS  # noqa: E402

BS2 = chr(92) * 2
LAYERS = ("thruster", "system", "sweep")
INTRO = (
    r"This appendix lists every automated check in the verification suites of the two models, "
    r"grouped by test level. \Cref{tab:app_verif_levels} summarises the levels and the number of "
    r"checks in each. Each check compares a model output, within a stated tolerance, with an "
    r"independently computed value, a limiting case, a required property, or a frozen value from an "
    r"earlier run; regression tests against frozen values guard against unintended change without "
    r"themselves establishing correctness (\Cref{sec:sys_verification}). Some checks repeat a design "
    r"with a composite-overwrapped \ce{N2O} tank in place of the all-metal tanks used for every "
    r"reported result; they confirm that the model handles that tank type, which is not used "
    r"elsewhere in this thesis."
)
# chapter 5 wording where the spacecraft suite differs from both suites
CH5_WHAT = {
    "unit": "A single relation or input agrees with an independent hand calculation or its specification"
}


def checks(path, kind):
    """(section, label, status) of every check line in a suite log."""
    out, sec = [], None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if kind == "sl" and line and not line.startswith((" ", "=")):
            sec = line.split(" — ")[0].strip()
        m = re.match(r"^\s*\[(OK|FAIL)\s*\]\s+(.*)$", line)
        if m:
            out.append((sec, m.group(2).split(" — ")[0].strip(), m.group(1)))
    return out


def build(check_log, sl_log):
    sl, ck = checks(sl_log, "sl"), checks(check_log, "check")
    assert len(sl) == len(SL), (len(sl), len(SL))
    assert len(ck) == len(CHK), (len(ck), len(CHK))
    rows = []  # (id, sentence, tolerance, layer, status)
    for n, ((_s, lab, st), (pre, sent, tol, lay)) in enumerate(zip(ck, CHK), 1):
        assert lab.startswith(pre), (lab, pre)
        rows.append((f"T.{n}", sent, tol, lay, st))
    gnum, cnt = OrderedDict(), Counter()
    for (sec, lab, st), (pre, sent, tol, lay) in zip(sl, SL):
        assert lab.startswith(pre), (sec, lab, pre)
        assert sec in GROUPS, sec
        g = sec
        gnum.setdefault(g, len(gnum) + 1)
        cnt[g] += 1
        rows.append((f"S{gnum[g]}.{cnt[g]}", sent, tol, lay, st))
    ids = [r[0] for r in rows]
    assert set(ids) == set(LEVEL) and len(ids) == len(LEVEL), sorted(
        set(ids) ^ set(LEVEL)
    )
    return rows


def longtable(caption, label, entries):
    t = [
        r"\begingroup",
        r"\small",
        r"\begin{longtable}{@{}l >{\raggedright\arraybackslash}p{0.53\textwidth} "
        r">{\raggedright\arraybackslash}p{0.2\textwidth} l@{}}",
        r"\caption[" + caption + "]{" + caption + r"}\label{" + label + r"} " + BS2,
        r"\toprule",
        r"ID & Assertion & Tolerance & Layer " + BS2,
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"ID & Assertion & Tolerance & Layer " + BS2,
        r"\midrule",
        r"\endhead",
        r"\bottomrule",
        r"\endlastfoot",
    ]
    for i, sent, tol, lay, _st in entries:
        t.append(f"{i} & {sent} & {tol} & {lay} " + BS2)
    t += [r"\end{longtable}", r"\endgroup", ""]
    return t


def summary_table(caption_short, caption, label, cols, body, total):
    """cols: list of (header, width fraction, align); body: list of cell lists."""
    spec = "\n".join(
        f"        >{{{a}\\arraybackslash}}p{{\\dimexpr {w}\\textwidth-2\\tabcolsep\\relax}}"
        for _h, w, a in cols
    )
    t = [
        r"\begin{table}[htbp]",
        r"    \centering",
        r"    \small",
        r"    \caption[" + caption_short + "]{" + caption + "}",
        r"    \label{" + label + "}",
        r"    \begin{tabular}{@{}",
        spec,
        r"    @{}}",
        r"        \toprule",
        "        " + " & ".join(h for h, _w, _a in cols) + " " + BS2,
        r"        \midrule",
    ]
    t += ["        " + " & ".join(c) + " " + BS2 for c in body]
    t += [
        r"        \midrule",
        "        " + " & ".join(total) + " " + BS2,
        r"        \bottomrule",
        r"    \end{tabular}",
        r"\end{table}",
        "",
    ]
    return t


def lint(tex):
    """LaTeX checks on the rendered text: no bare _ or ~, cell counts."""
    bad = []
    for ln in tex.splitlines():
        if ln.lstrip().startswith("%"):
            continue
        body = re.sub(r"\$[^$]*\$", "", ln)
        body = re.sub(r"\\(?:Cref|cref|ref|label)\{[^}]*\}", "", body)
        if "_" in body or "~" in body:
            bad.append(("_ or ~", ln[:80]))
        cells = body.count("&")
        if cells and cells not in (2, 3, 5):
            bad.append(("& count", ln[:80]))
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--thruster-log", required=True)
    ap.add_argument("--spacecraft-log", required=True)
    ap.add_argument("--out-dir", default=".")
    a = ap.parse_args()
    rows = build(a.thruster_log, a.spacecraft_log)
    out_dir = Path(a.out_dir)

    by_level = OrderedDict(
        (k, [r for r in rows if LEVEL[r[0]] == k]) for k, *_ in LEVELS
    )
    n_all = {k: len(v) for k, v in by_level.items()}
    n_sys = {k: sum(1 for r in v if r[0].startswith("S")) for k, v in by_level.items()}
    n_thr = {k: n_all[k] - n_sys[k] for k in n_all}
    lay = {
        k: [(x, sum(1 for r in v if r[3] == x)) for x in LAYERS]
        for k, v in by_level.items()
    }

    # Appendix B
    out = [
        "% Generated by scripts/make_appendix_b.py: Appendix B, the verification checks",
        f"% grouped by test level. Thruster suite ({len(CHK)}) and spacecraft suite ({len(SL)}).",
        r"\chapter{Verification Assertions}\label{app:verification_assertions}",
        INTRO,
        "",
    ]
    cnt = lambda n: rf"\num{{{n}}}" if n else "--"  # noqa: E731
    out += summary_table(
        "Test levels of the verification suites",
        "Test levels of the two verification suites, with what each establishes and the number of "
        "assertions at each model layer.",
        "tab:app_verif_levels",
        [
            ("Level", 0.14, r"\raggedright"),
            ("What it establishes", 0.42, r"\raggedright"),
        ]
        + [(x.capitalize(), 0.11, r"\centering") for x in LAYERS]
        + [("Total", 0.11, r"\centering")],
        [
            [title.replace(" tests", ""), what]
            + [cnt(n) for _x, n in lay[k]]
            + [cnt(n_all[k])]
            for k, title, what in LEVELS
        ],
        ["Total", ""]
        + [cnt(sum(dict(lay[k])[x] for k in lay)) for x in LAYERS]
        + [cnt(sum(n_all.values()))],
    )
    for k, title, _what in LEVELS:
        out += longtable(title, f"tab:app_verif_{k}", by_level[k])
    tex = "\n".join(out)

    # chapter 5 summary table
    ch5 = [
        "% Generated by scripts/make_appendix_b.py: the spacecraft suite by test level."
    ]
    ch5 += summary_table(
        "Test levels of the system-level verification suite",
        "Test levels of the system-level verification suite, with what each establishes and the "
        r"number of assertions at each level. The full list, with the thruster-level suite, is "
        r"given in Appendix \ref{app:verification_assertions}.",  # \Cref renders "Chapter B" here
        "tab:sys_gates",
        [
            ("Level", 0.2, r"\raggedright"),
            ("What it establishes", 0.62, r"\raggedright"),
            ("Assertions", 0.18, r"\centering"),
        ],
        [
            [title.replace(" tests", ""), CH5_WHAT.get(k, what), rf"\num{{{n_sys[k]}}}"]
            for k, title, what in LEVELS
        ],
        ["Total", "", rf"\num{{{sum(n_sys.values())}}}"],
    )
    ch5_tex = "\n".join(ch5)

    for name, t in (
        ("app_verification_assertions.tex", tex),
        ("tab_sys_test_levels.tex", ch5_tex),
    ):
        (out_dir / name).write_bytes((t + "\n").replace("\n", "\r\n").encode("utf-8"))

    # report
    names = {k: t.replace(" tests", "").lower() for k, t, _w in LEVELS}

    def seq(d):
        p = [f"{d[k]} {names[k]}" for k in names if d[k]]
        return ", ".join(p[:-1]) + " and " + p[-1]

    print("failing:", sum(1 for r in rows if r[4] != "OK"))
    print("all     ", n_all, sum(n_all.values()))
    print("system  ", n_sys, sum(n_sys.values()))
    print("thruster", n_thr, sum(n_thr.values()))
    print("layers  ", lay)
    print(
        "sentence: The "
        + str(sum(n_sys.values()))
        + " system-level assertions comprise "
        + seq(n_sys)
        + " tests; with the "
        + str(sum(n_thr.values()))
        + " thruster-level assertions ("
        + seq(n_thr)
        + ") there are "
        + str(sum(n_all.values()))
        + " in all."
    )
    print("lint:", lint(tex) + lint(ch5_tex) or "clean")
    # a sentence that leans on the rows above must keep its antecedent
    # chain (the back-referring rows before it and the row that heads them)
    # directly above it in its level table
    anaph = re.compile(r"(That |Its |The same |Every point of that )")
    stranded = []
    for v in by_level.values():
        for prev, cur in zip([None] + v[:-1], v):
            if not anaph.match(cur[1]):
                continue
            j, chain = rows.index(cur) - 1, set()
            while j >= 0:
                chain.add(rows[j][0])
                if not anaph.match(rows[j][1]):
                    break
                j -= 1
            if prev is None or prev[0] not in chain:
                stranded.append((cur[0], sorted(chain)))
    print("stranded back-references:", stranded or "none")


if __name__ == "__main__":
    main()
