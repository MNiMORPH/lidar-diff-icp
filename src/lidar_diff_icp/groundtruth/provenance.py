"""Parameters that carry their own origin, so the label on a number is generated.

``trust/provenance.py`` makes an analysis *script* declare where each parameter came
from. That only works if the library hands the script something to declare. Every
estimator in this package therefore returns its settings as :class:`Param` records --
value plus ``src`` (``"andy"`` / ``"repo"`` / ``"MINE"``) plus, for a repo default, the
file and line it was read from -- and :func:`declare` pushes them into a
:class:`trust.provenance.Run`.

The point is that a caller cannot print a tie number without also being able to print,
without typing them, the radius ladder, the ground source, the quantile and the geoid
grids that produced it.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Param:
    """One setting, with where it came from.

    ``src`` follows ``trust.provenance.Run.param``: ``"andy"`` (the human said so this
    session), ``"repo"`` (an established default elsewhere in this codebase, and ``why``
    names the file), ``"MINE"`` (chosen unasked -- ``why`` must say what it excludes).
    """

    name: str
    value: object
    src: str
    why: str = ""

    def __str__(self) -> str:  # pragma: no cover - display only
        tag = {"andy": "[asked]", "repo": "[repo default]",
               "MINE": "[** MINE, UNASKED **]"}.get(self.src, f"[{self.src}]")
        s = f"{self.name} = {self.value!r}  {tag}"
        return s + (f"\n      why: {self.why}" if self.why else "")


def declare(run, params, *, prefix: str = "") -> None:
    """Push ``params`` into a :class:`trust.provenance.Run` (``run.param``).

    ``prefix`` disambiguates when several estimators are declared into one run
    (e.g. one per checkpoint). Raises whatever ``Run.param`` raises -- an unexplained
    ``src="MINE"`` still fails there, which is the point.
    """
    for p in params:
        run.param(f"{prefix}{p.name}", p.value, src=p.src, why=p.why)
