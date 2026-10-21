from pathlib import Path
import re, json, pandas as pd
from .settings import DATA

D = DATA / "metadata"


def soft(ds):
    text = (D / (ds + "_official.soft.txt")).read_text(errors="replace")
    samples = {}
    current = None
    for line in text.splitlines():
        if line.startswith("^SAMPLE = "):
            current = line.split(" = ")[1]
            samples[current] = {}
        elif current and line.startswith("!Sample_") and " = " in line:
            k, v = line.split(" = ", 1)
            samples[current].setdefault(k, []).append(v)
    return text, samples
