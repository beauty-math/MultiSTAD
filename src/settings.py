from pathlib import Path
import os

PACKAGE = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("MULTISTAD_DATA_ROOT", PACKAGE / "data")).resolve()
OUTPUT = Path(os.environ.get("MULTISTAD_OUTPUT_ROOT", PACKAGE / "results")).resolve()
DEFINITIONS = PACKAGE / "data/definitions"
