"""Database access helpers shared by the pipeline and the dashboard."""

import sqlite3
from contextlib import closing
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DB_PATH = ROOT / "cell_counts.db"
OUTPUT_DIR = ROOT / "outputs"


def connect(db_path: Path = DB_PATH) -> sqlite3.Connection:
    if not Path(db_path).exists():
        raise FileNotFoundError(
            f"Database {Path(db_path).name} not found. "
            "Run `python load_data.py` or `make pipeline` first."
        )
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def query(sql: str, params: tuple | dict = (), db_path: Path = DB_PATH) -> pd.DataFrame:
    with closing(connect(db_path)) as conn:
        return pd.read_sql_query(sql, conn, params=params)
