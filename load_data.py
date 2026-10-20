"""Part 1: initialise the SQLite database and load cell-count.csv.

Usage:
    python load_data.py

Creates (or recreates) cell_counts.db in the repository root using the
schema defined in schema.sql, then loads every row of cell-count.csv.
"""

import csv
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
SCHEMA_PATH = ROOT / "schema.sql"
DB_PATH = ROOT / "cell_counts.db"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
SUBJECT_FIELDS = ["project", "condition", "age", "sex", "treatment", "response"]


def read_rows(csv_path: Path) -> list[dict]:
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = set(["subject", "sample", "sample_type",
                       "time_from_treatment_start", *SUBJECT_FIELDS,
                       *POPULATIONS]) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"cell-count.csv is missing columns: {sorted(missing)}")
        return list(reader)


def build_records(rows: list[dict]):
    """Split flat CSV rows into normalised records, validating consistency."""
    projects: set[str] = set()
    subjects: dict[str, tuple] = {}
    samples: list[tuple] = []
    counts: list[tuple] = []
    seen_samples: set[str] = set()

    for line_no, row in enumerate(rows, start=2):
        subject_id = row["subject"].strip()
        sample_id = row["sample"].strip()
        response = row["response"].strip() or None

        subject_record = (
            subject_id,
            row["project"].strip(),
            row["condition"].strip(),
            int(row["age"]) if row["age"].strip() else None,
            row["sex"].strip(),
            row["treatment"].strip(),
            response,
        )
        previous = subjects.get(subject_id)
        if previous is not None and previous != subject_record:
            raise ValueError(
                f"Line {line_no}: subject {subject_id} has inconsistent "
                f"metadata ({previous} vs {subject_record})"
            )
        subjects[subject_id] = subject_record
        projects.add(subject_record[1])

        if sample_id in seen_samples:
            raise ValueError(f"Line {line_no}: duplicate sample id {sample_id}")
        seen_samples.add(sample_id)
        samples.append((
            sample_id,
            subject_id,
            row["sample_type"].strip(),
            int(row["time_from_treatment_start"]),
        ))

        for population_id, population in enumerate(POPULATIONS, start=1):
            counts.append((sample_id, population_id, int(row[population])))

    return projects, subjects, samples, counts


def main() -> int:
    if not CSV_PATH.exists():
        print(f"Error: {CSV_PATH.name} not found in {ROOT}", file=sys.stderr)
        return 1

    rows = read_rows(CSV_PATH)
    projects, subjects, samples, counts = build_records(rows)

    if DB_PATH.exists():
        DB_PATH.unlink()

    conn = sqlite3.connect(DB_PATH)
    try:
        conn.executescript(SCHEMA_PATH.read_text(encoding="utf-8"))
        with conn:
            conn.executemany(
                "INSERT INTO projects (project_id) VALUES (?)",
                [(p,) for p in sorted(projects)],
            )
            conn.executemany(
                "INSERT INTO subjects (subject_id, project_id, condition, age, "
                "sex, treatment, response) VALUES (?, ?, ?, ?, ?, ?, ?)",
                list(subjects.values()),
            )
            conn.executemany(
                "INSERT INTO samples (sample_id, subject_id, sample_type, "
                "time_from_treatment_start) VALUES (?, ?, ?, ?)",
                samples,
            )
            conn.executemany(
                "INSERT INTO populations (population_id, name) VALUES (?, ?)",
                list(enumerate(POPULATIONS, start=1)),
            )
            conn.executemany(
                "INSERT INTO cell_counts (sample_id, population_id, count) "
                "VALUES (?, ?, ?)",
                counts,
            )
    finally:
        conn.close()

    print(f"Created {DB_PATH.name}")
    print(f"  projects:    {len(projects)}")
    print(f"  subjects:    {len(subjects)}")
    print(f"  samples:     {len(samples)}")
    print(f"  cell counts: {len(counts)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
