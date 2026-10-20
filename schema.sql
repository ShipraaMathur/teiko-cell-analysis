-- Relational schema for the Loblaw Bio immune cell count trial data.
--
-- Design notes
--   * Attributes that never vary within a subject (project, condition, age,
--     sex, treatment, response) are stored once on the subject row.
--   * Attributes that vary per biological sample (sample type, time point)
--     live on the sample row.
--   * Cell counts are stored in long format (one row per sample and
--     population). Adding a new population later only requires a new row in
--     the populations table, not a schema migration.

PRAGMA foreign_keys = ON;

DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS populations;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS subjects;
DROP TABLE IF EXISTS projects;

CREATE TABLE projects (
    project_id   TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id   TEXT PRIMARY KEY,
    project_id   TEXT NOT NULL REFERENCES projects(project_id),
    condition    TEXT NOT NULL,              -- indication, e.g. melanoma
    age          INTEGER,
    sex          TEXT CHECK (sex IN ('M', 'F')),
    treatment    TEXT NOT NULL,              -- e.g. miraclib, phauximab, none
    response     TEXT CHECK (response IN ('yes', 'no') OR response IS NULL)
);

CREATE TABLE samples (
    sample_id                  TEXT PRIMARY KEY,
    subject_id                 TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type                TEXT NOT NULL,   -- e.g. PBMC, WB
    time_from_treatment_start  INTEGER NOT NULL
);

CREATE TABLE populations (
    population_id  INTEGER PRIMARY KEY,
    name           TEXT NOT NULL UNIQUE      -- e.g. b_cell, cd8_t_cell
);

CREATE TABLE cell_counts (
    sample_id      TEXT    NOT NULL REFERENCES samples(sample_id),
    population_id  INTEGER NOT NULL REFERENCES populations(population_id),
    count          INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population_id)
);

-- Indexes supporting the filters used in the analysis.
CREATE INDEX idx_subjects_filter ON subjects (condition, treatment, response, sex);
CREATE INDEX idx_subjects_project ON subjects (project_id);
CREATE INDEX idx_samples_subject ON samples (subject_id);
CREATE INDEX idx_samples_filter ON samples (sample_type, time_from_treatment_start);
CREATE INDEX idx_counts_population ON cell_counts (population_id);
