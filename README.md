
# Take home assignment - Shipra Anil Mathur 


## Quick start (GitHub Codespaces)

```bash
make setup      # install dependencies from requirements.txt
make pipeline   # build cell_counts.db and write all results to outputs/
make dashboard  # start the dashboard on port 8501
```



## Dashboard

* **Local / Codespaces:** http://localhost:8501 (started by `make dashboard`)
* **Hosted version:** https://shipra-teiko-cell-analysis.streamlit.app/

The dashboard reads directly from `cell_counts.db` and has one tab per analysis:

| Tab | Contents |
|-----|----------|
| Part 2: Data Overview | Searchable and filterable relative frequency table, with CSV download |
| Part 3: Responder Analysis | Interactive boxplots (filterable by time point), plain-language findings, and sample-level, subject-level and per time point statistics |
| Part 4: Baseline Subset | Sample counts per project, subject counts by response and by sex, the baseline sample list, and the average B cell answer |

If the database is missing when the dashboard starts, it builds it from `cell-count.csv` first, so it also runs from a clean checkout on a hosting service.

## Repository structure

```
.
├── Makefile              setup, pipeline and dashboard targets
├── requirements.txt      Python dependencies
├── cell-count.csv        input data
├── schema.sql            SQLite schema (Part 1)
├── load_data.py          creates cell_counts.db and loads the CSV (Part 1)
├── run_pipeline.py       runs Parts 1 to 4 end to end and writes outputs/
├── app.py                Streamlit dashboard
├── teiko/
│   ├── db.py             database connection and query helpers
│   ├── analysis.py       SQL queries and statistics for Parts 2 to 4
│   └── plots.py          interactive (Plotly) and static (Matplotlib) boxplots
├── .devcontainer/        Codespaces configuration (Python 3.12, port 8501)
└── .streamlit/           Streamlit server configuration
```

Running `make pipeline` produces the following files in `outputs/`:

| File | Description |
|------|-------------|
| `part2_frequency_summary.csv` | One row per sample and population: sample, total_count, population, count, percentage |
| `part3_stats_sample_level.csv` | Responder vs non-responder tests, every sample as an observation |
| `part3_stats_subject_level.csv` | The same tests with time points averaged per subject |
| `part3_stats_by_timepoint.csv` | The same tests run separately at days 0, 7 and 14 |
| `part3_boxplots.png`, `part3_boxplots.html` | Static and interactive boxplots |
| `part4_*.csv` | Baseline sample list and the requested counts |
| `report.md` | Written summary of the Part 3 and Part 4 results |

## Part 1: Database schema

The schema (`schema.sql`) has five tables:

| Table | Grain | Columns |
|-------|-------|---------|
| `projects` | one row per project | project_id |
| `subjects` | one row per patient | subject_id, project_id, condition, age, sex, treatment, response |
| `samples` | one row per biological sample | sample_id, subject_id, sample_type, time_from_treatment_start |
| `populations` | one row per cell population | population_id, name |
| `cell_counts` | one row per sample and population | sample_id, population_id, count |

**Rationale.**

1. **No repeated metadata.** In the source file, project, condition, age, sex, treatment and response never vary within a subject, yet they are repeated on every sample row. Storing them once in `subjects` removes that redundancy and makes inconsistent values impossible. `load_data.py` verifies this property and stops with an error if a subject has conflicting metadata.
2. **Long format for counts.** Storing counts as (sample, population, count) rather than one column per population means a new population (for example regulatory T cells) is added by inserting a row in `populations`, with no schema change and no change to the analysis queries. Totals and percentages are computed with a single `GROUP BY`.
3. **Integrity.** Primary keys, foreign keys and `CHECK` constraints (non-negative counts, valid sex and response values) reject malformed data at load time.

**Scalability.** With hundreds of projects, thousands of samples and new analytics:

* Indexes on the common filter columns (condition, treatment, response, sex, sample type, time point) keep cohort selection fast as row counts grow.
* New sample-level attributes (collection site, assay batch) belong in `samples`, and new subject-level attributes (dosing arm, prior therapy) belong in `subjects`. Project-level details (sponsor, start date) can be added to `projects`. Each table grows independently.
* New measurement types (cytokine levels, gene expression) can be added as separate long-format tables keyed by `sample_id`, following the `cell_counts` pattern.
* Frequently used derived results, such as relative frequencies, can be materialised as a view or summary table. If data volume eventually exceeds what SQLite handles comfortably, the same schema migrates directly to PostgreSQL.

## Part 2: Relative frequencies

For each sample, `total_count` is the sum of the five population counts and `percentage = 100 * count / total_count`. The calculation is a single SQL query in `teiko/analysis.py` (`frequency_summary`). The full table has 52,500 rows (10,500 samples by 5 populations).

## Part 3: Responders vs non-responders

**Cohort.** Melanoma patients treated with miraclib, PBMC samples only: 1,968 samples from 656 subjects (331 responders, 325 non-responders). Each subject contributes samples at days 0, 7 and 14.

**Method.**

* Two-sided **Mann-Whitney U test** for each population. This test is non-parametric and does not assume normally distributed percentages.
* **Benjamini-Hochberg correction** across the five populations to control the false discovery rate (reported as q-values).
* **Rank-biserial correlation** as the effect size. Positive values mean responders tend to be higher.
* **Subject-level robustness check.** Samples from the same patient are not independent, which can overstate significance. Each subject's percentages are therefore also averaged across time points and the test is repeated with one observation per subject.

**Results (sample level).**

| Population | Median % (responders) | Median % (non-responders) | p-value | BH q-value | Subject-level p |
|------------|------:|------:|------:|------:|------:|
| b_cell | 9.43 | 9.79 | 0.0557 | 0.1393 | 0.3458 |
| cd8_t_cell | 24.73 | 24.60 | 0.6391 | 0.6391 | 0.6221 |
| cd4_t_cell | 30.22 | 29.66 | **0.0133** | 0.0667 | **0.0124** |
| nk_cell | 14.51 | 14.80 | 0.1211 | 0.2018 | 0.1267 |
| monocyte | 19.61 | 19.94 | 0.1632 | 0.2039 | 0.2645 |

**Conclusion.** CD4 T cells are the only population with a statistically significant difference at the 0.05 level (Mann-Whitney U, p = 0.013). Responders have a higher relative frequency of CD4 T cells than non-responders, and the result is unchanged when each subject is counted once (p = 0.012). The effect is small (median difference of about 0.6 percentage points, rank-biserial r = 0.06), and after Benjamini-Hochberg correction for five comparisons it is borderline (q = 0.067). CD4 T cell frequency is therefore the most promising candidate marker of miraclib response, but it should be confirmed in an independent cohort before it is used for prediction. No other population differs significantly.

## Part 4: Baseline subset

**Filter.** Melanoma, miraclib, PBMC, `time_from_treatment_start = 0`: **656 samples**.

| Question | Answer |
|----------|--------|
| Samples per project | prj1: 384, prj3: 272 |
| Subjects by response | responders: 331, non-responders: 325 |
| Subjects by sex | male: 344, female: 312 |

**Average B cell count** for melanoma male responders at time 0, across all sample types and treatments (485 samples): **10206.15**

## Reproducing the results without Make

```bash
pip install -r requirements.txt
python load_data.py          # creates cell_counts.db only
python run_pipeline.py       # runs Parts 1 to 4 and writes outputs/
streamlit run app.py         # starts the dashboard
```
