# Snakemake Pipeline for YAC RNA-seq Analysis

## Context

The current pipeline is 13 sequential shell scripts (`1_fastqc_check.sh` through `7_generate_numpy.sh`) that each submit individual SLURM jobs manually. This makes it error-prone (no automatic dependency tracking, manual re-runs on failure, hardcoded paths) and hard to scale to multiple samples. Converting to Snakemake gives us automatic dependency resolution, SLURM integration, multi-sample support, and reproducibility — all critical for running on Alliance/Compute Canada clusters.

## Pipeline DAG (Current Steps → Snakemake Rules)

```
fastqc_pretrim (1)        concatenate_genomes (4a)  [shared, runs once]
       |                          |
  trimmomatic (2)         star_index_pass1 (4b)
       |                          |
fastqc_posttrim (3)       star_align_pass1 (4c)
                                  |
                        find_yac_region (4d)  ← CHECKPOINT
                         /        |         \
              plot_coverage  build_refined_genome (4e)  extract_yac_fasta (6)
                                  |
                          star_index_pass2 (4f)
                                  |
                          star_align_pass2 (4g)
                                  |
                          samtools_index
                             /        \
              bamcoverage_fwd (5)   bamcoverage_rev (5)
                             \        /
                        bigwig_to_numpy (7)
```

## File Structure

```
yac-rnaseq-pipeline/
├── Snakefile                          # Entry point, includes rules, defines targets
├── config/
│   └── config.yaml                    # All configurable parameters
├── workflow/
│   ├── rules/
│   │   ├── qc.smk                     # fastqc_pretrim, trimmomatic, fastqc_posttrim
│   │   ├── genome_prep.smk            # concatenate_genomes
│   │   ├── alignment_pass1.smk        # star_index_pass1, star_align_pass1
│   │   ├── yac_discovery.smk          # find_yac_region (checkpoint), plot_coverage
│   │   ├── alignment_pass2.smk        # build_refined_genome, star_index_pass2, star_align_pass2
│   │   ├── coverage.smk               # samtools_index, bamcoverage (fwd/rev via wildcard)
│   │   └── final_outputs.smk          # extract_yac_fasta, bigwig_to_numpy
│   └── scripts/
│       ├── bw_to_np.py                # Existing, modified for --sample arg
│       └── graph_yac.py               # Existing, modified for headless execution
├── profiles/
│   └── slurm/
│       └── config.yaml                # Snakemake SLURM executor profile
└── README.md                          # How to configure and run
```

## config/config.yaml

```yaml
# --- Samples (add more entries to process multiple samples) ---
samples:
  yac_sample_1:
    R1: "data/1951f6-DL-8_10_S13_L007_R1_001.fastq.gz"
    R2: "data/1951f6-DL-8_10_S13_L007_R2_001.fastq.gz"

# --- Reference Genomes ---
ref:
  human_fasta: "data/hg38.fa"
  yeast_fasta: "data/sacCer3.fa"

# --- Output ---
results_dir: "results"

# --- Cluster ---
cluster:
  account: "def-cdeboer"
  account_large: "rrg-cdeboer"  # For STAR indexing/alignment (high memory)

# --- Python venv path (for deeptools, pyBigWig, etc.) ---
python_venv: "/path/to/venvs/yac-rnaseq"

# --- Tool Parameters ---
trimmomatic:
  adapters: "TruSeq3-PE.fa"
  clip: "2:30:10"
  leading: 3
  trailing: 3
  slidingwindow: "4:15"
  minlen: 36
  threads: 4

fastqc:
  threads: 4

star:
  threads: 8
  genome_sa_index_nbases: 14       # For full hybrid genome (~3.3GB)
  genome_sa_index_nbases_refined: 11  # For refined chimeric genome (~13MB)
  align_intron_max: 500000

bamcoverage:
  bin_size: 1
  normalize_using: "CPM"
  threads: 8

bigwig_to_numpy:
  bin_size: 1

# --- Module Versions (envmodules) ---
modules:
  fastqc: "fastqc/0.12.1"
  trimmomatic: "trimmomatic/0.39"
  star: "star/2.7.11b"
  samtools: "samtools/1.22.1"
  python: "python/3.11.5"
```

## Key Design Decisions

### 1. Dynamic YAC Discovery → Snakemake `checkpoint`

Step 4d discovers which human chromosome the YAC maps to. Downstream rules need this info to build the refined genome. This is handled with a Snakemake **checkpoint**:

- `find_yac_region` is declared as `checkpoint` (not `rule`)
- It outputs `my_yac.bed` containing `<chrom>\t<start>\t<end>`
- Downstream rules use an input function that calls `checkpoints.find_yac_region.get(sample=wc.sample).output.bed` to trigger DAG re-evaluation
- A helper `get_yac_bed_values(wildcards)` reads the BED file and returns `{"chrom": ..., "start": ..., "end": ...}` for use in `params:` blocks

### 2. Environment Management → `envmodules` + pip venv

- Cluster tools (STAR, samtools, FastQC, Trimmomatic) → `envmodules:` directive in each rule
- Python packages (deeptools, pyBigWig, numpy, matplotlib) → pip virtualenv, activated in shell blocks via `source {config[python_venv]}/bin/activate`
- This follows Alliance best practices (avoids conda conflicts with Lmod)

### 3. SLURM Integration → Snakemake 8+ Executor Plugin Profile

`profiles/slurm/config.yaml`:
```yaml
executor: slurm
default-resources:
  slurm_account: "def-cdeboer"
  mem_mb: 16000
  time: "2:00:00"
jobs: 20
latency-wait: 120    # Critical for NFS propagation on Alliance clusters
retries: 1
use-envmodules: true
```

Per-rule resource overrides (mem, time, account) are specified in each rule's `resources:` block. High-memory jobs (STAR index/align) use `rrg-cdeboer`.

### 4. Shared vs Per-Sample Rules

- **Shared** (run once): `concatenate_genomes`, `star_index_pass1` — the hybrid genome is sample-independent
- **Per-sample** (everything else): The refined genome is sample-specific since each sample's YAC may map to a different chromosome

## Implementation Steps

### Step 1: Create directory structure
Create `config/`, `workflow/rules/`, `workflow/scripts/`, `profiles/slurm/`.

### Step 2: Write `config/config.yaml`
As shown above — all hardcoded paths and parameters become configurable.

### Step 3: Write `profiles/slurm/config.yaml`
SLURM executor profile with Alliance-specific settings.

### Step 4: Write `Snakefile`
Top-level file that:
- Loads config
- Defines `SAMPLES` and `RESULTS` variables
- Defines `get_yac_bed_values()` helper
- Includes all `.smk` rule files
- Defines `rule all` with target outputs

### Step 5: Write `workflow/rules/qc.smk`
Rules: `fastqc_pretrim`, `trimmomatic`, `fastqc_posttrim`
- Fix bug: post-trim QC now correctly outputs to a `posttrim/` directory

### Step 6: Write `workflow/rules/genome_prep.smk`
Rule: `concatenate_genomes`
- Prefixes yeast chromosome headers with `yeast_`
- Concatenates with human genome, indexes result

### Step 7: Write `workflow/rules/alignment_pass1.smk`
Rules: `star_index_pass1`, `star_align_pass1`
- Standardized output naming with trailing underscore in STAR prefix

### Step 8: Write `workflow/rules/yac_discovery.smk`
Checkpoint: `find_yac_region` + Rule: `plot_coverage`
- Translates the complex shell logic from `4d_find_human_yac.sh`
- Outputs `my_yac.bed` for downstream consumption

### Step 9: Write `workflow/rules/alignment_pass2.smk`
Rules: `build_refined_genome`, `star_index_pass2`, `star_align_pass2`
- Uses `get_yac_bed_values()` to read checkpoint output
- Uses `genome_sa_index_nbases_refined: 11` (smaller genome needs smaller index)
- Uses STAR `--outTmpDir $SLURM_TMPDIR/STARtmp` to avoid stale temp dir issues on retries

### Step 10: Write `workflow/rules/coverage.smk`
Rules: `samtools_index`, `bamcoverage`
- Fix bug: uses second-pass BAM (not first-pass as in current `5_create_bw.sh`)
- Single `bamcoverage` rule with `{strand}` wildcard for forward/reverse

### Step 11: Write `workflow/rules/final_outputs.smk`
Rules: `extract_yac_fasta`, `bigwig_to_numpy`
- Both use checkpoint-derived coordinates

### Step 12: Modify existing Python scripts
- `bw_to_np.py`: Add `--sample` argument for flexible output naming
- `graph_yac.py`: Add `matplotlib.use('Agg')` for headless execution, add `--input`/`--output` CLI args, remove `plt.show()`

### Step 13: Write `README.md`
Setup instructions: venv creation, module availability check, dry-run, full-run commands.

## Bugs Fixed vs Current Pipeline
1. **`5_create_bw.sh`**: References first-pass alignment BAM → fixed to use second-pass
2. **`3_post_trim_check.sh`**: Outputs to `before_trim/` directory → fixed to `posttrim/`
3. **STAR prefix inconsistency**: `4c` uses no trailing underscore, `4g` uses one → standardized
4. **`graph_yac.py`**: Calls `plt.show()` which fails on headless clusters → removed

## How to Run

```bash
# One-time setup: install Snakemake 8+ and the SLURM plugin
module load python/3.11.5
python -m venv --system-site-packages ~/snakemake-env
source ~/snakemake-env/bin/activate
pip install snakemake snakemake-executor-plugin-slurm

# One-time setup: create pipeline Python venv
python -m venv --system-site-packages ~/project-venvs/yac-rnaseq
source ~/project-venvs/yac-rnaseq/bin/activate
pip install deeptools pyBigWig numpy matplotlib pandas

# Dry run (verify DAG)
snakemake --configfile config/config.yaml --profile profiles/slurm -n

# Visualize DAG
snakemake --configfile config/config.yaml --dag | dot -Tpng > dag.png

# Full run (submit all jobs to SLURM)
snakemake --configfile config/config.yaml --profile profiles/slurm

# Run specific sample only
snakemake --configfile config/config.yaml --profile profiles/slurm \
    results/yac_sample_1/final/yac_sample_1_fwd.npy
```

## Verification

1. **Dry run**: `snakemake -n` — verify DAG looks correct, all rules resolve
2. **DAG visualization**: Check the PNG for expected rule dependencies
3. **Single-sample run**: Execute full pipeline, verify outputs match current results
4. **Check final outputs exist**: `results/yac_sample_1/final/{human_yac_insert.fa, yac_sample_1_fwd.npy, yac_sample_1_rev.npy}`
5. **Multi-sample test**: Add a second sample to config.yaml, verify independent DAGs
