# Snakemake Guide for the YAC RNA-seq Pipeline

Everything in this guide uses real examples from our pipeline. If you understand this
document, you can build and modify the entire Snakefile yourself.

---

## Table of Contents

1. [What Snakemake Actually Does](#1-what-snakemake-actually-does)
2. [Rules: The Core Building Block](#2-rules-the-core-building-block)
3. [Wildcards: How One Rule Handles Many Samples](#3-wildcards-how-one-rule-handles-many-samples)
4. [Config Files: No More Hardcoded Paths](#4-config-files-no-more-hardcoded-paths)
5. [Rule Directives Reference](#5-rule-directives-reference)
6. [The DAG: How Snakemake Decides What to Run](#6-the-dag-how-snakemake-decides-what-to-run)
7. [Checkpoints: When You Don't Know the Output Until Runtime](#7-checkpoints-when-you-dont-know-the-output-until-runtime)
8. [Input Functions: Dynamic Logic in Rules](#8-input-functions-dynamic-logic-in-rules)
9. [Modular Snakefiles: include and .smk Files](#9-modular-snakefiles-include-and-smk-files)
10. [Running on SLURM Clusters](#10-running-on-slurm-clusters)
11. [envmodules: Loading Cluster Software](#11-envmodules-loading-cluster-software)
12. [Profiles: Packaging Cluster Config](#12-profiles-packaging-cluster-config)
13. [Common Commands You'll Use](#13-common-commands-youll-use)
14. [Translating a Shell Script to a Rule (Worked Example)](#14-translating-a-shell-script-to-a-rule-worked-example)
15. [Gotchas and Debugging](#15-gotchas-and-debugging)

---

## 1. What Snakemake Actually Does

Right now you have 13 scripts that you submit manually one at a time:

```bash
sbatch 1_fastqc_check.sh        # wait for it to finish...
sbatch 2_adapter_trimming.sh     # wait...
sbatch 3_post_trim_check.sh      # wait...
# ... 10 more times
```

Snakemake replaces all of this with one command:

```bash
snakemake --profile profiles/slurm
```

It figures out the order itself by looking at **which files each step needs** and
**which files each step produces**. If step 4c needs the trimmed FASTQ from step 2,
Snakemake knows to run step 2 first. If a file already exists and hasn't changed,
that step is skipped entirely.

**Key mental model:** Snakemake works **backwards** from the final files you want.
You say "I want `yac_sample_1_fwd.npy`" and it traces back through every rule needed
to produce it.

---

## 2. Rules: The Core Building Block

A **rule** is one step of your pipeline. It replaces one shell script. Here's
your `4b_genome_generation.sh`:

```bash
#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --account=rrg-cdeboer
#SBATCH --ntasks=8
#SBATCH --mem=40G

module load star/2.7.11b

GENOME_DIR=data/hybrid_star_index
YEAST_HUMAN_FILE=data/yeast_human_hybrid.fa

mkdir -p ${GENOME_DIR}

STAR --runMode genomeGenerate \
     --runThreadN 8 \
     --genomeDir ${GENOME_DIR} \
     --genomeFastaFiles ${YEAST_HUMAN_FILE} \
     --genomeSAindexNbases 14
```

Here's the same thing as a Snakemake rule:

```python
rule star_index_pass1:
    input:
        fasta = "results/genomes/yeast_human_hybrid.fa"
    output:
        genome_dir = directory("results/genomes/hybrid_star_index"),
        sa = "results/genomes/hybrid_star_index/SA"   # sentinel file
    threads: 8
    resources:
        mem_mb = 40000,
        time = "4:00:00",
        slurm_account = "rrg-cdeboer"
    envmodules:
        "star/2.7.11b"
    shell:
        """
        mkdir -p {output.genome_dir}
        STAR --runMode genomeGenerate \
             --runThreadN {threads} \
             --genomeDir {output.genome_dir} \
             --genomeFastaFiles {input.fasta} \
             --genomeSAindexNbases 14
        """
```

What changed:

| Shell script concept | Snakemake equivalent |
|---|---|
| `#SBATCH --time=4:00:00` | `resources: time = "4:00:00"` |
| `#SBATCH --mem=40G` | `resources: mem_mb = 40000` |
| `#SBATCH --ntasks=8` | `threads: 8` |
| `#SBATCH --account=...` | `resources: slurm_account = "..."` |
| `module load star/2.7.11b` | `envmodules: "star/2.7.11b"` |
| Hardcoded input path | `input: fasta = "..."` |
| Hardcoded output path | `output: genome_dir = directory("...")` |
| The actual commands | `shell: """..."""` |

Inside the `shell:` block, you reference inputs/outputs with curly braces:
`{input.fasta}`, `{output.genome_dir}`, `{threads}`.

---

## 3. Wildcards: How One Rule Handles Many Samples

This is the biggest upgrade over shell scripts. Instead of hardcoding
`yac_sample_1` everywhere, you use **wildcards** — placeholders wrapped in
`{curly_braces}`.

**Without wildcards** (what you have now — one script per sample):
```python
rule star_align_pass1:
    input:
        r1 = "results/yac_sample_1/trimmed/yac_sample_1_R1_paired.fq.gz",
        ...
    output:
        bam = "results/yac_sample_1/alignment_pass1/yac_sample_1_Aligned.sortedByCoord.out.bam"
```

**With wildcards** (one rule handles ALL samples):
```python
rule star_align_pass1:
    input:
        r1 = "results/{sample}/trimmed/{sample}_R1_paired.fq.gz",
        r2 = "results/{sample}/trimmed/{sample}_R2_paired.fq.gz",
        index_dir = "results/genomes/hybrid_star_index",
    output:
        bam = "results/{sample}/alignment_pass1/{sample}_Aligned.sortedByCoord.out.bam"
    params:
        prefix = "results/{sample}/alignment_pass1/{sample}_"
    threads: 8
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {input.index_dir} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outFileNamePrefix {params.prefix} \
             --outSAMtype BAM SortedByCoordinate
        """
```

When Snakemake needs the file
`results/yac_sample_1/alignment_pass1/yac_sample_1_Aligned.sortedByCoord.out.bam`,
it matches the output pattern and sets `{sample}` = `yac_sample_1`. That value
then fills in everywhere `{sample}` appears in input, output, params, and shell.

Add a second sample to your config and the same rule handles it automatically —
no code changes needed.

---

## 4. Config Files: No More Hardcoded Paths

Instead of editing every script when a path or parameter changes, you put
everything in `config/config.yaml`:

```yaml
samples:
  yac_sample_1:
    R1: "data/1951f6-DL-8_10_S13_L007_R1_001.fastq.gz"
    R2: "data/1951f6-DL-8_10_S13_L007_R2_001.fastq.gz"

ref:
  human_fasta: "data/hg38.fa"
  yeast_fasta: "data/sacCer3.fa"

star:
  threads: 8
  align_intron_max: 500000

modules:
  star: "star/2.7.11b"
  samtools: "samtools/1.22.1"
```

In the Snakefile, load it and access values:

```python
configfile: "config/config.yaml"

# Now config is a Python dict
SAMPLES = list(config["samples"].keys())   # ["yac_sample_1"]
RESULTS = config["results_dir"]            # "results"

rule star_align_pass1:
    input:
        r1 = lambda wc: config["samples"][wc.sample]["R1"],
        ...
    threads: config["star"]["threads"]
    envmodules:
        config["modules"]["star"]
    shell:
        """
        STAR --runThreadN {threads} \
             --alignIntronMax {config[star][align_intron_max]} \
             ...
        """
```

**To add a new sample**, just add lines to config.yaml. That's it:

```yaml
samples:
  yac_sample_1:
    R1: "data/1951f6-DL-8_10_S13_L007_R1_001.fastq.gz"
    R2: "data/1951f6-DL-8_10_S13_L007_R2_001.fastq.gz"
  yac_sample_2:                           # <-- new
    R1: "data/sample2_R1_001.fastq.gz"    # <-- new
    R2: "data/sample2_R2_001.fastq.gz"    # <-- new
```

---

## 5. Rule Directives Reference

Every directive you'll need for this pipeline:

```python
rule example_rule:
    input:                              # Files this rule NEEDS (dependencies)
        bam = "results/{sample}/file.bam"
    output:                             # Files this rule CREATES
        bw = "results/{sample}/file.bw"
    params:                             # Values passed to shell (not tracked as files)
        region = "chr1:100-200"
    threads: 8                          # CPUs (maps to SLURM --cpus-per-task)
    resources:                          # Cluster resources
        mem_mb = 32000,                 #   maps to SLURM --mem
        time = "4:00:00",              #   maps to SLURM --time
        slurm_account = "rrg-cdeboer"  #   maps to SLURM --account
    envmodules:                         # module load commands
        "samtools/1.22.1",
        "python/3.11.5"
    log:                                # Log file path (Snakemake tracks it)
        "logs/{sample}/example_rule.log"
    shell:                              # The actual commands
        """
        bamCoverage -b {input.bam} \
            -o {output.bw} \
            --numberOfProcessors {threads} \
            2> {log}
        """
```

### Key things to know:
- **`input`** and **`output`** are how Snakemake builds the dependency graph. If
  rule A's output matches rule B's input, A runs before B.
- **`output: directory("path")`** — use this when a tool creates a directory
  (like STAR genome indexes). Add a sentinel file too so Snakemake can verify it.
- **`params`** — for values that aren't files. Snakemake doesn't track these for
  dependency resolution, but you can use wildcards in them.
- **`threads`** — automatically passed to SLURM. Reference as `{threads}` in shell.
- **`resources`** — the SLURM executor maps these to sbatch flags.
- **`log`** — Snakemake keeps log files even when you re-run.

---

## 6. The DAG: How Snakemake Decides What to Run

DAG = Directed Acyclic Graph. It's just a fancy way of saying "dependency tree."

You define a **target rule** (usually called `rule all`) that lists every final
output you want:

```python
rule all:
    input:
        expand("results/{sample}/final/{sample}_fwd.npy", sample=SAMPLES),
        expand("results/{sample}/final/{sample}_rev.npy", sample=SAMPLES),
        expand("results/{sample}/final/human_yac_insert.fa", sample=SAMPLES),
```

`expand()` is a helper that generates a list from a pattern:
```python
expand("results/{sample}/final/{sample}_fwd.npy", sample=["yac_sample_1", "yac_sample_2"])
# produces:
# ["results/yac_sample_1/final/yac_sample_1_fwd.npy",
#  "results/yac_sample_2/final/yac_sample_2_fwd.npy"]
```

Snakemake sees `rule all` wants `yac_sample_1_fwd.npy`. It asks: "which rule
produces that?" It finds `bigwig_to_numpy`. That rule needs BigWig files, so it
asks: "which rule produces those?" It finds `bamcoverage`. And so on, all the
way back to the raw FASTQ files.

**Visualize it:**
```bash
snakemake --dag | dot -Tpng > dag.png
```

**Dry run** (show what would execute without running anything):
```bash
snakemake -n
```

**Key behavior:** If an output file already exists and is newer than all its
inputs, Snakemake **skips** that rule. This means if your alignment finished but
bamCoverage failed, re-running only re-runs bamCoverage and everything after it.

---

## 7. Checkpoints: When You Don't Know the Output Until Runtime

This is the hardest concept and the one most relevant to your pipeline.

**The problem:** Step 4d (`find_yac_region`) discovers which chromosome the YAC
is on. Before it runs, you literally don't know if it's chr1, chr5, or chr17.
But step 4e needs that chromosome name to build the refined genome.

In a normal rule, Snakemake resolves the entire DAG **before** running anything.
It can't do that here because downstream parameters depend on a file that doesn't
exist yet.

**The solution: `checkpoint`**

A checkpoint is a rule that tells Snakemake: "pause DAG resolution here, run me,
then re-evaluate."

```python
# Declare it as checkpoint instead of rule
checkpoint find_yac_region:
    input:
        bam = "results/{sample}/alignment_pass1/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        bed = "results/{sample}/yac_discovery/my_yac.bed",
        coverage = "results/{sample}/yac_discovery/coverage_data.txt"
    envmodules:
        config["modules"]["samtools"]
    shell:
        """
        samtools index {input.bam}

        TOP_CHR=$(samtools idxstats {input.bam} | \
            grep "^chr" | awk '$2 > 10000000' | \
            sort -k3,3rn | head -n 1 | cut -f1)

        NUM_READS=$(samtools view -c {input.bam} "$TOP_CHR")
        START_INDEX=$(( NUM_READS / 100 ))
        END_INDEX=$(( NUM_READS * 99 / 100 ))

        read START END <<< $(samtools view {input.bam} "$TOP_CHR" 2>/dev/null | \
            awk '{{print $4}}' | sort -n | \
            sed -n "${{START_INDEX}}p;${{END_INDEX}}p" | tr "\\n" " ")

        echo -e "${{TOP_CHR}}\\t${{START}}\\t${{END}}" > {output.bed}

        samtools depth -a -r "${{TOP_CHR}}:${{START}}-${{END}}" {input.bam} | \
            awk 'NR % 100 == 0' > {output.coverage}
        """
```

**IMPORTANT shell escaping:** Inside a Snakemake `shell:` block, curly braces
`{}` are interpreted as Snakemake placeholders. To use literal braces (like in
`awk '{print $4}'`), you **double them**: `awk '{{print $4}}'`.

Then downstream rules access the checkpoint result through a helper function:

```python
def get_yac_bed_values(wildcards):
    """Read my_yac.bed AFTER the checkpoint has run."""
    # This .get() call tells Snakemake to wait for the checkpoint
    bed_file = checkpoints.find_yac_region.get(sample=wildcards.sample).output.bed
    with open(bed_file) as f:
        parts = f.read().strip().split("\t")
        return {"chrom": parts[0], "start": parts[1], "end": parts[2]}


rule build_refined_genome:
    input:
        # This triggers the checkpoint — Snakemake won't run this rule
        # until find_yac_region has finished
        bed = lambda wc: checkpoints.find_yac_region.get(sample=wc.sample).output.bed,
        human = config["ref"]["human_fasta"],
        yeast = config["ref"]["yeast_fasta"]
    output:
        chimeric = "results/{sample}/genomes/yeast_human_chimeric.fa"
    params:
        # get_yac_bed_values reads the BED file and returns the coordinates
        region = lambda wc: get_yac_bed_values(wc)
    shell:
        """
        samtools faidx {input.human} \
            "{params.region[chrom]}:{params.region[start]}-{params.region[end]}" | \
            sed 's/^>.*/>human_yac/' > {output.chimeric}.tmp
        cat {input.yeast} {output.chimeric}.tmp > {output.chimeric}
        rm -f {output.chimeric}.tmp
        """
```

**The flow:**
1. Snakemake builds the DAG up to `find_yac_region` — everything before it is
   fully resolved
2. `find_yac_region` runs, produces `my_yac.bed` with (e.g.) `chr1  100  50000`
3. Snakemake re-evaluates the DAG — now `get_yac_bed_values()` can read the file
4. Downstream rules (`build_refined_genome`, `extract_yac_fasta`, `bigwig_to_numpy`)
   get their actual parameter values

---

## 8. Input Functions: Dynamic Logic in Rules

Sometimes you need Python logic to determine inputs. Use a **lambda** or a
**named function** in the `input:` block.

**Lambda example** — looking up FASTQ paths from config:
```python
rule trimmomatic:
    input:
        r1 = lambda wc: config["samples"][wc.sample]["R1"],
        r2 = lambda wc: config["samples"][wc.sample]["R2"]
    output:
        r1_paired = "results/{sample}/trimmed/{sample}_R1_paired.fq.gz",
        ...
```

`wc` (or `wildcards`) is an object with attributes for each wildcard. If
`{sample}` is `yac_sample_1`, then `wc.sample` is the string `"yac_sample_1"`.

**Named function example** — checkpoint access:
```python
def get_yac_bed_path(wildcards):
    return checkpoints.find_yac_region.get(sample=wildcards.sample).output.bed

rule extract_yac_fasta:
    input:
        bed = get_yac_bed_path,       # function reference, no parentheses
        human = config["ref"]["human_fasta"]
    ...
```

**When to use input functions:**
- Looking up sample-specific values from config (`lambda wc: config["samples"][wc.sample]["R1"]`)
- Accessing checkpoint outputs (`checkpoints.X.get(...)`)
- Any logic that depends on the wildcard value

---

## 9. Modular Snakefiles: include and .smk Files

With 14 rules, a single Snakefile gets unwieldy. Split it up:

```
Snakefile                       # Main entry point
workflow/rules/qc.smk           # QC rules
workflow/rules/genome_prep.smk  # Genome concatenation
workflow/rules/alignment_pass1.smk
...
```

In the main `Snakefile`:
```python
configfile: "config/config.yaml"

SAMPLES = list(config["samples"].keys())
RESULTS = config["results_dir"]

# Include rule files — they share the same namespace (SAMPLES, RESULTS, config)
include: "workflow/rules/qc.smk"
include: "workflow/rules/genome_prep.smk"
include: "workflow/rules/alignment_pass1.smk"
include: "workflow/rules/yac_discovery.smk"
include: "workflow/rules/alignment_pass2.smk"
include: "workflow/rules/coverage.smk"
include: "workflow/rules/final_outputs.smk"

rule all:
    input:
        expand(f"{RESULTS}/{{sample}}/final/{{sample}}_fwd.npy", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/final/{{sample}}_rev.npy", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/final/human_yac_insert.fa", sample=SAMPLES),
```

Each `.smk` file is just rules — no `configfile:` or other setup. They inherit
all variables from the main `Snakefile`.

---

## 10. Running on SLURM Clusters

Snakemake 8+ uses **executor plugins**. For SLURM, install the plugin:

```bash
pip install snakemake snakemake-executor-plugin-slurm
```

Then tell Snakemake to use it:
```bash
snakemake --executor slurm --configfile config/config.yaml
```

Snakemake translates your rule's `resources:` and `threads:` into `sbatch` flags
automatically:

| Rule directive | sbatch flag |
|---|---|
| `threads: 8` | `--cpus-per-task=8` |
| `resources: mem_mb=40000` | `--mem=40000M` |
| `resources: time="4:00:00"` | `--time=4:00:00` |
| `resources: slurm_account="rrg-cdeboer"` | `--account=rrg-cdeboer` |

**Snakemake itself runs on the login node** (or in a tmux/screen session). It's
lightweight — it just submits jobs and monitors them. Each rule execution is a
separate SLURM job.

**Alliance-specific gotcha — NFS latency:** When a job finishes, the output file
may not be visible on the login node for a few seconds (shared filesystem delay).
Use `--latency-wait 120` so Snakemake waits before checking for output files.

---

## 11. envmodules: Loading Cluster Software

On Alliance clusters, software is loaded with `module load`. Snakemake supports
this natively:

```python
rule star_align_pass1:
    envmodules:
        "star/2.7.11b",
        "samtools/1.22.1"
    shell:
        """
        # star and samtools are now available here
        STAR --runThreadN {threads} ...
        samtools index {output.bam}
        """
```

This is equivalent to adding `module load star/2.7.11b` at the top of your shell
script. You must enable it with `--use-envmodules` (or set it in your profile).

**For Python packages** (deeptools, pyBigWig) that aren't available as modules,
activate your pip venv in the shell block:

```python
rule bamcoverage:
    envmodules:
        config["modules"]["python"],
        config["modules"]["samtools"]
    shell:
        """
        source {config[python_venv]}/bin/activate
        bamCoverage -b {input.bam} -o {output.bw} ...
        """
```

---

## 12. Profiles: Packaging Cluster Config

Instead of typing `--executor slurm --latency-wait 120 --use-envmodules --jobs 20`
every time, save it in a profile:

**`profiles/slurm/config.yaml`:**
```yaml
executor: slurm

default-resources:
  slurm_account: "def-cdeboer"
  mem_mb: 16000
  time: "2:00:00"

jobs: 20                  # Max concurrent SLURM jobs
latency-wait: 120         # Wait for NFS propagation (seconds)
retries: 1                # Retry failed jobs once
use-envmodules: true       # Enable the envmodules: directive
printshellcmds: true       # Print commands for debugging
```

Now you just run:
```bash
snakemake --configfile config/config.yaml --profile profiles/slurm
```

The `default-resources` apply to every rule unless the rule overrides them in its
own `resources:` block. So a simple FastQC rule gets 16GB/2h by default, while
STAR alignment overrides to 40GB/8h.

---

## 13. Common Commands You'll Use

```bash
# Dry run — show what would execute (ALWAYS do this first)
snakemake -n --profile profiles/slurm

# Dry run with reasons — explains WHY each rule would run
snakemake -n -r --profile profiles/slurm

# Full run
snakemake --profile profiles/slurm

# Run just one specific target
snakemake results/yac_sample_1/final/yac_sample_1_fwd.npy --profile profiles/slurm

# Visualize the DAG
snakemake --dag | dot -Tpng > dag.png

# Visualize the rulegraph (simpler, no sample expansion)
snakemake --rulegraph | dot -Tpng > rulegraph.png

# Force re-run a specific rule (even if output exists)
snakemake --forcerun star_align_pass1 --profile profiles/slurm

# Force re-run a rule AND everything downstream of it
snakemake --forcerun star_align_pass1 -R --profile profiles/slurm

# See the status of all running/pending jobs
snakemake --profile profiles/slurm --summary

# Unlock after a crash (Snakemake locks the working directory)
snakemake --unlock

# Clean all outputs (careful!)
snakemake --delete-all-output
```

---

## 14. Translating a Shell Script to a Rule (Worked Example)

Let's walk through converting `5_create_bw.sh` step by step.

### Original script:
```bash
#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --account=rrg-cdeboer
#SBATCH --ntasks=8
#SBATCH --mem=32G

SAMPLE="yac_sample_1_"
BAM="outputs/alignment_2/${SAMPLE}Aligned.sortedByCoord.out.bam"
OUTDIR="outputs/bigwig/${SAMPLE}"
mkdir -p ${OUTDIR}

module load samtools
samtools index ${BAM}

bamCoverage -b ${BAM} \
    -o ${OUTDIR}/${SAMPLE}_forward.bw \
    -of bigwig \
    --filterRNAstrand forward \
    --binSize 1 \
    --normalizeUsing CPM \
    --numberOfProcessors 8

bamCoverage -b ${BAM} \
    -o ${OUTDIR}/${SAMPLE}_reverse.bw \
    -of bigwig \
    --filterRNAstrand reverse \
    --binSize 1 \
    --normalizeUsing CPM \
    --numberOfProcessors 8
```

### Step-by-step translation:

**a) Identify inputs and outputs:**
- Input: BAM file from alignment pass 2 + its index (.bai)
- Outputs: forward BigWig, reverse BigWig

**b) Replace hardcoded sample name with `{sample}` wildcard.**

**c) Replace SBATCH directives with `resources:`/`threads:`.**

**d) Replace `module load` with `envmodules:`.**

**e) Use a `{strand}` wildcard so one rule handles both forward and reverse:**

```python
rule samtools_index:
    input:
        bam = "results/{sample}/alignment_pass2/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        bai = "results/{sample}/alignment_pass2/{sample}_Aligned.sortedByCoord.out.bam.bai"
    envmodules:
        config["modules"]["samtools"]
    shell:
        "samtools index {input.bam}"


rule bamcoverage:
    input:
        bam = "results/{sample}/alignment_pass2/{sample}_Aligned.sortedByCoord.out.bam",
        bai = "results/{sample}/alignment_pass2/{sample}_Aligned.sortedByCoord.out.bam.bai"
    output:
        bw = "results/{sample}/bigwig/{sample}_{strand}.bw"
    wildcard_constraints:
        strand = "forward|reverse"     # Only match these two values
    threads: config["bamcoverage"]["threads"]
    resources:
        mem_mb = 32000,
        time = "4:00:00",
        slurm_account = config["cluster"]["account"]
    envmodules:
        config["modules"]["python"],
        config["modules"]["samtools"]
    shell:
        """
        source {config[python_venv]}/bin/activate
        bamCoverage -b {input.bam} \
            -o {output.bw} \
            -of bigwig \
            --filterRNAstrand {wildcards.strand} \
            --binSize {config[bamcoverage][bin_size]} \
            --normalizeUsing {config[bamcoverage][normalize_using]} \
            --numberOfProcessors {threads}
        """
```

The `{strand}` wildcard trick: when `bigwig_to_numpy` requests both
`{sample}_forward.bw` and `{sample}_reverse.bw` as inputs, Snakemake
instantiates the `bamcoverage` rule twice — once with `strand=forward`,
once with `strand=reverse`. They can run as two parallel SLURM jobs.

`wildcard_constraints: strand = "forward|reverse"` is a regex that prevents
Snakemake from accidentally matching other strings.

---

## 15. Gotchas and Debugging

### Curly brace escaping
Snakemake interprets `{` and `}` in `shell:` blocks. To use them literally
(awk, bash variables), double them:

```python
# WRONG — Snakemake tries to resolve {print $4} as a placeholder
shell: "awk '{print $4}' file.txt"

# CORRECT — escaped braces
shell: "awk '{{print $4}}' file.txt"

# CORRECT — bash variable in a shell block
shell: "echo ${{MY_VAR}}"
```

### Snakemake locks the directory
If Snakemake crashes or you Ctrl-C, it leaves a `.snakemake/locks` directory.
Next run will fail with "Directory is locked." Fix: `snakemake --unlock`

### "Missing input files" error
Usually means an upstream rule hasn't been defined yet, or there's a typo in the
filename pattern. Check that the output of rule A exactly matches the input of
rule B — character for character.

### NFS latency on Alliance
Output file is "missing" right after a SLURM job finishes. The fix is
`latency-wait: 120` in your profile. If you still see it, bump to 300.

### STAR _STARtmp conflicts
If STAR is interrupted, it leaves behind `_STARtmp/` directories. Snakemake
retries will fail because STAR refuses to overwrite them. Fix: use
`--outTmpDir $SLURM_TMPDIR/STARtmp_{wildcards.sample}` in the STAR command.
`$SLURM_TMPDIR` is local SSD on the compute node, fast and auto-cleaned.

### rule all doesn't run anything
`rule all` must be the **first** rule in the Snakefile (or the first included).
Snakemake uses the first rule as the default target.

### Debugging a specific rule
```bash
# Print the exact shell command that would be executed
snakemake -n -p results/yac_sample_1/bigwig/yac_sample_1_forward.bw

# Run just that one rule locally (no SLURM)
snakemake --cores 8 results/yac_sample_1/bigwig/yac_sample_1_forward.bw
```

### Checking SLURM job status
Snakemake monitors jobs automatically, but you can also check:
```bash
squeue -u $USER     # your running/pending jobs
sacct -j <jobid>    # details on a specific job
```

---

## Quick Reference: Our Pipeline's Rules

| Rule | Input | Output | Cluster resources |
|---|---|---|---|
| `fastqc_pretrim` | Raw FASTQ | HTML report | 4 CPU, 16GB, 2h |
| `trimmomatic` | Raw paired FASTQ | Trimmed paired FASTQ | 4 CPU, 16GB, 2h |
| `fastqc_posttrim` | Trimmed FASTQ | HTML report | 4 CPU, 16GB, 2h |
| `concatenate_genomes` | hg38 + sacCer3 | Hybrid FASTA | 1 CPU, 16GB, 2h |
| `star_index_pass1` | Hybrid FASTA | STAR index dir | 8 CPU, 40GB, 4h |
| `star_align_pass1` | Trimmed FASTQ + index | BAM | 8 CPU, 40GB, 8h |
| `find_yac_region` | BAM | my_yac.bed | 4 CPU, 40GB, 2h |
| `plot_coverage` | coverage_data.txt | PNG | 1 CPU, 4GB, 0.5h |
| `build_refined_genome` | BED + hg38 + sacCer3 | Chimeric FASTA | 1 CPU, 16GB, 1h |
| `star_index_pass2` | Chimeric FASTA | STAR index dir | 8 CPU, 32GB, 2h |
| `star_align_pass2` | Trimmed FASTQ + index | BAM | 8 CPU, 40GB, 4h |
| `samtools_index` | BAM | BAI | 1 CPU, 8GB, 1h |
| `bamcoverage` | BAM + BAI | BigWig (per strand) | 8 CPU, 32GB, 4h |
| `extract_yac_fasta` | BED + hg38 | FASTA | 1 CPU, 8GB, 1h |
| `bigwig_to_numpy` | BigWig (fwd+rev) + BED | NumPy arrays | 1 CPU, 16GB, 2h |