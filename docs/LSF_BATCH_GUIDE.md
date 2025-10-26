# LSF Batch Factorizer for noLZSS

## Overview

The LSF batch factorizer (`lsf_batch_factorize.py`) is a specialized script for submitting large-scale FASTA factorization jobs to LSF (Load Sharing Facility) clusters. It automatically estimates computational resources based on benchmarking data and efficiently manages job submission and monitoring.

## Features

- **Intelligent Resource Estimation**: Uses benchmark trend parameters to predict time, memory, and disk requirements
- **Automatic Resource Optimization**: Allocates resources efficiently to avoid waste while ensuring job completion
- **Job Monitoring**: Tracks job status and provides completion summaries
- **Minimal Logging Overhead**: Consolidates logs in a single directory without spamming user email
- **Failure Tracking**: Generates concise summaries of failed jobs for easy debugging
- **Dry-Run Mode**: Preview job specifications before submission

## Prerequisites

1. **LSF Cluster Access**: You must have access to an LSF cluster with `bsub` and `bjobs` commands
2. **noLZSS Installation**: The noLZSS package must be installed and accessible on compute nodes
3. **Benchmark Data** (optional): Run benchmarks to generate trend parameters for accurate resource estimation

## Quick Start

### 1. Create a File List

Create a text file containing paths to your FASTA files (one per line):

```bash
# files.txt
/path/to/genome1.fasta
/path/to/genome2.fasta
/path/to/genome3.fasta
```

### 2. Submit Jobs

```bash
python -m noLZSS.genomics.lsf_batch_factorize \
    --file-list files.txt \
    --output-dir results \
    --mode with_reverse_complement
```

### 3. Monitor Progress

The script will automatically monitor jobs and report when they complete. Use `Ctrl+C` if you want to exit early (jobs will continue running).

## Usage

### Basic Usage

```bash
python -m noLZSS.genomics.lsf_batch_factorize \
    --file-list <file_list.txt> \
    --output-dir <output_directory>
```

### Advanced Usage

```bash
python -m noLZSS.genomics.lsf_batch_factorize \
    --file-list files.txt \
    --output-dir results \
    --mode with_reverse_complement \
    --queue normal \
    --max-threads 4 \
    --safety-factor 2.0 \
    --log-dir lsf_logs \
    --trend-file custom_trends.pkl
```

## Command-Line Arguments

### Required Arguments

- `--file-list PATH`: Text file containing list of FASTA files to process (one per line)
- `--output-dir PATH`: Output directory for binary factorization results

### Factorization Options

- `--mode {with_reverse_complement,without_reverse_complement}`: Factorization mode (default: `with_reverse_complement`)
  - `with_reverse_complement`: Consider reverse complement in factorization (DNA/RNA)
  - `without_reverse_complement`: Standard factorization

### LSF Configuration

- `--queue QUEUE`: LSF queue name (uses cluster default if not specified)
- `--max-threads N`: Maximum number of threads per job (default: 1)
- `--job-name-prefix PREFIX`: Prefix for LSF job names (default: `nolzss`)
- `--log-dir PATH`: Directory for LSF job logs (default: `lsf_logs`)

### Resource Estimation

- `--trend-file PATH`: Path to trend parameters file (`.pkl` or `.json`)
  - If not specified, looks in default locations:
    - `benchmarks/fasta_results/trend_parameters.pkl`
    - `fasta_results/trend_parameters.pkl`
- `--safety-factor FACTOR`: Safety multiplier for resource estimates (default: 1.5)
  - Higher values provide more buffer but may wait longer for resources
  - Recommended: 1.3-2.0

### Execution Control

- `--no-wait`: Submit jobs and exit without waiting for completion
- `--check-interval SECONDS`: Interval between job status checks (default: 60)
- `--dry-run`: Print job specifications without submitting

### Logging

- `--log-level {DEBUG,INFO,WARNING,ERROR}`: Logging verbosity (default: INFO)
- `--log-file PATH`: Write logs to file (in addition to console)

## Output Structure

The script creates the following directory structure:

```
output_dir/
├── with_reverse_complement/    # Output for with_rc mode
│   ├── genome1.bin
│   ├── genome2.bin
│   └── genome3.bin
└── without_reverse_complement/ # Output for without_rc mode
    ├── genome1.bin
    ├── genome2.bin
    └── genome3.bin

lsf_logs/                      # LSF job logs
├── nolzss_1_genome1.log       # Individual job logs
├── nolzss_2_genome2.log
├── nolzss_3_genome3.log
├── job_tracking.json          # Job tracking information
└── job_summary.txt            # Final summary report
```

## Resource Estimation

### How It Works

The script uses benchmark trend parameters to predict:
- **Execution Time**: Based on input size (nucleotides)
- **Memory Usage**: Peak memory requirement
- **Disk Space**: Output file size

### Formula

For a FASTA file with `N` nucleotides:

```
Time (seconds) = trend_slope * N + trend_intercept
Memory (MB) = trend_slope * N + trend_intercept
Disk (MB) = trend_slope * N + trend_intercept

LSF Request = Estimate * safety_factor (rounded up)
```

### Example Estimates

| Input Size | Estimated Time | Estimated Memory | Cluster Memory | Cluster Time |
|------------|---------------|------------------|----------------|--------------|
| 10 kbp     | 0.3 sec       | 0.2 MB          | 1 GB           | 5 min        |
| 100 kbp    | 2.7 sec       | 1.4 MB          | 1 GB           | 5 min        |
| 1 Mbp      | 27 sec        | 14 MB           | 1 GB           | 5 min        |
| 10 Mbp     | 270 sec       | 140 MB          | 1 GB           | 10 min       |
| 100 Mbp    | 2700 sec      | 1.4 GB          | 2 GB           | 75 min       |

*Note: Actual values depend on your benchmark results*

## Generating Benchmark Data

To get accurate resource estimates for your cluster:

### 1. Run Benchmarks

```bash
cd benchmarks
python fasta_benchmark.py --output-dir fasta_results
```

This generates `fasta_results/trend_parameters.pkl` with scaling coefficients.

### 2. Verify Predictions

```bash
python fasta_predictor.py fasta_results/trend_parameters.pkl \
    --size 1000000 \
    --function write_factors_binary_file_fasta_multiple_dna_w_rc
```

### 3. Use Custom Trends

```bash
python -m noLZSS.genomics.lsf_batch_factorize \
    --file-list files.txt \
    --output-dir results \
    --trend-file benchmarks/fasta_results/trend_parameters.pkl
```

## Job Monitoring

### During Execution

The script displays:
- Job submission progress
- Job IDs and names
- Estimated resources
- Running job count

### After Completion

The script generates:
- **Console Summary**: Success/failure counts
- **job_summary.txt**: Detailed report with failed job information
- **job_tracking.json**: Complete job metadata

### Manual Monitoring

If you use `--no-wait`, monitor jobs manually:

```bash
# Check all jobs
bjobs

# Check specific jobs
bjobs <job_id>

# View logs
cat lsf_logs/nolzss_*.log
```

## Examples

### Example 1: Basic Submission

```bash
python -m noLZSS.genomics.lsf_batch_factorize \
    --file-list genomes.txt \
    --output-dir results
```

### Example 2: High Priority Queue

```bash
python -m noLZSS.genomics.lsf_batch_factorize \
    --file-list genomes.txt \
    --output-dir results \
    --queue priority \
    --safety-factor 2.0
```

### Example 3: Submit and Detach

```bash
python -m noLZSS.genomics.lsf_batch_factorize \
    --file-list genomes.txt \
    --output-dir results \
    --no-wait
```

### Example 4: Dry Run

```bash
python -m noLZSS.genomics.lsf_batch_factorize \
    --file-list genomes.txt \
    --output-dir results \
    --dry-run
```

### Example 5: Custom Logging

```bash
python -m noLZSS.genomics.lsf_batch_factorize \
    --file-list genomes.txt \
    --output-dir results \
    --log-level DEBUG \
    --log-file batch_run.log \
    --log-dir my_logs
```

## Troubleshooting

### Problem: Jobs fail with "Out of Memory"

**Solution**: Increase safety factor or check benchmark accuracy

```bash
--safety-factor 2.5
```

### Problem: Jobs timeout

**Solution**: Increase safety factor or verify time estimates

```bash
--safety-factor 2.0
```

### Problem: No trend parameters found

**Solution**: Run benchmarks or provide custom file

```bash
--trend-file /path/to/trend_parameters.pkl
```

Or the script will use conservative defaults.

### Problem: Jobs stuck in pending

**Solution**: Check cluster status and queue availability

```bash
bjobs -p  # Show pending jobs
bqueues   # Show queue status
```

### Problem: Import errors on compute nodes

**Solution**: Ensure noLZSS is installed in the environment used by LSF

```bash
# Test on a compute node
bsub -Is python -c "import noLZSS; print(noLZSS.__version__)"
```

## Best Practices

1. **Run Benchmarks First**: Generate accurate trend parameters for your cluster
2. **Use Appropriate Safety Factors**: 1.5-2.0 for most use cases
3. **Monitor First Few Jobs**: Verify resource estimates are accurate
4. **Consolidate Logs**: Use a dedicated log directory
5. **Use Dry-Run**: Preview specifications before large submissions
6. **Test Small Batches**: Start with a few files to verify setup

## Performance Tips

1. **Batch Size**: Submit 10-100 jobs at a time for optimal cluster utilization
2. **Queue Selection**: Use appropriate queues for different job sizes
3. **Safety Factor**: Balance between over-allocation and job failure risk
4. **File Organization**: Group similar-sized files together

## Integration with Existing Workflows

### With Snakemake

```python
rule lsf_factorize:
    input:
        files="file_list.txt"
    output:
        directory("results")
    shell:
        """
        python -m noLZSS.genomics.lsf_batch_factorize \
            --file-list {input.files} \
            --output-dir {output} \
            --no-wait
        """
```

### With Nextflow

```groovy
process factorize {
    executor 'lsf'
    
    script:
    """
    python -m noLZSS.genomics.lsf_batch_factorize \
        --file-list file_list.txt \
        --output-dir results \
        --no-wait
    """
}
```

## Comparison with batch_factorize.py

| Feature | batch_factorize.py | lsf_batch_factorize.py |
|---------|-------------------|------------------------|
| Execution | Local/single machine | LSF cluster |
| Resource Management | Manual | Automatic estimation |
| Parallelism | Python multiprocessing | LSF job scheduling |
| Monitoring | Built-in | LSF + script monitoring |
| Scalability | Limited by machine | Cluster-scale |
| Use Case | Small batches | Large batches |

## See Also

- [Batch Factorization Guide](batch_factorize.md)
- [Benchmarking Guide](../benchmarks/README.md)
- [Resource Predictor](../benchmarks/fasta_predictor.py)

## Support

For issues or questions:
1. Check logs in `lsf_logs/` directory
2. Review LSF documentation: `man bsub`, `man bjobs`
3. Open an issue on GitHub
