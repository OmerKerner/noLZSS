# LSF Batch Factorizer

A script for processing multiple FASTA files on LSF (Platform Load Sharing Facility) clusters with optimal resource allocation based on benchmarking results.

## Overview

The `lsf_batch_factorize.py` script extends the functionality of `batch_factorize.py` to work with LSF cluster environments. It:

1. **Estimates resource requirements** (time, memory, disk) based on file sizes and benchmark trends
2. **Optimally allocates threads** per job to avoid wasting resources
3. **Submits jobs** to LSF using `bsub` with appropriate resource requests
4. **Tracks job completion** and provides consolidated failure reports
5. **Minimizes log clutter** - organized logs without email spam

## Prerequisites

- LSF cluster with `bsub` and `bjobs` commands available
- noLZSS package installed (`pip install -e .`)
- Local FASTA files accessible to cluster nodes
- Optional: Benchmark trend parameters file (`.pkl` or `.json`)

## Basic Usage

### Process files with automatic resource estimation:
```bash
python -m noLZSS.genomics.lsf_batch_factorize \
    --file-list files.txt \
    --output-dir results \
    --max-threads 16 \
    --mode with_reverse_complement
```

### Process individual files:
```bash
python -m noLZSS.genomics.lsf_batch_factorize \
    file1.fasta file2.fasta file3.fasta \
    --output-dir results \
    --max-threads 8 \
    --mode with_reverse_complement
```

### Use custom benchmark trends and queue:
```bash
python -m noLZSS.genomics.lsf_batch_factorize \
    --file-list files.txt \
    --output-dir results \
    --max-threads 12 \
    --trend-file benchmarks/fasta_results/trend_parameters.pkl \
    --queue priority
```

### Specify additional LSF options:
```bash
python -m noLZSS.genomics.lsf_batch_factorize \
    --file-list files.txt \
    --output-dir results \
    --max-threads 16 \
    --bsub-args -P project_name -G group_name
```

## Command-Line Options

### Input Files
- `--file-list FILE_LIST`: Text file with one FASTA path per line
- `files`: Positional arguments for individual FASTA files

**Note**: All files must be local paths accessible to cluster nodes. Remote URLs are not supported (download them first using `batch_factorize.py`).

### Output Configuration
- `--output-dir OUTPUT_DIR`: Base directory for output files (required)
- `--mode MODE`: Factorization mode (default: `with_reverse_complement`)
  - `with_reverse_complement`: DNA factorization with reverse complement awareness
  - `without_reverse_complement`: DNA factorization without reverse complement

**Note**: The `both` mode is not supported as it would require submitting separate jobs for each mode.

### Resource Configuration
- `--max-threads N`: Maximum threads per job (required)
- `--trend-file FILE`: Benchmark trend parameters file (`.pkl` or `.json`)
- `--safety-factor FACTOR`: Safety factor for resource estimates (default: 1.5)

### LSF Configuration
- `--queue QUEUE`: LSF queue name
- `--bsub-args ARG [ARG ...]`: Additional bsub arguments
- `--check-interval N`: Job status check interval in seconds (default: 60)

### Processing Options
- `--force`: Overwrite existing output files (default: skip existing)
- `--log-level LEVEL`: Logging level (DEBUG, INFO, WARNING, ERROR)
- `--log-file FILE`: Log file path (default: console only)

## Resource Estimation

### With Benchmark Trends
If benchmark trends are available (from `benchmarks/fasta_results/trend_parameters.pkl`), the script uses them to accurately predict:
- **Time**: Based on input size with parallel speedup calculation
- **Memory**: Scaled for multi-threading with shared data structures
- **Disk space**: Based on expected factor output size

### Without Benchmark Trends
Falls back to conservative estimates:
- Time: ~0.5 ms per 1000 nucleotides
- Memory: ~30 bytes per nucleotide
- Disk: ~20% of input size

### Thread Allocation Strategy
The script automatically decides optimal thread count based on file size:

| File Size | Threads |
|-----------|---------|
| < 100 kb | 1 |
| < 1 Mbp | 2-4 |
| < 10 Mbp | up to 8 |
| ≥ 10 Mbp | max_threads |

This prevents overhead on small files while maximizing parallelism for large files.

### Memory Allocation
Memory requests are rounded to common cluster allocations: 2, 4, 8, 16, 32, 64, 128 GB.

## Output Structure

```
output-dir/
├── with_reverse_complement/   # or without_reverse_complement/
│   ├── file1.bin
│   ├── file2.bin
│   └── ...
├── lsf_scripts/
│   ├── factorize_file1.sh
│   ├── factorize_file2.sh
│   └── ...
├── lsf_logs/
│   ├── factorize_file1.out
│   ├── factorize_file1.err
│   ├── factorize_file2.out
│   ├── factorize_file2.err
│   └── ...
└── lsf_job_results.json
```

### Output Files
- **Binary factor files**: `{mode}/{basename}.bin` - factorization results
- **Job scripts**: `lsf_scripts/factorize_{basename}.sh` - submitted scripts
- **Job logs**: `lsf_logs/factorize_{basename}.{out,err}` - stdout/stderr
- **Results JSON**: `lsf_job_results.json` - job status and statistics

## Job Tracking

The script monitors jobs and provides real-time updates:

```
Job status: 5 pending, 10 running, 8 completed, 0 failed
Job factorize_chr1 (12345): PEND -> RUN
Job factorize_chr2 (12346): RUN -> DONE
...
```

After completion, it checks output files and error logs to determine actual success/failure.

## Error Handling

### Submission Failures
Jobs that fail to submit are reported immediately:
```
Failed to submit job factorize_file1
```

### Job Failures
Failed jobs are summarized at the end with error messages:
```
Failed jobs:
  - factorize_file1: Output file not created
    Input: /path/to/file1.fasta
    Error log: results/lsf_logs/factorize_file1.err
```

### Consolidated Reporting
- All job information saved to JSON for analysis
- Minimal log output to avoid clutter
- No email notifications per job (rely on LSF's own settings)

## Running Benchmark to Generate Trends

To generate benchmark trends for accurate resource estimation:

```bash
# Run FASTA benchmark
python benchmarks/fasta_benchmark.py --output-dir benchmarks/fasta_results

# This creates:
# - benchmarks/fasta_results/trend_parameters.pkl
# - benchmarks/fasta_results/trend_parameters.json
```

Then use the trends file:
```bash
python -m noLZSS.genomics.lsf_batch_factorize \
    --file-list files.txt \
    --output-dir results \
    --max-threads 16 \
    --trend-file benchmarks/fasta_results/trend_parameters.pkl
```

## Comparison with batch_factorize.py

| Feature | batch_factorize.py | lsf_batch_factorize.py |
|---------|-------------------|------------------------|
| Execution | Local multi-threading | LSF cluster jobs |
| Remote files | Yes (downloads) | No (must be local) |
| Resource estimation | No | Yes (from benchmarks) |
| Thread allocation | Fixed | Dynamic per file |
| Job tracking | N/A | Yes (bjobs) |
| Suitable for | Local machines | HPC clusters |

## Example Workflow

1. **Prepare file list**:
```bash
echo "/data/genome1.fasta" > files.txt
echo "/data/genome2.fasta" >> files.txt
echo "/data/genome3.fasta" >> files.txt
```

2. **Generate benchmarks** (optional but recommended):
```bash
python benchmarks/fasta_benchmark.py --output-dir benchmarks/fasta_results
```

3. **Submit jobs**:
```bash
python -m noLZSS.genomics.lsf_batch_factorize \
    --file-list files.txt \
    --output-dir results \
    --max-threads 16 \
    --trend-file benchmarks/fasta_results/trend_parameters.pkl \
    --queue normal \
    --log-file lsf_batch.log
```

4. **Monitor progress**: The script will wait and report status updates.

5. **Check results**:
```bash
# View summary
cat results/lsf_job_results.json

# Check individual logs if needed
less results/lsf_logs/factorize_genome1.err
```

## Troubleshooting

### "bsub command not found"
- Ensure LSF is installed and configured
- Check that `bsub` is in your PATH

### "File not found" errors
- Verify all input files exist and are accessible to cluster nodes
- Use absolute paths in file list

### Jobs timing out
- Increase safety factor: `--safety-factor 2.0`
- Check benchmark trends are accurate for your data

### High memory usage
- Reduce max threads for smaller cluster nodes
- Check if files are larger than expected

### Jobs stuck in PEND
- Check queue limits and availability
- Verify resource requests are within queue limits

## Advanced Usage

### Process only new files
```bash
# Skip files with existing output
python -m noLZSS.genomics.lsf_batch_factorize \
    --file-list files.txt \
    --output-dir results \
    --max-threads 16
```

### Reprocess all files
```bash
# Overwrite existing output
python -m noLZSS.genomics.lsf_batch_factorize \
    --file-list files.txt \
    --output-dir results \
    --max-threads 16 \
    --force
```

### Debug mode
```bash
# Enable debug logging
python -m noLZSS.genomics.lsf_batch_factorize \
    --file-list files.txt \
    --output-dir results \
    --max-threads 16 \
    --log-level DEBUG \
    --log-file debug.log
```

## See Also

- `batch_factorize.py` - Local batch processing with file download support
- `fasta_benchmark.py` - Generate benchmark trends
- `fasta_predictor.py` - Predict resource requirements
- `parallel_benchmark.py` - Benchmark parallel functions
