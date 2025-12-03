#!/usr/bin/env python3
"""
LSF Cluster batch FASTA file factorization script.

This script processes multiple FASTA files and factorizes each one on an LSF cluster
using bsub. It uses benchmarking results to estimate time, memory, and disk space
needed for each job, and optimally allocates threads and memory to avoid waste.

The script:
- Estimates resource requirements based on file sizes and benchmark trends
- Submits jobs to LSF cluster with optimal thread and memory allocation
- Tracks job completion and provides consolidated failure reports
- Avoids creating excessive logs or email notifications
"""

import argparse
import logging
import os
import sys
import subprocess
import tempfile
import time
import json
import pickle
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Dict, Union, Optional, Tuple, Any
from urllib.error import URLError, HTTPError

from ..utils import NoLZSSError
from .batch_factorize import (
    setup_logging,
    is_url,
    download_file,
    is_gzipped,
    decompress_gzip,
    read_file_list,
    FactorizationMode,
    BatchFactorizeError
)


class LSFBatchFactorizeError(NoLZSSError):
    """Raised when LSF batch factorization encounters an error."""
    pass


def get_file_size(file_path: Path) -> int:
    """
    Get the size of a file in bytes.
    
    Args:
        file_path: Path to the file
        
    Returns:
        File size in bytes
    """
    return file_path.stat().st_size


def estimate_fasta_nucleotides(file_path: Path) -> int:
    """
    Estimate the number of nucleotides in a FASTA file by sampling.
    
    Args:
        file_path: Path to the FASTA file
        
    Returns:
        Estimated number of nucleotides
    """
    # Read a sample to estimate overhead from headers
    sample_size = min(100000, file_path.stat().st_size)
    
    with open(file_path, 'r') as f:
        sample = f.read(sample_size)
    
    # Count headers and sequence characters
    lines = sample.split('\n')
    header_chars = sum(len(line) + 1 for line in lines if line.startswith('>'))
    total_chars = len(sample)
    sequence_chars = total_chars - header_chars
    
    # Estimate total nucleotides based on ratio
    if total_chars > 0:
        sequence_ratio = sequence_chars / total_chars
        file_size = file_path.stat().st_size
        estimated_nucleotides = int(file_size * sequence_ratio)
    else:
        # Fallback: assume 80% of file is sequence data
        estimated_nucleotides = int(file_path.stat().st_size * 0.8)
    
    return estimated_nucleotides


def load_benchmark_trends(trend_file: Optional[Path] = None) -> Optional[Dict[str, Any]]:
    """
    Load benchmark trends from file.
    
    Args:
        trend_file: Path to trend parameters file (default: benchmarks/fasta_results/trend_parameters.pkl)
        
    Returns:
        Trend parameters dictionary or None if file not found
    """
    if trend_file is None:
        # Try default locations
        possible_paths = [
            Path(__file__).parent.parent.parent.parent.parent / "benchmarks" / "fasta_results" / "trend_parameters.pkl",
            Path("benchmarks") / "fasta_results" / "trend_parameters.pkl",
            Path.cwd() / "benchmarks" / "fasta_results" / "trend_parameters.pkl"
        ]
        
        for path in possible_paths:
            if path.exists():
                trend_file = path
                break
    
    if trend_file is None or not trend_file.exists():
        return None
    
    try:
        if trend_file.suffix == '.pkl':
            with open(trend_file, 'rb') as f:
                return pickle.load(f)
        elif trend_file.suffix == '.json':
            with open(trend_file, 'r') as f:
                return json.load(f)
    except Exception as e:
        logging.warning(f"Failed to load benchmark trends from {trend_file}: {e}")
        return None


def estimate_resources_from_trends(
    file_size_nucleotides: int,
    trends: Dict[str, Any],
    mode: str,
    num_threads: int = 1,
    safety_factor: float = 1.5
) -> Dict[str, Any]:
    """
    Estimate resource requirements based on benchmark trends.
    
    Args:
        file_size_nucleotides: Input size in nucleotides
        trends: Benchmark trend parameters
        mode: Factorization mode (with_reverse_complement or without_reverse_complement)
        num_threads: Number of threads to use
        safety_factor: Safety factor for resource estimates
        
    Returns:
        Dictionary with resource estimates
    """
    # Import predictor functions
    import sys
    from pathlib import Path
    
    # Add benchmarks to path if not already there
    benchmarks_path = Path(__file__).parent.parent.parent.parent.parent / "benchmarks"
    if benchmarks_path.exists() and str(benchmarks_path) not in sys.path:
        sys.path.insert(0, str(benchmarks_path))
    
    try:
        from fasta_predictor import predict_from_trend
    except ImportError:
        # Fallback: use simple linear estimates
        def predict_from_trend(size, trend_params):
            if trend_params.get('log_scale', False):
                import numpy as np
                log_size = np.log10(size)
                log_prediction = trend_params['slope'] * log_size + trend_params['intercept']
                return 10 ** log_prediction
            else:
                return trend_params['slope'] * size + trend_params['intercept']
    
    # Select the appropriate function based on mode
    if mode == FactorizationMode.WITH_REVERSE_COMPLEMENT:
        func_name = "parallel_write_factors_binary_file_fasta_multiple_dna_w_rc"
    else:
        func_name = "parallel_write_factors_binary_file_fasta_multiple_dna_no_rc"
    
    # Get trend parameters for single-threaded execution (baseline)
    # Note: parallel benchmarks may not be in the trends file, so we use single-threaded as baseline
    # and apply a scaling factor for parallelism
    base_func_name = func_name.replace("parallel_", "")
    
    time_key = f"{base_func_name}_time"
    memory_key = f"{base_func_name}_memory"
    disk_key = f"{base_func_name}_disk_space"
    
    estimate = {
        'input_size_nucleotides': file_size_nucleotides,
        'input_size_mbp': file_size_nucleotides / 1_000_000,
        'num_threads': num_threads,
        'safety_factor': safety_factor
    }
    
    # Estimate time (parallel speedup is not linear, use Amdahl's law approximation)
    if time_key in trends:
        base_time_ms = predict_from_trend(file_size_nucleotides, trends[time_key])
        # Assume 80% of the work can be parallelized
        parallel_fraction = 0.8
        serial_fraction = 1 - parallel_fraction
        speedup = 1 / (serial_fraction + parallel_fraction / num_threads)
        parallel_time_ms = base_time_ms / speedup
        
        estimate['estimated_time_seconds'] = parallel_time_ms / 1000
        estimate['estimated_time_minutes'] = parallel_time_ms / 60000
        estimate['estimated_time_hours'] = parallel_time_ms / 3600000
        estimate['safe_time_seconds'] = (parallel_time_ms / 1000) * safety_factor
        estimate['safe_time_minutes'] = (parallel_time_ms / 60000) * safety_factor
        estimate['safe_time_hours'] = (parallel_time_ms / 3600000) * safety_factor
    
    # Estimate memory (scales with threads and data size)
    if memory_key in trends:
        base_memory_mb = predict_from_trend(file_size_nucleotides, trends[memory_key])
        # Memory increases with threads but not linearly (shared data structures)
        thread_memory_factor = 1 + 0.3 * (num_threads - 1)  # 30% increase per additional thread
        parallel_memory_mb = base_memory_mb * thread_memory_factor
        
        estimate['estimated_memory_mb'] = parallel_memory_mb
        estimate['estimated_memory_gb'] = parallel_memory_mb / 1024
        estimate['safe_memory_gb'] = (parallel_memory_mb / 1024) * safety_factor
        
        # Round up to common cluster memory allocations (powers of 2 or multiples of 4)
        safe_mem_gb = estimate['safe_memory_gb']
        if safe_mem_gb <= 2:
            cluster_mem = 2
        elif safe_mem_gb <= 4:
            cluster_mem = 4
        elif safe_mem_gb <= 8:
            cluster_mem = 8
        elif safe_mem_gb <= 16:
            cluster_mem = 16
        elif safe_mem_gb <= 32:
            cluster_mem = 32
        elif safe_mem_gb <= 64:
            cluster_mem = 64
        elif safe_mem_gb <= 128:
            cluster_mem = 128
        else:
            # Round up to nearest 64 GB
            cluster_mem = int((safe_mem_gb + 63) // 64 * 64)
        
        estimate['cluster_memory_gb'] = cluster_mem
    
    # Estimate disk space (doesn't change with parallelism)
    if disk_key in trends:
        disk_mb = predict_from_trend(file_size_nucleotides, trends[disk_key])
        estimate['estimated_disk_mb'] = disk_mb
        estimate['estimated_disk_gb'] = disk_mb / 1024
        estimate['safe_disk_gb'] = (disk_mb / 1024) * safety_factor
    
    return estimate


def estimate_resources_fallback(
    file_size_bytes: int,
    num_threads: int = 1,
    safety_factor: float = 1.5
) -> Dict[str, Any]:
    """
    Fallback resource estimation when benchmark trends are not available.
    
    Uses conservative estimates based on empirical observations.
    
    Args:
        file_size_bytes: Input file size in bytes
        num_threads: Number of threads to use
        safety_factor: Safety factor for resource estimates
        
    Returns:
        Dictionary with resource estimates
    """
    # Estimate nucleotides (assume ~1 byte per nucleotide with overhead)
    nucleotides = int(file_size_bytes * 0.8)
    
    # Conservative estimates based on typical LZSS behavior
    # Time: roughly 0.1-1 ms per 1000 nucleotides depending on complexity
    base_time_seconds = (nucleotides / 1000) * 0.0005  # 0.5 ms per 1000 nt
    
    # Apply parallel speedup (conservative estimate)
    speedup = num_threads * 0.7  # 70% efficiency
    parallel_time_seconds = base_time_seconds / speedup
    
    # Memory: CST requires ~20-40 bytes per input character
    base_memory_mb = (nucleotides * 30) / (1024 * 1024)  # 30 bytes per char
    thread_memory_factor = 1 + 0.3 * (num_threads - 1)
    parallel_memory_mb = base_memory_mb * thread_memory_factor
    
    # Disk: binary factors are typically 10-30% of input size
    disk_mb = (file_size_bytes * 0.2) / (1024 * 1024)
    
    estimate = {
        'input_size_nucleotides': nucleotides,
        'input_size_mbp': nucleotides / 1_000_000,
        'num_threads': num_threads,
        'safety_factor': safety_factor,
        'estimated_time_seconds': parallel_time_seconds,
        'estimated_time_minutes': parallel_time_seconds / 60,
        'estimated_time_hours': parallel_time_seconds / 3600,
        'safe_time_seconds': parallel_time_seconds * safety_factor,
        'safe_time_minutes': (parallel_time_seconds / 60) * safety_factor,
        'safe_time_hours': (parallel_time_seconds / 3600) * safety_factor,
        'estimated_memory_mb': parallel_memory_mb,
        'estimated_memory_gb': parallel_memory_mb / 1024,
        'safe_memory_gb': (parallel_memory_mb / 1024) * safety_factor,
        'estimated_disk_mb': disk_mb,
        'estimated_disk_gb': disk_mb / 1024,
        'safe_disk_gb': (disk_mb / 1024) * safety_factor,
    }
    
    # Round up memory
    safe_mem_gb = estimate['safe_memory_gb']
    if safe_mem_gb <= 2:
        cluster_mem = 2
    elif safe_mem_gb <= 4:
        cluster_mem = 4
    elif safe_mem_gb <= 8:
        cluster_mem = 8
    elif safe_mem_gb <= 16:
        cluster_mem = 16
    elif safe_mem_gb <= 32:
        cluster_mem = 32
    elif safe_mem_gb <= 64:
        cluster_mem = 64
    elif safe_mem_gb <= 128:
        cluster_mem = 128
    else:
        cluster_mem = int((safe_mem_gb + 63) // 64 * 64)
    
    estimate['cluster_memory_gb'] = cluster_mem
    
    return estimate


def decide_num_threads(
    file_size_nucleotides: int,
    max_threads: int,
    trends: Optional[Dict[str, Any]] = None
) -> int:
    """
    Decide optimal number of threads for a job based on file size and max_threads.
    
    Smaller files don't benefit from many threads due to overhead.
    
    Args:
        file_size_nucleotides: Input size in nucleotides
        max_threads: Maximum threads allowed by user
        trends: Optional benchmark trends for better decision
        
    Returns:
        Optimal number of threads
    """
    # Size thresholds (nucleotides)
    SMALL_FILE = 100_000      # 100kb
    MEDIUM_FILE = 1_000_000   # 1Mbp
    LARGE_FILE = 10_000_000   # 10Mbp
    
    if file_size_nucleotides < SMALL_FILE:
        # Small files: single thread is often faster
        optimal_threads = 1
    elif file_size_nucleotides < MEDIUM_FILE:
        # Medium files: use 2-4 threads
        optimal_threads = min(4, max_threads)
    elif file_size_nucleotides < LARGE_FILE:
        # Large files: use up to 8 threads
        optimal_threads = min(8, max_threads)
    else:
        # Very large files: use all available threads
        optimal_threads = max_threads
    
    return optimal_threads


def submit_lsf_job(
    job_name: str,
    script_path: Path,
    num_threads: int,
    memory_gb: int,
    time_minutes: int,
    output_log: Path,
    error_log: Path,
    queue: Optional[str] = None,
    extra_bsub_args: Optional[List[str]] = None,
    logger: Optional[logging.Logger] = None
) -> Optional[str]:
    """
    Submit a job to LSF cluster using bsub.
    
    Args:
        job_name: Name for the job
        script_path: Path to the script to execute
        num_threads: Number of threads/cores to request
        memory_gb: Memory in GB to request
        time_minutes: Time limit in minutes
        output_log: Path for stdout log
        error_log: Path for stderr log
        queue: Optional queue name
        extra_bsub_args: Optional additional bsub arguments
        logger: Logger instance
        
    Returns:
        Job ID if submission successful, None otherwise
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Build bsub command
    bsub_cmd = [
        'bsub',
        '-J', job_name,
        '-n', str(num_threads),
        '-R', f'rusage[mem={memory_gb * 1024}]',  # LSF typically uses MB
        '-R', f'span[hosts=1]',  # Keep all cores on one host
        '-W', str(time_minutes),
        '-o', str(output_log),
        '-e', str(error_log),
    ]
    
    # Add queue if specified
    if queue:
        bsub_cmd.extend(['-q', queue])
    
    # Add extra arguments if provided
    if extra_bsub_args:
        bsub_cmd.extend(extra_bsub_args)
    
    # Add the script to execute
    bsub_cmd.append(str(script_path))
    
    try:
        logger.debug(f"Submitting job with command: {' '.join(bsub_cmd)}")
        result = subprocess.run(
            bsub_cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Parse job ID from output
        # LSF typically returns: Job <12345> is submitted to queue <queuename>
        output = result.stdout.strip()
        logger.debug(f"bsub output: {output}")
        
        import re
        match = re.search(r'Job <(\d+)>', output)
        if match:
            job_id = match.group(1)
            logger.info(f"Successfully submitted job {job_name} with ID {job_id}")
            return job_id
        else:
            logger.warning(f"Could not parse job ID from bsub output: {output}")
            return None
            
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to submit job {job_name}: {e.stderr}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error submitting job {job_name}: {e}")
        return None


def create_job_script(
    input_file: Path,
    output_file: Path,
    mode: str,
    num_threads: int,
    script_path: Path,
    python_executable: Optional[str] = None
) -> bool:
    """
    Create a bash script for the LSF job.
    
    Args:
        input_file: Path to input FASTA file
        output_file: Path to output binary file
        mode: Factorization mode
        num_threads: Number of threads to use
        script_path: Path where to save the script
        python_executable: Optional Python executable path (default: sys.executable)
        
    Returns:
        True if script created successfully
    """
    if python_executable is None:
        python_executable = sys.executable
    
    # Determine which function to call
    if mode == FactorizationMode.WITH_REVERSE_COMPLEMENT:
        function_name = "parallel_write_factors_binary_file_fasta_multiple_dna_w_rc"
    elif mode == FactorizationMode.WITHOUT_REVERSE_COMPLEMENT:
        function_name = "parallel_write_factors_binary_file_fasta_multiple_dna_no_rc"
    else:
        return False
    
    # Create Python script content
    python_code = f"""#!/usr/bin/env python3
import sys
from noLZSS._noLZSS import {function_name}

input_file = "{input_file}"
output_file = "{output_file}"
num_threads = {num_threads}

try:
    factor_count = {function_name}(input_file, output_file, num_threads)
    print(f"Successfully factorized {{input_file}}: {{factor_count}} factors")
    sys.exit(0)
except Exception as e:
    print(f"Error factorizing {{input_file}}: {{e}}", file=sys.stderr)
    sys.exit(1)
"""
    
    # Create bash wrapper script
    bash_script = f"""#!/bin/bash
set -e
set -u

# Run the Python factorization
{python_executable} << 'PYTHON_EOF'
{python_code}
PYTHON_EOF

exit $?
"""
    
    try:
        script_path.parent.mkdir(parents=True, exist_ok=True)
        with open(script_path, 'w') as f:
            f.write(bash_script)
        
        # Make script executable
        script_path.chmod(0o755)
        return True
        
    except Exception as e:
        logging.error(f"Failed to create job script {script_path}: {e}")
        return False


def check_job_status(job_id: str, logger: Optional[logging.Logger] = None) -> str:
    """
    Check the status of an LSF job.
    
    Args:
        job_id: Job ID to check
        logger: Logger instance
        
    Returns:
        Job status string: 'PEND', 'RUN', 'DONE', 'EXIT', 'UNKNOWN'
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    try:
        result = subprocess.run(
            ['bjobs', '-noheader', '-o', 'stat', job_id],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if result.returncode == 0:
            status = result.stdout.strip()
            return status if status else 'UNKNOWN'
        else:
            # Job not found or bjobs error
            logger.debug(f"bjobs returned non-zero for job {job_id}: {result.stderr}")
            return 'DONE'  # Assume completed if not in queue
            
    except subprocess.TimeoutExpired:
        logger.warning(f"Timeout checking status for job {job_id}")
        return 'UNKNOWN'
    except Exception as e:
        logger.warning(f"Error checking status for job {job_id}: {e}")
        return 'UNKNOWN'


def wait_for_jobs(
    job_ids: Dict[str, str],
    check_interval: int = 60,
    logger: Optional[logging.Logger] = None
) -> Dict[str, str]:
    """
    Wait for all jobs to complete and return their final statuses.
    
    Args:
        job_ids: Dictionary mapping job names to job IDs
        check_interval: Interval between status checks in seconds
        logger: Logger instance
        
    Returns:
        Dictionary mapping job names to final status
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    logger.info(f"Waiting for {len(job_ids)} jobs to complete...")
    
    job_statuses = {name: 'PEND' for name in job_ids.keys()}
    completed_statuses = {'DONE', 'EXIT'}
    
    while True:
        all_done = True
        
        for job_name, job_id in job_ids.items():
            if job_statuses[job_name] not in completed_statuses:
                status = check_job_status(job_id, logger)
                
                if status != job_statuses[job_name]:
                    logger.info(f"Job {job_name} ({job_id}): {job_statuses[job_name]} -> {status}")
                    job_statuses[job_name] = status
                
                if status not in completed_statuses:
                    all_done = False
        
        if all_done:
            break
        
        # Count jobs in different states
        pend_count = sum(1 for s in job_statuses.values() if s == 'PEND')
        run_count = sum(1 for s in job_statuses.values() if s == 'RUN')
        done_count = sum(1 for s in job_statuses.values() if s == 'DONE')
        exit_count = sum(1 for s in job_statuses.values() if s == 'EXIT')
        
        logger.info(f"Job status: {pend_count} pending, {run_count} running, "
                   f"{done_count} completed, {exit_count} failed")
        
        time.sleep(check_interval)
    
    logger.info("All jobs completed")
    return job_statuses


def check_job_output(
    output_file: Path,
    error_log: Path,
    logger: Optional[logging.Logger] = None
) -> Tuple[bool, Optional[str]]:
    """
    Check if a job completed successfully by examining output and error logs.
    
    Args:
        output_file: Expected output binary file
        error_log: Path to error log
        logger: Logger instance
        
    Returns:
        Tuple of (success, error_message)
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Check if output file exists and is non-empty
    if not output_file.exists():
        return False, "Output file not created"
    
    if output_file.stat().st_size == 0:
        return False, "Output file is empty"
    
    # Check error log for errors
    if error_log.exists() and error_log.stat().st_size > 0:
        try:
            with open(error_log, 'r') as f:
                error_content = f.read()
            
            # Look for common error patterns
            if 'Error' in error_content or 'Exception' in error_content:
                # Extract first few lines of error
                error_lines = error_content.split('\n')[:5]
                error_msg = '\n'.join(error_lines)
                return False, f"Errors in log: {error_msg}"
        except Exception as e:
            logger.debug(f"Could not read error log {error_log}: {e}")
    
    return True, None


def compute_sequence_complexity_table_on_cluster(
    fasta_path: Union[str, Path],
    output_tsv: Union[str, Path],
    scripts_dir: Optional[Path] = None,
    logs_dir: Optional[Path] = None,
    queue: Optional[str] = None,
    check_interval: int = 60,
    extra_bsub_args: Optional[List[str]] = None,
    python_executable: Optional[str] = None,
    logger: Optional[logging.Logger] = None
) -> bool:
    """
    Compute sequence complexity table using LSF cluster.
    
    Parses FASTA and submits individual jobs for each sequence (both RC and no-RC counting)
    to maximize parallelization on the cluster. Combines results into a TSV file.
    
    Args:
        fasta_path: Path to the FASTA file
        output_tsv: Path to output TSV file
        scripts_dir: Directory for job scripts (default: temp directory)
        logs_dir: Directory for job logs (default: temp directory)
        queue: LSF queue name (optional)
        check_interval: Interval for checking job status (seconds)
        extra_bsub_args: Additional bsub arguments
        python_executable: Python executable path (default: sys.executable)
        logger: Logger instance
        
    Returns:
        True if successful, False otherwise
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    if python_executable is None:
        python_executable = sys.executable
    
    fasta_path = Path(fasta_path)
    output_tsv = Path(output_tsv)
    
    # Create temporary directories if not provided
    cleanup_temp = False
    if scripts_dir is None or logs_dir is None:
        temp_dir = Path(tempfile.mkdtemp(prefix="complexity_cluster_"))
        scripts_dir = temp_dir / "scripts"
        logs_dir = temp_dir / "logs"
        cleanup_temp = True
    
    scripts_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        # Parse FASTA file to get all sequences
        logger.info(f"Parsing FASTA file: {fasta_path}")
        with open(fasta_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        from .fasta import _parse_fasta_content
        sequences = _parse_fasta_content(content)
        
        logger.info(f"Found {len(sequences)} sequences")
        
        # Create temporary output directory for individual results
        temp_results_dir = output_tsv.parent / f"{output_tsv.stem}_temp_results"
        temp_results_dir.mkdir(parents=True, exist_ok=True)
        
        # Create and submit jobs for each sequence (both RC and no-RC)
        logger.info("Creating job scripts for all sequences...")
        job_ids = {}
        job_info = {}  # Maps job_name -> (seq_id, mode)
        
        for seq_idx, (seq_id, sequence) in enumerate(sequences.items()):
            # Job for counting with RC
            job_name_w_rc = f"complexity_w_rc_{seq_idx}"
            script_w_rc = scripts_dir / f"{job_name_w_rc}.sh"
            temp_output_w_rc = temp_results_dir / f"{job_name_w_rc}.json"
            
            python_code_w_rc = f"""#!/usr/bin/env python3
import sys
import json
from noLZSS._noLZSS import count_factors_dna_w_rc

seq_id = "{seq_id}"
sequence = "{sequence}"
output_file = "{temp_output_w_rc}"

try:
    seq_bytes = sequence.encode('ascii')
    count = count_factors_dna_w_rc(seq_bytes)
    result = {{'seq_id': seq_id, 'count': count}}
    with open(output_file, 'w') as f:
        json.dump(result, f)
    print(f"Successfully counted factors (w/ RC) for {{seq_id}}: {{count}} factors")
    sys.exit(0)
except Exception as e:
    print(f"Error counting factors (w/ RC) for {{seq_id}}: {{e}}", file=sys.stderr)
    sys.exit(1)
"""
            
            with open(script_w_rc, 'w') as f:
                f.write(f"""#!/bin/bash
set -e
set -u

{python_executable} << 'PYTHON_EOF'
{python_code_w_rc}
PYTHON_EOF

exit $?
""")
            script_w_rc.chmod(0o755)
            
            # Job for counting without RC
            job_name_no_rc = f"complexity_no_rc_{seq_idx}"
            script_no_rc = scripts_dir / f"{job_name_no_rc}.sh"
            temp_output_no_rc = temp_results_dir / f"{job_name_no_rc}.json"
            
            python_code_no_rc = f"""#!/usr/bin/env python3
import sys
import json
from noLZSS._noLZSS import count_factors

seq_id = "{seq_id}"
sequence = "{sequence}"
output_file = "{temp_output_no_rc}"

try:
    seq_bytes = sequence.encode('ascii')
    count = count_factors(seq_bytes)
    result = {{'seq_id': seq_id, 'count': count}}
    with open(output_file, 'w') as f:
        json.dump(result, f)
    print(f"Successfully counted factors (no RC) for {{seq_id}}: {{count}} factors")
    sys.exit(0)
except Exception as e:
    print(f"Error counting factors (no RC) for {{seq_id}}: {{e}}", file=sys.stderr)
    sys.exit(1)
"""
            
            with open(script_no_rc, 'w') as f:
                f.write(f"""#!/bin/bash
set -e
set -u

{python_executable} << 'PYTHON_EOF'
{python_code_no_rc}
PYTHON_EOF

exit $?
""")
            script_no_rc.chmod(0o755)
            
            # Submit both jobs
            job_id_w_rc = submit_lsf_job(
                job_name=job_name_w_rc,
                script_path=script_w_rc,
                num_threads=1,
                memory_gb=4,  # Per-sequence jobs need less memory
                time_minutes=30,
                output_log=logs_dir / f"{job_name_w_rc}.out",
                error_log=logs_dir / f"{job_name_w_rc}.err",
                queue=queue,
                extra_bsub_args=extra_bsub_args,
                logger=logger
            )
            
            job_id_no_rc = submit_lsf_job(
                job_name=job_name_no_rc,
                script_path=script_no_rc,
                num_threads=1,
                memory_gb=4,
                time_minutes=30,
                output_log=logs_dir / f"{job_name_no_rc}.out",
                error_log=logs_dir / f"{job_name_no_rc}.err",
                queue=queue,
                extra_bsub_args=extra_bsub_args,
                logger=logger
            )
            
            if job_id_w_rc:
                job_ids[job_name_w_rc] = job_id_w_rc
                job_info[job_name_w_rc] = (seq_id, 'w_rc', temp_output_w_rc)
            else:
                logger.warning(f"Failed to submit job for {seq_id} (w/ RC)")
            
            if job_id_no_rc:
                job_ids[job_name_no_rc] = job_id_no_rc
                job_info[job_name_no_rc] = (seq_id, 'no_rc', temp_output_no_rc)
            else:
                logger.warning(f"Failed to submit job for {seq_id} (no RC)")
        
        if not job_ids:
            logger.error("Failed to submit any jobs")
            return False
        
        logger.info(f"Submitted {len(job_ids)} jobs to LSF cluster")
        
        # Wait for all jobs to complete
        logger.info("Waiting for all jobs to complete...")
        job_statuses = wait_for_jobs(job_ids, check_interval, logger)
        
        # Check for failures
        failed_jobs = [name for name, status in job_statuses.items() if status != 'DONE']
        if failed_jobs:
            logger.error(f"{len(failed_jobs)} jobs failed: {failed_jobs}")
            return False
        
        # Read all results and combine
        logger.info("Reading results and creating TSV...")
        results_by_seq = {}  # seq_id -> {'w_rc': count, 'no_rc': count, 'length': length}
        
        for job_name, (seq_id, mode, temp_output) in job_info.items():
            if job_statuses.get(job_name) != 'DONE':
                continue
            
            try:
                with open(temp_output, 'r') as f:
                    result = json.load(f)
                
                if seq_id not in results_by_seq:
                    results_by_seq[seq_id] = {}
                results_by_seq[seq_id][mode] = result['count']
            except Exception as e:
                logger.error(f"Failed to read result for {seq_id} ({mode}): {e}")
                return False
        
        # Verify we have both counts for all sequences
        for seq_id in sequences.keys():
            if seq_id not in results_by_seq:
                logger.error(f"Missing results for sequence: {seq_id}")
                return False
            if 'w_rc' not in results_by_seq[seq_id] or 'no_rc' not in results_by_seq[seq_id]:
                logger.error(f"Incomplete results for sequence: {seq_id}")
                return False
        
        # Extract full headers from original content
        headers_map = {}  # seq_id -> full_header
        for line in content.split('\n'):
            if line.startswith('>'):
                full_header = line[1:].strip()  # Remove '>' and strip whitespace
                # Extract sequence_id (first word before space or tab)
                seq_id = full_header.split()[0] if full_header else full_header
                headers_map[seq_id] = full_header
        
        # Write combined TSV (preserve original sequence order)
        with open(output_tsv, 'w', encoding='utf-8') as f:
            f.write("sequence_id\theader\tlength\tcomplexity_w_rc\tcomplexity_no_rc\n")
            for seq_id, sequence in sequences.items():
                seq_length = len(sequence)
                full_header = headers_map.get(seq_id, seq_id)
                count_w_rc = results_by_seq[seq_id]['w_rc']
                count_no_rc = results_by_seq[seq_id]['no_rc']
                f.write(f"{seq_id}\t{full_header}\t{seq_length}\t{count_w_rc}\t{count_no_rc}\n")
        
        logger.info(f"Successfully created complexity table: {output_tsv}")
        logger.info(f"Total sequences: {len(sequences)}")
        
        # Clean up temp results
        import shutil
        shutil.rmtree(temp_results_dir, ignore_errors=True)
        
        return True
        
    except Exception as e:
        logger.error(f"Error computing complexity table on cluster: {e}")
        return False
    
    finally:
        if cleanup_temp and 'temp_dir' in locals():
            try:
                import shutil
                shutil.rmtree(temp_dir)
                logger.debug(f"Cleaned up temporary directory: {temp_dir}")
            except Exception:
                pass


def process_files_on_cluster(
    file_list: List[str],
    output_dir: Path,
    mode: str,
    max_threads: int,
    trends: Optional[Dict[str, Any]],
    download_dir: Optional[Path] = None,
    queue: Optional[str] = None,
    safety_factor: float = 1.5,
    check_interval: int = 60,
    extra_bsub_args: Optional[List[str]] = None,
    skip_existing: bool = True,
    logger: Optional[logging.Logger] = None
) -> Dict[str, Dict[str, Any]]:
    """
    Process files on LSF cluster by submitting jobs with optimal resource allocation.
    
    Args:
        file_list: List of file paths (must be local files accessible to cluster)
        output_dir: Base output directory
        mode: Factorization mode
        max_threads: Maximum threads per job
        trends: Benchmark trend parameters (optional)
        download_dir: Directory for downloaded/prepared files
        queue: LSF queue name (optional)
        safety_factor: Safety factor for resource estimates
        check_interval: Interval for checking job status (seconds)
        extra_bsub_args: Additional bsub arguments
        skip_existing: Skip files with existing output
        logger: Logger instance
        
    Returns:
        Dictionary with job results and statistics
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Create directories
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create subdirectories for scripts and logs
    scripts_dir = output_dir / "lsf_scripts"
    logs_dir = output_dir / "lsf_logs"
    scripts_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    
    # Prepare jobs
    jobs = []
    job_info = {}
    
    logger.info(f"Preparing {len(file_list)} jobs...")
    
    for file_path_str in file_list:
        file_path = Path(file_path_str)
        
        # Check if file exists
        if not file_path.exists():
            logger.error(f"File not found: {file_path}")
            continue
        
        # Determine output path
        base_name = file_path.stem
        if mode == FactorizationMode.WITH_REVERSE_COMPLEMENT:
            mode_dir = output_dir / "with_reverse_complement"
        elif mode == FactorizationMode.WITHOUT_REVERSE_COMPLEMENT:
            mode_dir = output_dir / "without_reverse_complement"
        else:  # BOTH mode
            logger.error(f"BOTH mode not supported for individual job submission. Use separate jobs.")
            continue
        
        mode_dir.mkdir(parents=True, exist_ok=True)
        output_file = mode_dir / f"{base_name}.bin"
        
        # Skip if output exists and skip_existing is True
        if skip_existing and output_file.exists():
            logger.info(f"Skipping {file_path.name} (output exists)")
            continue
        
        # Estimate nucleotides
        try:
            nucleotides = estimate_fasta_nucleotides(file_path)
        except Exception as e:
            logger.warning(f"Could not estimate nucleotides for {file_path}: {e}, using file size")
            nucleotides = int(file_path.stat().st_size * 0.8)
        
        # Decide number of threads
        num_threads = decide_num_threads(nucleotides, max_threads, trends)
        
        # Estimate resources
        if trends:
            estimates = estimate_resources_from_trends(
                nucleotides, trends, mode, num_threads, safety_factor
            )
        else:
            estimates = estimate_resources_fallback(
                file_path.stat().st_size, num_threads, safety_factor
            )
        
        # Round up time to minutes (minimum 5 minutes)
        time_minutes = max(5, int(estimates.get('safe_time_minutes', 60) + 1))
        memory_gb = estimates.get('cluster_memory_gb', 8)
        
        # Create job script
        job_name = f"factorize_{base_name}"
        script_path = scripts_dir / f"{job_name}.sh"
        
        if not create_job_script(file_path, output_file, mode, num_threads, script_path):
            logger.error(f"Failed to create script for {file_path}")
            continue
        
        # Prepare job info
        output_log = logs_dir / f"{job_name}.out"
        error_log = logs_dir / f"{job_name}.err"
        
        job = {
            'name': job_name,
            'input_file': file_path,
            'output_file': output_file,
            'script_path': script_path,
            'num_threads': num_threads,
            'memory_gb': memory_gb,
            'time_minutes': time_minutes,
            'output_log': output_log,
            'error_log': error_log,
            'estimates': estimates
        }
        
        jobs.append(job)
    
    if not jobs:
        logger.warning("No jobs to submit")
        return {}
    
    # Submit jobs
    logger.info(f"Submitting {len(jobs)} jobs to LSF cluster...")
    
    job_ids = {}
    failed_submissions = []
    
    for job in jobs:
        job_id = submit_lsf_job(
            job_name=job['name'],
            script_path=job['script_path'],
            num_threads=job['num_threads'],
            memory_gb=job['memory_gb'],
            time_minutes=job['time_minutes'],
            output_log=job['output_log'],
            error_log=job['error_log'],
            queue=queue,
            extra_bsub_args=extra_bsub_args,
            logger=logger
        )
        
        if job_id:
            job_ids[job['name']] = job_id
            job_info[job['name']] = job
        else:
            failed_submissions.append(job['name'])
            logger.error(f"Failed to submit job {job['name']}")
    
    if not job_ids:
        logger.error("No jobs were successfully submitted")
        return {'failed_submissions': failed_submissions}
    
    logger.info(f"Successfully submitted {len(job_ids)} jobs")
    
    # Wait for jobs to complete
    job_statuses = wait_for_jobs(job_ids, check_interval, logger)
    
    # Check results
    results = {
        'total_jobs': len(jobs),
        'submitted_jobs': len(job_ids),
        'failed_submissions': failed_submissions,
        'completed': [],
        'failed': [],
        'job_details': {}
    }
    
    for job_name, status in job_statuses.items():
        job = job_info[job_name]
        
        # Check output
        success, error_msg = check_job_output(
            job['output_file'],
            job['error_log'],
            logger
        )
        
        job_result = {
            'name': job_name,
            'input_file': str(job['input_file']),
            'output_file': str(job['output_file']),
            'status': status,
            'success': success,
            'error_message': error_msg,
            'num_threads': job['num_threads'],
            'memory_gb': job['memory_gb'],
            'time_minutes': job['time_minutes'],
            'estimates': job['estimates'],
            'error_log': str(job['error_log'])
        }
        
        results['job_details'][job_name] = job_result
        
        if success:
            results['completed'].append(job_name)
            logger.info(f"Job {job_name} completed successfully")
        else:
            results['failed'].append(job_name)
            logger.error(f"Job {job_name} failed: {error_msg}")
    
    return results


def print_summary(results: Dict[str, Any], logger: Optional[logging.Logger] = None):
    """
    Print a summary of job results.
    
    Args:
        results: Job results dictionary
        logger: Logger instance
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    logger.info("="*80)
    logger.info("LSF BATCH FACTORIZATION SUMMARY")
    logger.info("="*80)
    
    total = results.get('total_jobs', 0)
    submitted = results.get('submitted_jobs', 0)
    completed = len(results.get('completed', []))
    failed = len(results.get('failed', []))
    failed_submissions = len(results.get('failed_submissions', []))
    
    logger.info(f"Total jobs prepared: {total}")
    logger.info(f"Successfully submitted: {submitted}")
    logger.info(f"Failed to submit: {failed_submissions}")
    logger.info(f"Completed successfully: {completed}")
    logger.info(f"Failed: {failed}")
    
    if failed_submissions:
        logger.info("\nFailed to submit:")
        for name in results['failed_submissions']:
            logger.info(f"  - {name}")
    
    if failed > 0:
        logger.info("\nFailed jobs:")
        for job_name in results['failed']:
            job_detail = results['job_details'][job_name]
            logger.info(f"  - {job_name}: {job_detail['error_message']}")
            logger.info(f"    Input: {job_detail['input_file']}")
            logger.info(f"    Error log: {job_detail.get('error_log', 'N/A')}")
    
    # Resource statistics
    if 'job_details' in results and results['job_details']:
        logger.info("\nResource allocation summary:")
        
        total_threads = sum(j['num_threads'] for j in results['job_details'].values())
        total_memory = sum(j['memory_gb'] for j in results['job_details'].values())
        total_time_minutes = sum(j['time_minutes'] for j in results['job_details'].values())
        avg_threads = total_threads / len(results['job_details'])
        avg_memory = total_memory / len(results['job_details'])
        
        logger.info(f"  Average threads per job: {avg_threads:.1f}")
        logger.info(f"  Average memory per job: {avg_memory:.1f} GB")
        logger.info(f"  Total thread-hours requested: {total_threads * total_time_minutes / 60:.1f}")


def save_results(results: Dict[str, Any], output_dir: Path, logger: Optional[logging.Logger] = None):
    """
    Save job results to file.
    
    Args:
        results: Job results dictionary
        output_dir: Output directory
        logger: Logger instance
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    results_file = output_dir / "lsf_job_results.json"
    
    try:
        # Convert Path objects to strings for JSON serialization
        serializable_results = {}
        for key, value in results.items():
            if isinstance(value, dict):
                serializable_results[key] = {
                    k: str(v) if isinstance(v, Path) else v
                    for k, v in value.items()
                }
            else:
                serializable_results[key] = value
        
        with open(results_file, 'w') as f:
            json.dump(serializable_results, f, indent=2)
        
        logger.info(f"Results saved to {results_file}")
        
    except Exception as e:
        logger.error(f"Failed to save results: {e}")


def main():
    """Main entry point for the LSF batch factorization script."""
    parser = argparse.ArgumentParser(
        description="Batch factorize FASTA files on LSF cluster with optimal resource allocation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process files with automatic resource estimation
  python -m noLZSS.genomics.lsf_batch_factorize --file-list files.txt --output-dir results --max-threads 16 --mode with_reverse_complement
  
  # Use custom benchmark trends and queue
  python -m noLZSS.genomics.lsf_batch_factorize file1.fasta file2.fasta --output-dir results --max-threads 8 --trend-file benchmarks/trends.pkl --queue priority
  
  # Process with custom safety factor and check interval
  python -m noLZSS.genomics.lsf_batch_factorize --file-list files.txt --output-dir results --max-threads 12 --safety-factor 2.0 --check-interval 120

Note: All input files must be local paths accessible to the LSF cluster nodes.
Remote URLs are not supported in this script (download them first using batch_factorize.py).
        """
    )
    
    # Input specification
    parser.add_argument(
        "--file-list", type=Path,
        help="Text file containing list of FASTA file paths (one per line)"
    )
    parser.add_argument(
        "files", nargs="*",
        help="FASTA file paths to process (must be local paths)"
    )
    
    # Output configuration
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory for binary factorization results"
    )
    parser.add_argument(
        "--mode", 
        choices=[FactorizationMode.WITHOUT_REVERSE_COMPLEMENT, FactorizationMode.WITH_REVERSE_COMPLEMENT],
        default=FactorizationMode.WITH_REVERSE_COMPLEMENT,
        help="Factorization mode (default: with_reverse_complement). Note: 'both' mode not supported."
    )
    
    # Resource configuration
    parser.add_argument(
        "--max-threads", type=int, required=True,
        help="Maximum number of threads per job"
    )
    parser.add_argument(
        "--trend-file", type=Path,
        help="Path to benchmark trend parameters file (.pkl or .json). If not provided, will search default locations."
    )
    parser.add_argument(
        "--safety-factor", type=float, default=1.5,
        help="Safety factor for resource estimates (default: 1.5)"
    )
    
    # LSF configuration
    parser.add_argument(
        "--queue", type=str,
        help="LSF queue name (optional)"
    )
    parser.add_argument(
        "--bsub-args", nargs="+",
        help="Additional bsub arguments (e.g., --bsub-args -P project_name -G group_name)"
    )
    parser.add_argument(
        "--check-interval", type=int, default=60,
        help="Interval for checking job status in seconds (default: 60)"
    )
    
    # Processing options
    parser.add_argument(
        "--force", action="store_true",
        help="Overwrite existing output files (default: skip existing)"
    )
    
    # Logging configuration
    parser.add_argument(
        "--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO",
        help="Logging level (default: INFO)"
    )
    parser.add_argument(
        "--log-file", type=Path,
        help="Log file path (logs to console if not specified)"
    )
    
    args = parser.parse_args()
    
    # Set up logging
    logger = setup_logging(args.log_level, args.log_file)
    
    try:
        # Check if bsub is available
        try:
            subprocess.run(['which', 'bsub'], capture_output=True, check=True)
        except subprocess.CalledProcessError:
            raise LSFBatchFactorizeError(
                "bsub command not found. This script requires LSF to be installed and configured."
            )
        
        # Get file list
        if args.file_list and args.files:
            raise LSFBatchFactorizeError("Cannot specify both --file-list and individual files")
        elif args.file_list:
            file_list = read_file_list(args.file_list, logger)
        elif args.files:
            file_list = args.files
        else:
            raise LSFBatchFactorizeError("Must specify either --file-list or individual files")
        
        # Validate that all files are local
        for file_path in file_list:
            if is_url(file_path):
                raise LSFBatchFactorizeError(
                    f"Remote URLs not supported: {file_path}. "
                    "Please download files first using batch_factorize.py"
                )
            if not Path(file_path).exists():
                raise LSFBatchFactorizeError(f"File not found: {file_path}")
        
        logger.info(f"Starting LSF batch factorization of {len(file_list)} files")
        logger.info(f"Mode: {args.mode}")
        logger.info(f"Output directory: {args.output_dir}")
        logger.info(f"Max threads per job: {args.max_threads}")
        logger.info(f"Safety factor: {args.safety_factor}")
        
        # Load benchmark trends
        trends = load_benchmark_trends(args.trend_file)
        if trends:
            logger.info(f"Loaded benchmark trends from {args.trend_file or 'default location'}")
        else:
            logger.warning("Benchmark trends not available, using fallback estimates")
        
        # Process files
        results = process_files_on_cluster(
            file_list=file_list,
            output_dir=args.output_dir,
            mode=args.mode,
            max_threads=args.max_threads,
            trends=trends,
            queue=args.queue,
            safety_factor=args.safety_factor,
            check_interval=args.check_interval,
            extra_bsub_args=args.bsub_args,
            skip_existing=not args.force,
            logger=logger
        )
        
        # Save results
        save_results(results, args.output_dir, logger)
        
        # Print summary
        print_summary(results, logger)
        
        # Exit with appropriate code
        failed_count = len(results.get('failed', [])) + len(results.get('failed_submissions', []))
        if failed_count > 0:
            logger.warning(f"Completed with {failed_count} failures")
            sys.exit(1)
        else:
            logger.info("All jobs completed successfully")
            sys.exit(0)
        
    except LSFBatchFactorizeError as e:
        logger.error(f"LSF batch factorization error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        logger.info("Interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
