"""
Batch FASTA file factorization script with support for local and remote files.

This script processes multiple FASTA files and factorizes each one using optimized
C++ functions that handle validation and binary output. It supports parallel downloads
and factorization for improved performance. Also supports per-sequence complexity
TSV generation mode.
"""

import argparse
import gzip
import logging
import os
import random
import sys
import tempfile
import time
import urllib.request
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import List, Dict, Union, Optional, Tuple
from urllib.error import URLError, HTTPError

from ..utils import NoLZSSError
from .._noLZSS import (
    write_factors_binary_file_fasta_multiple_dna_w_rc,
    write_factors_binary_file_fasta_multiple_dna_no_rc,
    count_factors,
    count_factors_dna_w_rc,
)
from .fasta import _parse_fasta_content


class BatchFactorizeError(NoLZSSError):
    """Raised when batch factorization encounters an error."""
    pass


class FactorizationMode:
    """Enumeration of factorization modes."""
    WITHOUT_REVERSE_COMPLEMENT = "without_reverse_complement"
    WITH_REVERSE_COMPLEMENT = "with_reverse_complement"
    BOTH = "both"


def setup_logging(log_level: str = "INFO", log_file: Optional[Path] = None) -> logging.Logger:
    """
    Set up logging for the batch factorization process.
    
    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        log_file: Optional log file path
        
    Returns:
        Configured logger instance
    """
    logger = logging.getLogger("batch_factorize")
    logger.setLevel(getattr(logging, log_level.upper()))
    
    # Remove existing handlers to avoid duplication
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_format = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)
    
    # File handler if specified
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_format = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s'
        )
        file_handler.setFormatter(file_format)
        logger.addHandler(file_handler)
    
    return logger


def is_url(path: str) -> bool:
    """
    Check if a path is a URL.
    
    Args:
        path: Path string to check
        
    Returns:
        True if path appears to be a URL
    """
    return path.startswith(('http://', 'https://', 'ftp://'))


def is_gzipped(file_path: Path) -> bool:
    """
    Check if a file is gzipped by reading the first few bytes.
    
    Args:
        file_path: Path to the file to check
        
    Returns:
        True if the file appears to be gzipped
    """
    try:
        with open(file_path, 'rb') as f:
            # Gzip files start with 0x1f 0x8b
            magic_bytes = f.read(2)
            return magic_bytes == b'\x1f\x8b'
    except (OSError, IOError):
        return False


def decompress_gzip(input_path: Path, output_path: Path, logger: Optional[logging.Logger] = None) -> bool:
    """
    Decompress a gzipped file.
    
    Args:
        input_path: Path to the gzipped input file
        output_path: Path where to save the decompressed file
        logger: Logger instance for progress reporting
        
    Returns:
        True if decompression successful, False otherwise
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    try:
        logger.info(f"Decompressing {input_path} to {output_path}")
        
        # Create output directory if it doesn't exist
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with gzip.open(input_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                # Read in chunks to handle large files
                chunk_size = 8192
                while True:
                    chunk = f_in.read(chunk_size)
                    if not chunk:
                        break
                    f_out.write(chunk)
        
        logger.info(f"Successfully decompressed {input_path}")
        return True
        
    except Exception as e:
        logger.error(f"Failed to decompress {input_path}: {e}")
        return False


def download_file(url: str, output_path: Path, max_retries: int = 3, 
                 timeout: int = 30, logger: Optional[logging.Logger] = None) -> bool:
    """
    Download a file from URL with retry logic.
    
    Args:
        url: URL to download from
        output_path: Local path to save the file
        max_retries: Maximum number of retry attempts
        timeout: Download timeout in seconds
        logger: Logger instance for progress reporting
        
    Returns:
        True if download successful, False otherwise
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    for attempt in range(max_retries):
        try:
            logger.info(f"Downloading {url} (attempt {attempt + 1}/{max_retries})")
            
            # Create directory if it doesn't exist
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Download with timeout
            with urllib.request.urlopen(url, timeout=timeout) as response:
                with open(output_path, 'wb') as f:
                    # Read in chunks to handle large files
                    chunk_size = 8192
                    while True:
                        chunk = response.read(chunk_size)
                        if not chunk:
                            break
                        f.write(chunk)
            
            logger.info(f"Successfully downloaded {url} to {output_path}")
            return True
            
        except (URLError, HTTPError, OSError, TimeoutError) as e:
            logger.warning(f"Download attempt {attempt + 1} failed for {url}: {e}")
            if attempt < max_retries - 1:
                wait_time = (attempt + 1) * 2  # Exponential backoff
                logger.info(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                logger.error(f"Failed to download {url} after {max_retries} attempts")
    
    return False


def shuffle_fasta_sequences(input_path: Path, output_path: Path, 
                            seed: Optional[int] = None,
                            logger: Optional[logging.Logger] = None) -> bool:
    """
    Create a shuffled version of a FASTA file.
    
    Randomizes each sequence separately while maintaining original headers.
    This is useful for creating control datasets for statistical analysis.
    
    Args:
        input_path: Path to input FASTA file
        output_path: Path to output shuffled FASTA file
        seed: Random seed for reproducibility (optional)
        logger: Logger instance for progress reporting
        
    Returns:
        True if shuffling successful, False otherwise
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    try:
        logger.info(f"Creating shuffled version of {input_path}")
        
        # Read the FASTA file
        with open(input_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Parse FASTA content
        sequences = _parse_fasta_content(content)
        
        if not sequences:
            logger.error(f"No sequences found in {input_path}")
            return False
        
        # Set random seed if provided
        if seed is not None:
            random.seed(seed)
        
        # Create output directory if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Write shuffled sequences
        with open(output_path, 'w', encoding='utf-8') as f:
            for seq_id, sequence in sequences.items():
                # Shuffle the sequence
                seq_list = list(sequence)
                random.shuffle(seq_list)
                shuffled_seq = ''.join(seq_list)
                
                # Write with original header
                f.write(f">{seq_id}\n")
                # Write in lines of 80 characters (standard FASTA format)
                for i in range(0, len(shuffled_seq), 80):
                    f.write(shuffled_seq[i:i+80] + '\n')
        
        logger.info(f"Successfully created shuffled FASTA at {output_path}")
        return True
        
    except Exception as e:
        logger.error(f"Failed to shuffle {input_path}: {e}")
        return False


def plot_factor_comparison(original_factors_file: Path, shuffled_factors_file: Path,
                          output_plot_path: Path, 
                          original_label: str = "Original",
                          shuffled_label: str = "Shuffled Control",
                          logger: Optional[logging.Logger] = None) -> bool:
    """
    Create a comparison plot showing original factors vs shuffled control.
    
    Creates a plot with:
    - Solid line for original sequence factorization
    - Dotted black line for shuffled sequence control
    
    Args:
        original_factors_file: Path to binary factors file for original sequence
        shuffled_factors_file: Path to binary factors file for shuffled sequence
        output_plot_path: Path to save the comparison plot
        original_label: Label for original data (default: "Original")
        shuffled_label: Label for shuffled data (default: "Shuffled Control")
        logger: Logger instance for progress reporting
        
    Returns:
        True if plot created successfully, False otherwise
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib is required for plotting. Skipping plot generation.")
        return False
    
    from ..utils import read_factors_binary_file
    
    try:
        logger.info(f"Creating comparison plot: {output_plot_path}")
        
        # Read factors from binary files
        original_factors = read_factors_binary_file(original_factors_file)
        shuffled_factors = read_factors_binary_file(shuffled_factors_file)
        
        if not original_factors:
            logger.error(f"No factors found in original file: {original_factors_file}")
            return False
        
        if not shuffled_factors:
            logger.error(f"No factors found in shuffled file: {shuffled_factors_file}")
            return False
        
        # Compute cumulative lengths for both
        def compute_cumulative(factors):
            cumulative = []
            current_sum = 0
            for i, factor in enumerate(factors):
                length = factor[1]  # (position, length, ref)
                current_sum += length
                cumulative.append((i + 1, current_sum))
            return cumulative
        
        original_cumulative = compute_cumulative(original_factors)
        shuffled_cumulative = compute_cumulative(shuffled_factors)
        
        # Extract x and y
        orig_y, orig_x = zip(*original_cumulative) if original_cumulative else ([], [])
        shuf_y, shuf_x = zip(*shuffled_cumulative) if shuffled_cumulative else ([], [])
        
        # Create plot
        plt.figure(figsize=(10, 6))
        
        # Plot original (solid line)
        plt.step(orig_x, orig_y, where='post', linewidth=1.5, label=original_label, color='blue')
        plt.plot(orig_x, orig_y, linestyle='', marker='o', markersize=3, alpha=0.6, color='blue')
        
        # Plot shuffled control (dotted black line)
        plt.step(shuf_x, shuf_y, where='post', linewidth=1.5, linestyle=':', 
                label=shuffled_label, color='black')
        plt.plot(shuf_x, shuf_y, linestyle='', marker='s', markersize=3, alpha=0.4, color='black')
        
        plt.xlabel('Cumulative Factor Length')
        plt.ylabel('Factor Index')
        plt.title('Factor Length Accumulation: Original vs Shuffled Control')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Save plot
        output_plot_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Successfully created comparison plot: {output_plot_path}")
        return True
        
    except Exception as e:
        logger.error(f"Failed to create comparison plot: {e}")
        return False


def _count_sequence_factors(seq_info: Tuple[str, str, str]) -> Tuple[str, str, int, int, int]:
    """
    Count factors for a single sequence (both with and without reverse complement).
    
    This is a module-level function to support multiprocessing (must be picklable).
    
    Args:
        seq_info: Tuple of (sequence_id, full_header, sequence)
        
    Returns:
        Tuple of (sequence_id, full_header, length, count_w_rc, count_no_rc)
    """
    seq_id, full_header, sequence = seq_info
    seq_bytes = sequence.encode('ascii')
    seq_length = len(sequence)
    count_w_rc = count_factors_dna_w_rc(seq_bytes)
    count_no_rc = count_factors(seq_bytes)
    return (seq_id, full_header, seq_length, count_w_rc, count_no_rc)


def compute_sequence_complexity_table(fasta_path: Union[str, Path], 
                                     num_processes: Optional[int] = None) -> List[Tuple[str, str, int, int, int]]:
    """
    Compute per-sequence complexity table with both RC and no-RC factor counts.
    
    Uses Python multiprocessing to process all sequences in parallel.
    Each sequence is counted both with and without reverse complement.
    
    Args:
        fasta_path: Path to the FASTA file
        num_processes: Number of processes to use (None = use CPU count)
        
    Returns:
        List of tuples: (sequence_id, full_header, length, complexity_w_rc, complexity_no_rc)
    """
    fasta_path = Path(fasta_path)
    
    # Read and parse FASTA file and extract full headers
    with open(fasta_path, 'r', encoding='utf-8') as f:
        content = f.read()
    sequences = _parse_fasta_content(content)
    
    # Extract full headers from original content
    headers_map = {}  # seq_id -> full_header
    for line in content.split('\n'):
        if line.startswith('>'):
            full_header = line[1:].strip()  # Remove '>' and strip whitespace
            # Extract sequence_id (first word before space or tab)
            seq_id = full_header.split()[0] if full_header else full_header
            headers_map[seq_id] = full_header
    
    # Prepare input for parallel processing: (seq_id, full_header, sequence)
    processing_input = [(seq_id, headers_map.get(seq_id, seq_id), sequence) 
                       for seq_id, sequence in sequences.items()]
    
    # Process all sequences in parallel using module-level function
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        results = list(executor.map(_count_sequence_factors, processing_input))
    
    return results


def write_sequence_complexity_tsv(fasta_path: Union[str, Path], 
                                  output_path: Union[str, Path],
                                  num_processes: Optional[int] = None) -> int:
    """
    Write per-sequence complexity table to TSV file.
    
    Args:
        fasta_path: Path to the FASTA file
        output_path: Path to output TSV file
        num_processes: Number of processes to use (None = use 2)
        
    Returns:
        Number of sequences written
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Compute complexity table
    rows = compute_sequence_complexity_table(fasta_path, num_processes)
    
    # Write to TSV
    with open(output_path, 'w', encoding='utf-8') as f:
        # Header
        f.write("sequence_id\theader\tlength\tcomplexity_w_rc\tcomplexity_no_rc\n")
        
        # Data rows
        for seq_id, full_header, length, count_w_rc, count_no_rc in rows:
            f.write(f"{seq_id}\t{full_header}\t{length}\t{count_w_rc}\t{count_no_rc}\n")
    
    return len(rows)


def get_output_paths(input_path: Path, output_dir: Path, mode: str) -> Dict[str, Path]:
    """
    Generate output file paths based on input file and mode.
    
    Args:
        input_path: Input FASTA file path
        output_dir: Base output directory
        mode: Factorization mode
        
    Returns:
        Dictionary mapping mode names to output file paths
    """
    base_name = input_path.stem  # File name without extension
    
    paths = {}
    if mode in [FactorizationMode.WITHOUT_REVERSE_COMPLEMENT, FactorizationMode.BOTH]:
        without_rc_dir = output_dir / "without_reverse_complement"
        without_rc_dir.mkdir(parents=True, exist_ok=True)
        paths["without_reverse_complement"] = without_rc_dir / f"{base_name}.bin"
    
    if mode in [FactorizationMode.WITH_REVERSE_COMPLEMENT, FactorizationMode.BOTH]:
        with_rc_dir = output_dir / "with_reverse_complement"
        with_rc_dir.mkdir(parents=True, exist_ok=True)
        paths["with_reverse_complement"] = with_rc_dir / f"{base_name}.bin"
    
    return paths


def factorize_single_file(input_path: Path, output_paths: Dict[str, Path],
                         skip_existing: bool = True, 
                         logger: Optional[logging.Logger] = None) -> Dict[str, bool]:
    """
    Factorize a single FASTA file with specified modes using optimized C++ functions.
    
    Args:
        input_path: Path to input FASTA file
        output_paths: Dictionary mapping mode names to output paths
        skip_existing: Whether to skip if output already exists
        logger: Logger instance
        
    Returns:
        Dictionary mapping mode names to success status
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    results = {}
    
    # Define a helper function for factorizing a single mode
    def factorize_mode(mode: str, output_path: Path) -> Tuple[str, bool]:
        try:
            # Check if output already exists and skip if requested
            if skip_existing and output_path.exists():
                logger.info(f"Skipping {mode} factorization for {input_path.name} "
                           f"(output already exists: {output_path})")
                return mode, True
            
            logger.info(f"Starting {mode} factorization for {input_path.name}")
            
            # Create output directory
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            if mode == "without_reverse_complement":
                # Factorization without reverse complement
                factor_count = write_factors_binary_file_fasta_multiple_dna_no_rc(
                    str(input_path), str(output_path)
                )
            elif mode == "with_reverse_complement":
                # Factorization with reverse complement awareness
                factor_count = write_factors_binary_file_fasta_multiple_dna_w_rc(
                    str(input_path), str(output_path)
                )
            
            logger.info(f"Successfully completed {mode} factorization for {input_path.name} "
                       f"({factor_count} factors)")
            return mode, True
            
        except Exception as e:
            logger.error(f"Failed {mode} factorization for {input_path.name}: {e}")
            # Clean up partial output file
            if output_path.exists():
                try:
                    output_path.unlink()
                    logger.debug(f"Cleaned up partial output file: {output_path}")
                except OSError:
                    pass
            return mode, False
    
    # Use ThreadPoolExecutor to parallelize mode processing
    with ThreadPoolExecutor(max_workers=len(output_paths)) as executor:
        futures = {executor.submit(factorize_mode, mode, output_path): mode 
                  for mode, output_path in output_paths.items()}
        
        for future in as_completed(futures):
            mode = futures[future]
            try:
                mode_name, success = future.result()
                results[mode_name] = success
            except Exception as e:
                logger.error(f"Unexpected error in {mode} factorization: {e}")
                results[mode] = False
    
    return results


def download_file_worker(file_info: Tuple[str, Path, int, str]) -> Tuple[str, bool, Optional[Path]]:
    """
    Download a single file. This function is used for parallel processing.
    
    Args:
        file_info: Tuple of (file_path_or_url, download_dir, max_retries, logger_name)
        
    Returns:
        Tuple of (original_path, success, local_path)
    """
    file_path, download_dir, max_retries, logger_name = file_info
    
    # Create a logger for this process
    logger = logging.getLogger(logger_name)
    
    if is_url(file_path):
        # Download remote file
        file_name = Path(urllib.parse.urlparse(file_path).path).name
        if not file_name:
            file_name = f"downloaded_{hash(file_path) % 10000}.fasta"
        
        local_path = download_dir / file_name
        
        if not download_file(file_path, local_path, max_retries=max_retries, logger=logger):
            logger.error(f"Failed to download {file_path}")
            return file_path, False, None
    else:
        # Local file
        local_path = Path(file_path)
        if not local_path.exists():
            logger.error(f"Local file not found: {file_path}")
            return file_path, False, None
    
    # Check if file is gzipped and decompress if needed
    if is_gzipped(local_path):
        logger.info(f"Detected gzipped file: {local_path}")
        decompressed_path = local_path.with_suffix('')  # Remove .gz extension if present
        if decompressed_path.suffix == '.gz':
            decompressed_path = decompressed_path.with_suffix('')
        
        # If decompressed file already exists, use it
        if decompressed_path.exists():
            logger.info(f"Decompressed file already exists: {decompressed_path}")
            return file_path, True, decompressed_path
        
        # Decompress the file
        if decompress_gzip(local_path, decompressed_path, logger):
            # Clean up the compressed file if it was downloaded
            if is_url(file_path):
                try:
                    local_path.unlink()
                    logger.debug(f"Cleaned up compressed file: {local_path}")
                except OSError:
                    pass
            return file_path, True, decompressed_path
        else:
            logger.error(f"Failed to decompress {local_path}")
            return file_path, False, None
    
    return file_path, True, local_path


def factorize_file_worker(job_info: Tuple[str, Path, Dict[str, Path], bool, str]) -> Tuple[str, Dict[str, bool]]:
    """
    Worker function for parallel factorization.
    
    Args:
        job_info: Tuple of (original_path, input_path, output_paths, skip_existing, logger_name)
        
    Returns:
        Tuple of (original_path, factorization_results)
    """
    original_path, input_path, output_paths, skip_existing, logger_name = job_info
    
    # Create a logger for this process
    logger = logging.getLogger(logger_name)
    
    # Factorize the file
    factorization_results = factorize_single_file(
        input_path, output_paths, skip_existing=skip_existing, logger=logger
    )
    
    return original_path, factorization_results


def process_file_list(file_list: List[str], output_dir: Path, mode: str,
                     download_dir: Optional[Path] = None, skip_existing: bool = True,
                     max_retries: int = 3, max_workers: Optional[int] = None, 
                     logger: Optional[logging.Logger] = None) -> Dict[str, Dict[str, bool]]:
    """
    Process a list of FASTA files (local or remote) with parallel download and factorization.
    
    Args:
        file_list: List of file paths or URLs
        output_dir: Base output directory
        mode: Factorization mode
        download_dir: Directory for downloaded files (uses temp if None)
        skip_existing: Whether to skip existing output files
        max_retries: Maximum download retry attempts
        max_workers: Maximum number of worker threads/processes (None = auto)
        logger: Logger instance
        
    Returns:
        Dictionary mapping file names to their processing results
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    results = {}
    decompressed_files = []  # Track files that were decompressed for cleanup
    
    # Use provided download directory or create temp directory
    if download_dir is None:
        download_dir = Path(tempfile.mkdtemp(prefix="batch_factorize_"))
        cleanup_temp = True
        logger.info(f"Using temporary download directory: {download_dir}")
    else:
        download_dir.mkdir(parents=True, exist_ok=True)
        cleanup_temp = False
    
    try:
        # Step 1: Parallel download
        logger.info(f"Starting parallel download of {len(file_list)} files")
        
        download_jobs = []
        for file_path in file_list:
            download_jobs.append((file_path, download_dir, max_retries, logger.name))
        
        prepared_files = []
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_path = {executor.submit(download_file_worker, job): job[0] 
                            for job in download_jobs}
            
            for future in as_completed(future_to_path):
                original_path = future_to_path[future]
                try:
                    file_path, success, local_path = future.result()
                    if success and local_path:
                        prepared_files.append((file_path, local_path))
                        logger.info(f"Successfully prepared {file_path}")
                        
                        # Check if this file was decompressed (different from original download path)
                        original_download_path = download_dir / Path(urllib.parse.urlparse(file_path).path).name
                        if not is_url(file_path):
                            original_download_path = Path(file_path)
                        
                        if local_path != original_download_path and local_path.exists():
                            decompressed_files.append(local_path)
                        
                    else:
                        # Determine error type
                        if not success:
                            if is_url(file_path):
                                results[file_path] = {"error": "download_failed"}
                            else:
                                results[file_path] = {"error": "file_not_found"}
                        
                except Exception as e:
                    logger.error(f"Unexpected error processing {original_path}: {e}")
                    results[original_path] = {"error": "processing_error"}
        
        logger.info(f"Download complete: {len(prepared_files)} files ready for factorization")
        
        # Step 2: Parallel factorization
        if prepared_files:
            logger.info(f"Starting parallel factorization of {len(prepared_files)} files")
            
            factorization_jobs = []
            for original_path, local_path in prepared_files:
                output_paths = get_output_paths(local_path, output_dir, mode)
                factorization_jobs.append((original_path, local_path, output_paths, skip_existing, logger.name))
            
            # Use ProcessPoolExecutor for CPU-intensive factorization work
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                future_to_path = {executor.submit(factorize_file_worker, job): job[0] 
                                for job in factorization_jobs}
                
                for future in as_completed(future_to_path):
                    original_path = future_to_path[future]
                    try:
                        file_path, factorization_results = future.result()
                        results[file_path] = factorization_results
                        logger.info(f"Completed factorization for {file_path}")
                        
                    except Exception as e:
                        logger.error(f"Factorization error for {original_path}: {e}")
                        results[original_path] = {"error": "factorization_error"}
        
    finally:
        # Clean up decompressed files
        for decompressed_file in decompressed_files:
            try:
                if decompressed_file.exists():
                    decompressed_file.unlink()
                    logger.debug(f"Cleaned up decompressed file: {decompressed_file}")
            except OSError:
                logger.warning(f"Failed to clean up decompressed file: {decompressed_file}")
        
        # Clean up temporary download directory if we created it
        if cleanup_temp:
            try:
                import shutil
                shutil.rmtree(download_dir)
                logger.debug(f"Cleaned up temporary directory: {download_dir}")
            except OSError:
                logger.warning(f"Failed to clean up temporary directory: {download_dir}")
    
    return results


def read_file_list(list_file: Path, logger: Optional[logging.Logger] = None) -> List[str]:
    """
    Read a list of file paths/URLs from a text file.
    
    Args:
        list_file: Path to file containing list of paths/URLs
        logger: Logger instance
        
    Returns:
        List of file paths/URLs
        
    Raises:
        BatchFactorizeError: If file cannot be read or is empty
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    try:
        with open(list_file, 'r') as f:
            lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        
        if not lines:
            raise BatchFactorizeError(f"No valid file paths found in {list_file}")
        
        logger.info(f"Read {len(lines)} file paths from {list_file}")
        return lines
        
    except IOError as e:
        raise BatchFactorizeError(f"Failed to read file list {list_file}: {e}")


def print_summary(results: Dict[str, Dict[str, bool]], logger: Optional[logging.Logger] = None):
    """
    Print a summary of processing results.
    
    Args:
        results: Processing results dictionary
        logger: Logger instance
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    total_files = len(results)
    successful_files = 0
    failed_files = 0
    skipped_files = 0
    
    mode_stats = {}
    
    for file_path, file_results in results.items():
        if "error" in file_results:
            failed_files += 1
            continue
        
        file_success = True
        for mode, success in file_results.items():
            if mode not in mode_stats:
                mode_stats[mode] = {"success": 0, "failed": 0}
            
            if success:
                mode_stats[mode]["success"] += 1
            else:
                mode_stats[mode]["failed"] += 1
                file_success = False
        
        if file_success:
            successful_files += 1
        else:
            failed_files += 1
    
    logger.info("="*60)
    logger.info("BATCH FACTORIZATION SUMMARY")
    logger.info("="*60)
    logger.info(f"Total files: {total_files}")
    logger.info(f"Successfully processed: {successful_files}")
    logger.info(f"Failed: {failed_files}")
    
    if mode_stats:
        logger.info("\nMode-specific results:")
        for mode, stats in mode_stats.items():
            logger.info(f"  {mode}: {stats['success']} successful, {stats['failed']} failed")
    
    if failed_files > 0:
        logger.info("\nFailed files:")
        for file_path, file_results in results.items():
            if "error" in file_results:
                logger.info(f"  {file_path}: {file_results['error']}")
            elif not all(file_results.values()):
                failed_modes = [mode for mode, success in file_results.items() if not success]
                logger.info(f"  {file_path}: failed modes: {failed_modes}")


def process_with_shuffle_analysis(file_list: List[str], output_dir: Path, mode: str,
                                  download_dir: Optional[Path] = None, skip_existing: bool = True,
                                  max_retries: int = 3, max_workers: Optional[int] = None,
                                  shuffle_seed: Optional[int] = None,
                                  logger: Optional[logging.Logger] = None) -> Dict[str, Dict[str, any]]:
    """
    Process files with shuffle analysis: create shuffled versions, factorize both, and prepare for plotting.
    
    Args:
        file_list: List of file paths or URLs
        output_dir: Base output directory
        mode: Factorization mode
        download_dir: Directory for downloaded files (uses temp if None)
        skip_existing: Whether to skip existing output files
        max_retries: Maximum download retry attempts
        max_workers: Maximum number of worker threads/processes (None = auto)
        shuffle_seed: Random seed for shuffling (for reproducibility)
        logger: Logger instance
        
    Returns:
        Dictionary with original and shuffled factorization results
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Create a subdirectory for shuffled files
    shuffled_dir = output_dir / "shuffled_sequences"
    shuffled_dir.mkdir(parents=True, exist_ok=True)
    
    # First, process the original files
    logger.info("Processing original FASTA files...")
    original_results = process_file_list(
        file_list=file_list,
        output_dir=output_dir,
        mode=mode,
        download_dir=download_dir,
        skip_existing=skip_existing,
        max_retries=max_retries,
        max_workers=max_workers,
        logger=logger
    )
    
    # Now create and process shuffled versions
    logger.info("Creating and processing shuffled versions...")
    shuffled_file_list = []
    shuffled_file_mapping = {}  # Maps shuffled path to original path
    
    # Use temporary directory for shuffled files
    temp_shuffle_dir = Path(tempfile.mkdtemp(prefix="shuffled_fasta_"))
    
    try:
        # Download original files if needed and create shuffled versions
        for file_path in file_list:
            # Determine local path for the original file
            if is_url(file_path):
                # For URLs, we need to download first
                file_name = Path(urllib.parse.urlparse(file_path).path).name
                if download_dir:
                    local_path = download_dir / file_name
                else:
                    # Download to temp directory
                    local_path = temp_shuffle_dir / file_name
                    if not local_path.exists():
                        success = download_file(file_path, local_path, max_retries=max_retries, logger=logger)
                        if not success:
                            logger.warning(f"Skipping shuffle for {file_path} - download failed")
                            continue
            else:
                local_path = Path(file_path)
                if not local_path.exists():
                    logger.warning(f"Skipping shuffle for {file_path} - file not found")
                    continue
            
            # Check if file is gzipped and decompress if needed
            if is_gzipped(local_path):
                logger.info(f"Detected gzipped file for shuffling: {local_path}")
                decompressed_path = local_path.with_suffix('')  # Remove .gz extension if present
                if decompressed_path.suffix == '.gz':
                    decompressed_path = decompressed_path.with_suffix('')
                
                # Use temp directory for decompressed file
                decompressed_path = temp_shuffle_dir / decompressed_path.name
                
                # Decompress the file
                if not decompress_gzip(local_path, decompressed_path, logger):
                    logger.warning(f"Skipping shuffle for {file_path} - failed to decompress")
                    continue
                
                # Use decompressed file for shuffling
                local_path = decompressed_path
            
            # Create shuffled version
            base_name = local_path.stem
            shuffled_path = temp_shuffle_dir / f"{base_name}_shuffled.fasta"
            
            # Check if shuffled version already exists
            if skip_existing and shuffled_path.exists():
                logger.info(f"Shuffled version already exists: {shuffled_path}")
            else:
                success = shuffle_fasta_sequences(local_path, shuffled_path, seed=shuffle_seed, logger=logger)
                if not success:
                    logger.warning(f"Failed to create shuffled version of {file_path}")
                    continue
            
            shuffled_file_list.append(str(shuffled_path))
            shuffled_file_mapping[str(shuffled_path)] = file_path
        
        # Process shuffled files
        if shuffled_file_list:
            # Create separate output directory for shuffled factorizations
            shuffled_output_dir = output_dir / "shuffled"
            shuffled_output_dir.mkdir(parents=True, exist_ok=True)
            
            logger.info(f"Factorizing {len(shuffled_file_list)} shuffled files...")
            shuffled_results = process_file_list(
                file_list=shuffled_file_list,
                output_dir=shuffled_output_dir,
                mode=mode,
                download_dir=None,  # Files are already local
                skip_existing=skip_existing,
                max_retries=max_retries,
                max_workers=max_workers,
                logger=logger
            )
        else:
            logger.warning("No shuffled files to process")
            shuffled_results = {}
        
        # Generate comparison plots
        if shuffled_results:
            logger.info("Generating comparison plots...")
            plots_dir = output_dir / "comparison_plots"
            plots_dir.mkdir(parents=True, exist_ok=True)
            
            # For each original file that has a shuffled counterpart, create a comparison plot
            for shuffled_path, original_path in shuffled_file_mapping.items():
                # Get the base name to find the corresponding output files
                shuffled_local_path = Path(shuffled_path)
                base_name = shuffled_local_path.stem.replace('_shuffled', '')
                
                # Construct paths to the binary factor files
                # We need to check each mode that was processed
                for mode_name in ['with_reverse_complement', 'without_reverse_complement']:
                    original_bin_path = output_dir / mode_name / f"{base_name}.bin"
                    shuffled_bin_path = shuffled_output_dir / mode_name / f"{shuffled_local_path.stem}.bin"
                    
                    if original_bin_path.exists() and shuffled_bin_path.exists():
                        plot_path = plots_dir / f"{base_name}_{mode_name}_comparison.png"
                        
                        plot_label = mode_name.replace('_', ' ').title()
                        success = plot_factor_comparison(
                            original_factors_file=original_bin_path,
                            shuffled_factors_file=shuffled_bin_path,
                            output_plot_path=plot_path,
                            original_label=f"Original ({plot_label})",
                            shuffled_label="Shuffled Control",
                            logger=logger
                        )
                        
                        if success:
                            logger.info(f"Created comparison plot: {plot_path}")
                        else:
                            logger.warning(f"Failed to create comparison plot for {base_name} ({mode_name})")
    
    finally:
        # Clean up temporary shuffle directory
        try:
            import shutil
            shutil.rmtree(temp_shuffle_dir)
            logger.debug(f"Cleaned up temporary shuffle directory: {temp_shuffle_dir}")
        except OSError:
            logger.warning(f"Failed to clean up temporary directory: {temp_shuffle_dir}")
    
    # Combine results
    combined_results = {
        'original': original_results,
        'shuffled': shuffled_results,
        'mapping': shuffled_file_mapping
    }
    
    return combined_results


def complexity_file_worker(job_info: Tuple[str, Path, Path, Optional[int], str]) -> Tuple[str, Dict[str, any]]:
    """
    Worker function for parallel complexity TSV generation.
    
    Args:
        job_info: Tuple of (original_path, input_path, output_path, num_processes, logger_name)
        
    Returns:
        Tuple of (original_path, result_dict)
    """
    original_path, input_path, output_path, num_processes, logger_name = job_info
    
    # Create a logger for this process
    logger = logging.getLogger(logger_name)
    
    try:
        logger.info(f"Computing complexity TSV for {input_path.name}")
        
        # Generate complexity TSV
        num_sequences = write_sequence_complexity_tsv(
            fasta_path=input_path,
            output_path=output_path,
            num_processes=num_processes
        )
        
        logger.info(f"Successfully generated complexity TSV for {input_path.name} ({num_sequences} sequences)")
        return original_path, {"success": True, "num_sequences": num_sequences}
        
    except Exception as e:
        logger.error(f"Failed to generate complexity TSV for {input_path.name}: {e}")
        # Clean up partial output file
        if output_path.exists():
            try:
                output_path.unlink()
                logger.debug(f"Cleaned up partial output file: {output_path}")
            except OSError:
                pass
        return original_path, {"error": str(e)}


def process_file_list_complexity(file_list: List[str], output_dir: Path,
                                download_dir: Optional[Path] = None, skip_existing: bool = True,
                                max_retries: int = 3, max_workers: Optional[int] = None,
                                num_processes: Optional[int] = None,
                                logger: Optional[logging.Logger] = None) -> Dict[str, Dict[str, any]]:
    """
    Process a list of FASTA files to generate per-sequence complexity TSV files.
    
    Args:
        file_list: List of file paths or URLs
        output_dir: Base output directory for TSV files
        download_dir: Directory for downloaded files (uses temp if None)
        skip_existing: Whether to skip existing output files
        max_retries: Maximum download retry attempts
        max_workers: Maximum number of worker processes (None = auto)
        num_processes: Number of processes for complexity computation per file
        logger: Logger instance
        
    Returns:
        Dictionary mapping file names to their processing results
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    results = {}
    decompressed_files = []  # Track files that were decompressed for cleanup
    
    # Create complexity output directory
    complexity_dir = output_dir / "complexity"
    complexity_dir.mkdir(parents=True, exist_ok=True)
    
    # Use provided download directory or create temp directory
    if download_dir is None:
        download_dir = Path(tempfile.mkdtemp(prefix="batch_complexity_"))
        cleanup_temp = True
        logger.info(f"Using temporary download directory: {download_dir}")
    else:
        download_dir.mkdir(parents=True, exist_ok=True)
        cleanup_temp = False
    
    try:
        # Step 1: Parallel download
        logger.info(f"Starting parallel download of {len(file_list)} files")
        
        download_jobs = []
        for file_path in file_list:
            download_jobs.append((file_path, download_dir, max_retries, logger.name))
        
        prepared_files = []
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_path = {executor.submit(download_file_worker, job): job[0] 
                            for job in download_jobs}
            
            for future in as_completed(future_to_path):
                original_path = future_to_path[future]
                try:
                    file_path, success, local_path = future.result()
                    if success and local_path:
                        prepared_files.append((file_path, local_path))
                        logger.info(f"Successfully prepared {file_path}")
                        
                        # Check if this file was decompressed (different from original download path)
                        original_download_path = download_dir / Path(urllib.parse.urlparse(file_path).path).name
                        if not is_url(file_path):
                            original_download_path = Path(file_path)
                        
                        if local_path != original_download_path and local_path.exists():
                            decompressed_files.append(local_path)
                        
                    else:
                        # Determine error type
                        if not success:
                            if is_url(file_path):
                                results[file_path] = {"error": "download_failed"}
                            else:
                                results[file_path] = {"error": "file_not_found"}
                        
                except Exception as e:
                    logger.error(f"Unexpected error processing {original_path}: {e}")
                    results[original_path] = {"error": "processing_error"}
        
        logger.info(f"Download complete: {len(prepared_files)} files ready for complexity analysis")
        
        # Step 2: Parallel complexity TSV generation
        if prepared_files:
            logger.info(f"Starting parallel complexity analysis of {len(prepared_files)} files")
            
            complexity_jobs = []
            for original_path, local_path in prepared_files:
                base_name = local_path.stem
                output_path = complexity_dir / f"{base_name}_complexity.tsv"
                
                # Check if output already exists and skip if requested
                if skip_existing and output_path.exists():
                    logger.info(f"Skipping {original_path} (output already exists: {output_path})")
                    results[original_path] = {"success": True, "skipped": True}
                    continue
                
                complexity_jobs.append((original_path, local_path, output_path, num_processes, logger.name))
            
            # Use ProcessPoolExecutor for CPU-intensive complexity computation
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                future_to_path = {executor.submit(complexity_file_worker, job): job[0] 
                                for job in complexity_jobs}
                
                for future in as_completed(future_to_path):
                    original_path = future_to_path[future]
                    try:
                        file_path, complexity_result = future.result()
                        results[file_path] = complexity_result
                        logger.info(f"Completed complexity analysis for {file_path}")
                        
                    except Exception as e:
                        logger.error(f"Complexity analysis error for {original_path}: {e}")
                        results[original_path] = {"error": "complexity_error"}
        
    finally:
        # Clean up decompressed files
        for decompressed_file in decompressed_files:
            try:
                if decompressed_file.exists():
                    decompressed_file.unlink()
                    logger.debug(f"Cleaned up decompressed file: {decompressed_file}")
            except OSError:
                logger.warning(f"Failed to clean up decompressed file: {decompressed_file}")
        
        # Clean up temporary download directory if we created it
        if cleanup_temp:
            try:
                import shutil
                shutil.rmtree(download_dir)
                logger.debug(f"Cleaned up temporary directory: {download_dir}")
            except OSError:
                logger.warning(f"Failed to clean up temporary directory: {download_dir}")
    
    return results


def main():
    """Main entry point for the batch factorization script."""
    parser = argparse.ArgumentParser(
        description="Batch factorize FASTA files with support for local and remote files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process files listed in a text file with both modes
  python -m noLZSS.genomics.batch_factorize --file-list files.txt --output-dir results --mode both
  
  # Process individual files with reverse complement only and 4 parallel workers
  python -m noLZSS.genomics.batch_factorize file1.fasta file2.fasta --output-dir results --mode with_reverse_complement --max-workers 4
  
  # Process remote files with custom download directory
  python -m noLZSS.genomics.batch_factorize --file-list urls.txt --output-dir results --download-dir downloads --mode without_reverse_complement
  
  # Process files with shuffle analysis (creates shuffled controls and comparison plots)
  python -m noLZSS.genomics.batch_factorize file.fasta --output-dir results --mode with_reverse_complement --shuffle-analysis --shuffle-seed 42
  
  # Generate per-sequence complexity TSV files
  python -m noLZSS.genomics.batch_factorize --file-list files.txt --output-dir results --complexity-tsv --complexity-processes 4
        """
    )
    
    # Input specification
    parser.add_argument(
        "--file-list", type=Path,
        help="Text file containing list of FASTA file paths/URLs (one per line)"
    )
    parser.add_argument(
        "files", nargs="*",
        help="FASTA file paths/URLs to process"
    )
    
    # Output configuration
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory for binary factorization results"
    )
    parser.add_argument(
        "--mode", choices=[FactorizationMode.WITHOUT_REVERSE_COMPLEMENT, FactorizationMode.WITH_REVERSE_COMPLEMENT, FactorizationMode.BOTH],
        default=FactorizationMode.BOTH,
        help="Factorization mode (default: both)"
    )
    
    # Download configuration
    parser.add_argument(
        "--download-dir", type=Path,
        help="Directory for downloaded files (uses temp directory if not specified)"
    )
    parser.add_argument(
        "--max-retries", type=int, default=3,
        help="Maximum download retry attempts (default: 3)"
    )
    
    # Processing options
    parser.add_argument(
        "--force", action="store_true",
        help="Overwrite existing output files (default: skip existing)"
    )
    parser.add_argument(
        "--max-workers", type=int,
        help="Maximum number of parallel workers for downloads and factorization (default: CPU count)"
    )
    parser.add_argument(
        "--shuffle-analysis", action="store_true",
        help="Create shuffled versions of sequences as controls, factorize them, and include in plots"
    )
    parser.add_argument(
        "--shuffle-seed", type=int, default=None,
        help="Random seed for shuffling (for reproducibility)"
    )
    parser.add_argument(
        "--complexity-tsv", action="store_true",
        help="Generate per-sequence complexity TSV files instead of binary factorization"
    )
    parser.add_argument(
        "--complexity-processes", type=int, default=None,
        help="Number of processes for complexity computation per file (default: CPU count)"
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
        # Get file list
        if args.file_list and args.files:
            raise BatchFactorizeError("Cannot specify both --file-list and individual files")
        elif args.file_list:
            file_list = read_file_list(args.file_list, logger)
        elif args.files:
            file_list = args.files
        else:
            raise BatchFactorizeError("Must specify either --file-list or individual files")
        
        # Check for mutually exclusive modes
        if args.complexity_tsv and args.shuffle_analysis:
            raise BatchFactorizeError("Cannot specify both --complexity-tsv and --shuffle-analysis")
        
        logger.info(f"Starting batch processing of {len(file_list)} files")
        logger.info(f"Output directory: {args.output_dir}")
        
        # Process files based on mode
        if args.complexity_tsv:
            logger.info("Mode: Complexity TSV generation")
            if args.complexity_processes:
                logger.info(f"Complexity processes per file: {args.complexity_processes}")
            
            results = process_file_list_complexity(
                file_list=file_list,
                output_dir=args.output_dir,
                download_dir=args.download_dir,
                skip_existing=not args.force,
                max_retries=args.max_retries,
                max_workers=args.max_workers,
                num_processes=args.complexity_processes,
                logger=logger
            )
            
            # Print summary
            total_files = len(results)
            successful_files = sum(1 for r in results.values() if r.get("success", False))
            skipped_files = sum(1 for r in results.values() if r.get("skipped", False))
            failed_files = sum(1 for r in results.values() if "error" in r)
            
            logger.info("="*60)
            logger.info("COMPLEXITY TSV GENERATION SUMMARY")
            logger.info("="*60)
            logger.info(f"Total files: {total_files}")
            logger.info(f"Successfully processed: {successful_files}")
            logger.info(f"Skipped (already exist): {skipped_files}")
            logger.info(f"Failed: {failed_files}")
            
            if failed_files > 0:
                logger.info("\nFailed files:")
                for file_path, result in results.items():
                    if "error" in result:
                        logger.info(f"  {file_path}: {result['error']}")
            
            # Exit with appropriate code
            if failed_files > 0:
                logger.warning(f"Completed with {failed_files} failures")
                sys.exit(1)
            else:
                logger.info("All files processed successfully")
                sys.exit(0)
        
        elif args.shuffle_analysis:
            logger.info(f"Mode: Factorization with shuffle analysis ({args.mode})")
            combined_results = process_with_shuffle_analysis(
                file_list=file_list,
                output_dir=args.output_dir,
                mode=args.mode,
                download_dir=args.download_dir,
                skip_existing=not args.force,
                max_retries=args.max_retries,
                max_workers=args.max_workers,
                shuffle_seed=args.shuffle_seed,
                logger=logger
            )
            
            # Print summaries for both original and shuffled
            logger.info("\n=== ORIGINAL FILES SUMMARY ===")
            print_summary(combined_results['original'], logger)
            
            logger.info("\n=== SHUFFLED FILES SUMMARY ===")
            print_summary(combined_results['shuffled'], logger)
            
            # Exit with appropriate code
            original_failed = sum(1 for r in combined_results['original'].values() 
                                 if "error" in r or not all(r.values()))
            shuffled_failed = sum(1 for r in combined_results['shuffled'].values() 
                                 if "error" in r or not all(r.values()))
            total_failed = original_failed + shuffled_failed
            
            if total_failed > 0:
                logger.warning(f"Completed with {total_failed} failures "
                             f"({original_failed} original, {shuffled_failed} shuffled)")
                sys.exit(1)
            else:
                logger.info("All files (original and shuffled) processed successfully")
                sys.exit(0)
        
        else:
            logger.info(f"Mode: Factorization ({args.mode})")
            results = process_file_list(
                file_list=file_list,
                output_dir=args.output_dir,
                mode=args.mode,
                download_dir=args.download_dir,
                skip_existing=not args.force,
                max_retries=args.max_retries,
                max_workers=args.max_workers,
                logger=logger
            )
            
            # Print summary
            print_summary(results, logger)
            
            # Exit with appropriate code
            failed_count = sum(1 for r in results.values() if "error" in r or not all(r.values()))
            if failed_count > 0:
                logger.warning(f"Completed with {failed_count} failures")
                sys.exit(1)
            else:
                logger.info("All files processed successfully")
                sys.exit(0)
            
    except BatchFactorizeError as e:
        logger.error(f"Batch factorization error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        logger.info("Interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()