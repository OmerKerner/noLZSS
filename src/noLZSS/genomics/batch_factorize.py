"""
Batch FASTA file factorization script with support for local and remote files.

This script processes multiple FASTA files and factorizes each one using optimized
C++ functions that handle validation and binary output. It supports parallel downloads
and factorization for improved performance.
"""

import argparse
import gzip
import logging
import os
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
)


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
        
        logger.info(f"Starting batch factorization of {len(file_list)} files")
        logger.info(f"Mode: {args.mode}")
        logger.info(f"Output directory: {args.output_dir}")
        
        # Process files
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