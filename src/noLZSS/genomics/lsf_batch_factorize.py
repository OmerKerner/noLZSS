#!/usr/bin/env python3
"""
LSF batch FASTA file factorization script with intelligent resource estimation.

This script processes multiple FASTA files on LSF clusters with optimized resource
allocation based on benchmarking results. It estimates time, memory, and disk space
requirements for each job, monitors job completion, and provides failure summaries.
"""

import argparse
import json
import logging
import os
import pickle
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import tempfile

# Import the predictor and batch utilities
try:
    from ..utils import NoLZSSError
    from . import batch_factorize
except ImportError:
    # Fallback for development/testing
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from utils import NoLZSSError
    from genomics import batch_factorize


class LSFBatchError(NoLZSSError):
    """Raised when LSF batch submission encounters an error."""
    pass


class ResourceEstimator:
    """Estimates computational resources based on benchmark trends."""
    
    def __init__(self, trend_file: Optional[Path] = None, safety_factor: float = 1.5):
        """
        Initialize resource estimator.
        
        Args:
            trend_file: Path to trend parameters file (.pkl or .json)
            safety_factor: Safety factor to apply to estimates (default: 1.5)
        """
        self.safety_factor = safety_factor
        self.trends = None
        
        if trend_file is None:
            # Try default locations
            default_locations = [
                Path(__file__).parent.parent.parent.parent / "benchmarks" / "fasta_results" / "trend_parameters.pkl",
                Path("benchmarks/fasta_results/trend_parameters.pkl"),
                Path("fasta_results/trend_parameters.pkl"),
            ]
            
            for loc in default_locations:
                if loc.exists():
                    trend_file = loc
                    break
        
        if trend_file and trend_file.exists():
            self.load_trends(trend_file)
    
    def load_trends(self, trend_file: Path):
        """Load trend parameters from file."""
        if trend_file.suffix == '.pkl':
            with open(trend_file, 'rb') as f:
                self.trends = pickle.load(f)
        elif trend_file.suffix == '.json':
            with open(trend_file, 'r') as f:
                self.trends = json.load(f)
        else:
            raise ValueError(f"Unsupported file format: {trend_file.suffix}")
    
    def predict_from_trend(self, size: float, trend_params: Dict[str, float]) -> float:
        """Predict value from trend line parameters."""
        if not trend_params:
            return 0.0
        
        if trend_params.get('log_scale', False):
            import numpy as np
            log_size = np.log10(size)
            log_prediction = trend_params['slope'] * log_size + trend_params['intercept']
            return 10 ** log_prediction
        else:
            return trend_params['slope'] * size + trend_params['intercept']
    
    def estimate_resources(self, input_size: int, mode: str, 
                          max_threads: int = 1) -> Dict[str, Any]:
        """
        Estimate resources for a factorization job.
        
        Args:
            input_size: Input size in nucleotides
            mode: Factorization mode (with_reverse_complement or without_reverse_complement)
            max_threads: Maximum number of threads available
            
        Returns:
            Dictionary with resource estimates
        """
        if self.trends is None:
            # Return conservative defaults if no trends available
            return self._default_estimates(input_size, max_threads)
        
        # Determine function name based on mode
        if mode == "with_reverse_complement":
            func_name = "write_factors_binary_file_fasta_multiple_dna_w_rc"
        else:
            func_name = "write_factors_binary_file_fasta_multiple_dna_no_rc"
        
        # Get trend parameters
        time_key = f"{func_name}_time"
        memory_key = f"{func_name}_memory"
        disk_key = f"{func_name}_disk_space"
        
        estimate = {
            'input_size': input_size,
            'input_size_mbp': input_size / 1_000_000,
            'mode': mode,
            'threads': min(max_threads, 1),  # Single-threaded for now
        }
        
        # Estimate time (in seconds)
        if time_key in self.trends:
            time_ms = self.predict_from_trend(input_size, self.trends[time_key])
            base_time_sec = time_ms / 1000
            safe_time_sec = base_time_sec * self.safety_factor
            estimate['estimated_time_seconds'] = base_time_sec
            estimate['safe_time_seconds'] = safe_time_sec
            estimate['safe_time_minutes'] = safe_time_sec / 60
            estimate['safe_time_hours'] = safe_time_sec / 3600
            # Round up to next 5 minutes for LSF
            estimate['lsf_time_minutes'] = int((safe_time_sec / 60) + 4) // 5 * 5 + 5
        else:
            estimate['lsf_time_minutes'] = 60  # Default 1 hour
        
        # Estimate memory (in GB)
        if memory_key in self.trends:
            memory_mb = self.predict_from_trend(input_size, self.trends[memory_key])
            memory_gb = memory_mb / 1024
            safe_memory_gb = memory_gb * self.safety_factor
            estimate['estimated_memory_gb'] = memory_gb
            estimate['safe_memory_gb'] = safe_memory_gb
            # Round up to nearest GB
            estimate['lsf_memory_gb'] = max(1, int(safe_memory_gb + 0.99))
        else:
            estimate['lsf_memory_gb'] = 4  # Default 4 GB
        
        # Estimate disk space (in GB)
        if disk_key in self.trends:
            disk_mb = self.predict_from_trend(input_size, self.trends[disk_key])
            disk_gb = disk_mb / 1024
            safe_disk_gb = disk_gb * self.safety_factor
            estimate['estimated_disk_gb'] = disk_gb
            estimate['safe_disk_gb'] = safe_disk_gb
            # Round up to nearest GB
            estimate['lsf_disk_gb'] = max(1, int(safe_disk_gb + 0.99))
        else:
            estimate['lsf_disk_gb'] = 1  # Default 1 GB
        
        return estimate
    
    def _default_estimates(self, input_size: int, max_threads: int) -> Dict[str, Any]:
        """Return conservative default estimates when trends not available."""
        # Very rough estimates: 1 second per 10kb, 10 bytes per nucleotide
        size_mb = input_size / 1_000_000
        
        return {
            'input_size': input_size,
            'input_size_mbp': size_mb,
            'threads': min(max_threads, 1),
            'lsf_time_minutes': max(30, int(size_mb * 10)),
            'lsf_memory_gb': max(4, int(size_mb * 10)),
            'lsf_disk_gb': max(1, int(size_mb)),
        }


def get_fasta_size(fasta_path: Path, logger: Optional[logging.Logger] = None) -> int:
    """
    Estimate the size of a FASTA file in nucleotides.
    
    Args:
        fasta_path: Path to FASTA file
        logger: Logger instance
        
    Returns:
        Approximate number of nucleotides
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    try:
        total_size = 0
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line.startswith('>'):
                    # Count nucleotides (ignoring whitespace)
                    total_size += len(line.replace(' ', '').replace('\t', ''))
        
        logger.debug(f"FASTA file {fasta_path.name} has approximately {total_size:,} nucleotides")
        return total_size
        
    except Exception as e:
        logger.warning(f"Could not determine size of {fasta_path}: {e}")
        # Use file size as rough estimate (assume ~1 char per nucleotide)
        file_size = fasta_path.stat().st_size
        return file_size // 2  # Divide by 2 to account for headers


class LSFJobManager:
    """Manages LSF job submission and monitoring."""
    
    def __init__(self, queue: Optional[str] = None, log_dir: Optional[Path] = None,
                 logger: Optional[logging.Logger] = None):
        """
        Initialize LSF job manager.
        
        Args:
            queue: LSF queue name
            log_dir: Directory for job logs
            logger: Logger instance
        """
        self.queue = queue
        self.log_dir = log_dir or Path("lsf_logs")
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logger or logging.getLogger(__name__)
        self.submitted_jobs = {}  # job_id -> job_info
    
    def submit_job(self, job_name: str, command: str, 
                   time_minutes: int, memory_gb: int, 
                   output_file: Path, error_file: Optional[Path] = None,
                   depends_on: Optional[List[str]] = None) -> Optional[str]:
        """
        Submit a job to LSF.
        
        Args:
            job_name: Name for the job
            command: Command to execute
            time_minutes: Wall time limit in minutes
            memory_gb: Memory requirement in GB
            output_file: Path for stdout log
            error_file: Path for stderr log (uses output_file if None)
            depends_on: List of job IDs this job depends on
            
        Returns:
            Job ID if successful, None otherwise
        """
        if error_file is None:
            error_file = output_file
        
        # Build bsub command
        bsub_cmd = ["bsub"]
        
        # Job name
        bsub_cmd.extend(["-J", job_name])
        
        # Queue
        if self.queue:
            bsub_cmd.extend(["-q", self.queue])
        
        # Resources
        bsub_cmd.extend(["-W", str(time_minutes)])
        bsub_cmd.extend(["-M", f"{memory_gb}GB"])
        bsub_cmd.extend(["-R", f"rusage[mem={memory_gb}GB]"])
        
        # Output files
        bsub_cmd.extend(["-o", str(output_file)])
        bsub_cmd.extend(["-e", str(error_file)])
        
        # Dependencies
        if depends_on:
            dep_str = " && ".join([f"done({jid})" for jid in depends_on])
            bsub_cmd.extend(["-w", dep_str])
        
        # Command
        bsub_cmd.append(command)
        
        # Submit job
        try:
            result = subprocess.run(
                bsub_cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            # Parse job ID from output
            # LSF output format: "Job <12345> is submitted to queue <normal>."
            output = result.stdout.strip()
            if "Job <" in output:
                job_id = output.split("Job <")[1].split(">")[0]
                self.logger.info(f"Submitted job {job_id}: {job_name}")
                
                self.submitted_jobs[job_id] = {
                    'name': job_name,
                    'command': command,
                    'time_minutes': time_minutes,
                    'memory_gb': memory_gb,
                    'output_file': output_file,
                    'error_file': error_file,
                }
                
                return job_id
            else:
                self.logger.error(f"Could not parse job ID from bsub output: {output}")
                return None
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to submit job {job_name}: {e}")
            self.logger.error(f"bsub stderr: {e.stderr}")
            return None
        except Exception as e:
            self.logger.error(f"Unexpected error submitting job {job_name}: {e}")
            return None
    
    def check_job_status(self, job_id: str) -> Optional[str]:
        """
        Check the status of a job.
        
        Args:
            job_id: Job ID to check
            
        Returns:
            Job status (PEND, RUN, DONE, EXIT, etc.) or None if not found
        """
        try:
            result = subprocess.run(
                ["bjobs", "-noheader", job_id],
                capture_output=True,
                text=True,
                check=False
            )
            
            if result.returncode == 0 and result.stdout.strip():
                # Parse status from bjobs output
                # Format: JOBID USER STAT QUEUE FROM_HOST EXEC_HOST JOB_NAME SUBMIT_TIME
                parts = result.stdout.strip().split()
                if len(parts) >= 3:
                    return parts[2]  # STAT column
            
            return None
            
        except Exception as e:
            self.logger.warning(f"Error checking job status for {job_id}: {e}")
            return None
    
    def wait_for_jobs(self, job_ids: List[str], check_interval: int = 60) -> Dict[str, str]:
        """
        Wait for jobs to complete and return their final status.
        
        Args:
            job_ids: List of job IDs to monitor
            check_interval: Time between status checks in seconds
            
        Returns:
            Dictionary mapping job_id to final status
        """
        pending_jobs = set(job_ids)
        final_status = {}
        
        self.logger.info(f"Monitoring {len(job_ids)} jobs...")
        
        while pending_jobs:
            for job_id in list(pending_jobs):
                status = self.check_job_status(job_id)
                
                if status is None or status in ["DONE", "EXIT"]:
                    # Job completed or not found
                    final_status[job_id] = status or "UNKNOWN"
                    pending_jobs.remove(job_id)
                    
                    job_info = self.submitted_jobs.get(job_id, {})
                    job_name = job_info.get('name', job_id)
                    
                    if status == "DONE":
                        self.logger.info(f"Job {job_id} ({job_name}) completed successfully")
                    elif status == "EXIT":
                        self.logger.warning(f"Job {job_id} ({job_name}) failed")
                    else:
                        self.logger.warning(f"Job {job_id} ({job_name}) status unknown")
            
            if pending_jobs:
                self.logger.info(f"{len(pending_jobs)} jobs still running...")
                time.sleep(check_interval)
        
        return final_status
    
    def generate_summary(self, job_status: Dict[str, str], output_file: Optional[Path] = None):
        """
        Generate a summary of job results.
        
        Args:
            job_status: Dictionary mapping job_id to status
            output_file: Optional file to write summary to
        """
        summary_lines = []
        summary_lines.append("=" * 80)
        summary_lines.append("LSF BATCH FACTORIZATION SUMMARY")
        summary_lines.append("=" * 80)
        
        total_jobs = len(job_status)
        successful = sum(1 for s in job_status.values() if s == "DONE")
        failed = sum(1 for s in job_status.values() if s == "EXIT")
        unknown = sum(1 for s in job_status.values() if s not in ["DONE", "EXIT"])
        
        summary_lines.append(f"Total jobs: {total_jobs}")
        summary_lines.append(f"Successful: {successful}")
        summary_lines.append(f"Failed: {failed}")
        summary_lines.append(f"Unknown: {unknown}")
        
        if failed > 0 or unknown > 0:
            summary_lines.append("\nFailed/Unknown Jobs:")
            summary_lines.append("-" * 80)
            
            for job_id, status in job_status.items():
                if status != "DONE":
                    job_info = self.submitted_jobs.get(job_id, {})
                    job_name = job_info.get('name', job_id)
                    output_file_path = job_info.get('output_file', 'N/A')
                    
                    summary_lines.append(f"Job ID: {job_id}")
                    summary_lines.append(f"  Name: {job_name}")
                    summary_lines.append(f"  Status: {status}")
                    summary_lines.append(f"  Log: {output_file_path}")
                    summary_lines.append("")
        
        summary_text = "\n".join(summary_lines)
        
        # Print to logger
        self.logger.info("\n" + summary_text)
        
        # Write to file if requested
        if output_file:
            with open(output_file, 'w') as f:
                f.write(summary_text)
            self.logger.info(f"Summary written to {output_file}")


def prepare_jobs(file_list: List[str], output_dir: Path, mode: str,
                estimator: ResourceEstimator, max_threads: int = 1,
                logger: Optional[logging.Logger] = None) -> List[Dict[str, Any]]:
    """
    Prepare job specifications for all files.
    
    Args:
        file_list: List of FASTA file paths
        output_dir: Output directory
        mode: Factorization mode
        estimator: Resource estimator
        max_threads: Maximum threads per job
        logger: Logger instance
        
    Returns:
        List of job specifications
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    jobs = []
    
    for file_path in file_list:
        fasta_path = Path(file_path)
        
        if not fasta_path.exists():
            logger.warning(f"File not found, skipping: {file_path}")
            continue
        
        # Get file size
        input_size = get_fasta_size(fasta_path, logger)
        
        # Estimate resources
        resources = estimator.estimate_resources(input_size, mode, max_threads)
        
        # Determine output paths
        output_paths = batch_factorize.get_output_paths(fasta_path, output_dir, mode)
        
        # Create job specification
        job = {
            'input_file': fasta_path,
            'output_paths': output_paths,
            'mode': mode,
            'resources': resources,
        }
        
        jobs.append(job)
        
        logger.debug(f"Job for {fasta_path.name}: "
                    f"{resources['lsf_time_minutes']}min, "
                    f"{resources['lsf_memory_gb']}GB memory")
    
    return jobs


def main():
    """Main entry point for LSF batch factorization."""
    parser = argparse.ArgumentParser(
        description="Submit FASTA factorization jobs to LSF cluster with intelligent resource allocation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Submit jobs with default settings
  python -m noLZSS.genomics.lsf_batch_factorize --file-list files.txt --output-dir results
  
  # Specify LSF queue and max threads
  python -m noLZSS.genomics.lsf_batch_factorize --file-list files.txt --output-dir results --queue normal --max-threads 4
  
  # Use custom trend parameters file
  python -m noLZSS.genomics.lsf_batch_factorize --file-list files.txt --output-dir results --trend-file my_trends.pkl
  
  # Submit without waiting for completion
  python -m noLZSS.genomics.lsf_batch_factorize --file-list files.txt --output-dir results --no-wait
        """
    )
    
    # Input specification
    parser.add_argument(
        "--file-list", type=Path, required=True,
        help="Text file containing list of FASTA file paths (one per line)"
    )
    parser.add_argument(
        "files", nargs="*",
        help="Additional FASTA file paths to process"
    )
    
    # Output configuration
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory for binary factorization results"
    )
    parser.add_argument(
        "--mode", 
        choices=[batch_factorize.FactorizationMode.WITHOUT_REVERSE_COMPLEMENT, 
                batch_factorize.FactorizationMode.WITH_REVERSE_COMPLEMENT],
        default=batch_factorize.FactorizationMode.WITH_REVERSE_COMPLEMENT,
        help="Factorization mode (default: with_reverse_complement)"
    )
    
    # LSF configuration
    parser.add_argument(
        "--queue", type=str,
        help="LSF queue name (uses default if not specified)"
    )
    parser.add_argument(
        "--max-threads", type=int, default=1,
        help="Maximum number of threads per job (default: 1)"
    )
    parser.add_argument(
        "--job-name-prefix", type=str, default="nolzss",
        help="Prefix for LSF job names (default: nolzss)"
    )
    parser.add_argument(
        "--log-dir", type=Path, default=Path("lsf_logs"),
        help="Directory for LSF job logs (default: lsf_logs)"
    )
    
    # Resource estimation
    parser.add_argument(
        "--trend-file", type=Path,
        help="Path to trend parameters file (.pkl or .json). Uses default location if not specified."
    )
    parser.add_argument(
        "--safety-factor", type=float, default=1.5,
        help="Safety factor for resource estimates (default: 1.5)"
    )
    
    # Execution options
    parser.add_argument(
        "--no-wait", action="store_true",
        help="Submit jobs and exit without waiting for completion"
    )
    parser.add_argument(
        "--check-interval", type=int, default=60,
        help="Interval between job status checks in seconds (default: 60)"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print job specifications without submitting"
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
    logger = batch_factorize.setup_logging(args.log_level, args.log_file)
    
    try:
        # Get file list
        file_list = batch_factorize.read_file_list(args.file_list, logger)
        
        if args.files:
            file_list.extend(args.files)
        
        logger.info(f"Processing {len(file_list)} FASTA files")
        logger.info(f"Mode: {args.mode}")
        logger.info(f"Output directory: {args.output_dir}")
        
        # Initialize resource estimator
        logger.info("Loading resource estimation parameters...")
        estimator = ResourceEstimator(args.trend_file, args.safety_factor)
        
        if estimator.trends is None:
            logger.warning("No trend parameters found. Using conservative defaults.")
        else:
            logger.info("Loaded benchmark trend parameters successfully")
        
        # Prepare jobs
        logger.info("Preparing job specifications...")
        jobs = prepare_jobs(
            file_list, args.output_dir, args.mode,
            estimator, args.max_threads, logger
        )
        
        if not jobs:
            logger.error("No valid jobs to submit")
            sys.exit(1)
        
        # Print job summary
        total_time_hours = sum(j['resources'].get('safe_time_hours', 0) for j in jobs)
        max_memory_gb = max(j['resources'].get('lsf_memory_gb', 0) for j in jobs)
        total_disk_gb = sum(j['resources'].get('lsf_disk_gb', 0) for j in jobs)
        
        logger.info(f"\nJob Summary:")
        logger.info(f"  Total jobs: {len(jobs)}")
        logger.info(f"  Estimated total time: {total_time_hours:.1f} hours")
        logger.info(f"  Max memory per job: {max_memory_gb} GB")
        logger.info(f"  Total disk space: {total_disk_gb} GB")
        
        if args.dry_run:
            logger.info("\n=== DRY RUN MODE - Job Specifications ===")
            for i, job in enumerate(jobs, 1):
                res = job['resources']
                logger.info(f"\nJob {i}: {job['input_file'].name}")
                logger.info(f"  Input size: {res.get('input_size_mbp', 0):.2f} Mbp")
                logger.info(f"  Time: {res.get('lsf_time_minutes', 0)} minutes")
                logger.info(f"  Memory: {res.get('lsf_memory_gb', 0)} GB")
                logger.info(f"  Disk: {res.get('lsf_disk_gb', 0)} GB")
            
            logger.info("\nDry run completed. No jobs submitted.")
            sys.exit(0)
        
        # Initialize job manager
        job_manager = LSFJobManager(args.queue, args.log_dir, logger)
        
        # Submit jobs
        logger.info("\nSubmitting jobs to LSF...")
        submitted_job_ids = []
        
        for i, job in enumerate(jobs, 1):
            input_file = job['input_file']
            output_paths = job['output_paths']
            resources = job['resources']
            
            # Create job name
            job_name = f"{args.job_name_prefix}_{i}_{input_file.stem}"
            
            # Get the output path for this mode
            if args.mode == batch_factorize.FactorizationMode.WITH_REVERSE_COMPLEMENT:
                output_file = output_paths.get("with_reverse_complement")
            else:
                output_file = output_paths.get("without_reverse_complement")
            
            # Create command using the worker script
            mode_arg = "with_rc" if args.mode == batch_factorize.FactorizationMode.WITH_REVERSE_COMPLEMENT else "no_rc"
            command = f"python -m noLZSS.genomics.lsf_job_worker '{input_file}' '{output_file}' --mode {mode_arg}"
            
            # Log files
            log_file = args.log_dir / f"{job_name}.log"
            
            # Submit job
            job_id = job_manager.submit_job(
                job_name=job_name,
                command=command,
                time_minutes=resources['lsf_time_minutes'],
                memory_gb=resources['lsf_memory_gb'],
                output_file=log_file,
            )
            
            if job_id:
                submitted_job_ids.append(job_id)
            else:
                logger.error(f"Failed to submit job for {input_file.name}")
        
        logger.info(f"\nSubmitted {len(submitted_job_ids)} jobs successfully")
        
        if not submitted_job_ids:
            logger.error("No jobs were submitted successfully")
            sys.exit(1)
        
        # Save job tracking info
        tracking_file = args.log_dir / "job_tracking.json"
        with open(tracking_file, 'w') as f:
            json.dump({
                'submitted_jobs': job_manager.submitted_jobs,
                'job_ids': submitted_job_ids,
            }, f, indent=2, default=str)
        
        logger.info(f"Job tracking information saved to {tracking_file}")
        
        # Wait for completion or exit
        if args.no_wait:
            logger.info("\nJobs submitted. Exiting without waiting for completion.")
            logger.info(f"Use 'bjobs' to monitor jobs or check logs in {args.log_dir}")
            sys.exit(0)
        
        # Monitor jobs
        logger.info("\nWaiting for jobs to complete...")
        job_status = job_manager.wait_for_jobs(submitted_job_ids, args.check_interval)
        
        # Generate summary
        summary_file = args.log_dir / "job_summary.txt"
        job_manager.generate_summary(job_status, summary_file)
        
        # Exit with appropriate code
        failed_count = sum(1 for s in job_status.values() if s != "DONE")
        if failed_count > 0:
            logger.warning(f"Completed with {failed_count} failures")
            sys.exit(1)
        else:
            logger.info("All jobs completed successfully")
            sys.exit(0)
    
    except LSFBatchError as e:
        logger.error(f"LSF batch error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        logger.info("Interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
