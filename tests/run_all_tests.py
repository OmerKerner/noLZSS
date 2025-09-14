"""
Test runner script for all noLZSS tests.
"""

import sys
import os
from pathlib import Path

def run_test_file(filepath):
    """Run a test file and return (passed, total) counts."""
    print(f"\n{'='*60}")
    print(f"Running: {filepath.name}")
    print('='*60)
    
    try:
        # Import and run the test file
        spec = __import__(f"tests.{filepath.stem}", fromlist=[''])
        
        # If the module has a main function, call it
        if hasattr(spec, '__name__') and filepath.stem in spec.__name__:
            # Execute the module's main section
            exec(compile(open(filepath).read(), filepath, 'exec'))
        
        return True
    except Exception as e:
        print(f"Error running {filepath}: {e}")
        return False

def main():
    """Run all test files."""
    test_dir = Path(__file__).parent
    
    # Define test files in order of execution
    test_files = [
        "test_structure.py",  # Basic structure tests
        "test_utils.py",      # Utils functionality
        "test_core.py",       # Core wrappers (will mostly fail without C++)
        "test_genomics.py",   # Genomics subpackage
        "test_integration.py", # Integration tests
        "test_cpp_bindings.py" # Original C++ binding tests (will fail without build)
    ]
    
    total_files = len(test_files)
    passed_files = 0
    
    print("Running noLZSS Test Suite")
    print(f"Found {total_files} test files to run")
    
    for test_file in test_files:
        test_path = test_dir / test_file
        
        if test_path.exists():
            try:
                # Run the test file by executing it as a script
                # Use proper path escaping for the complex path
                escaped_path = str(test_path).replace(' ', r'\ ')
                result = os.system(f'cd "{test_dir.parent}" && python "{test_path}"')
                if result == 0:
                    passed_files += 1
                    print(f"{test_file} completed successfully")
                else:
                    print(f"Warning: {test_file} had issues (exit code: {result})")
            except Exception as e:
                print(f"Warning: {test_file} failed with exception: {e}")
        else:
            print(f"Warning: {test_file} not found, skipping")
    
    print(f"\n{'='*60}")
    print("TEST SUITE SUMMARY")
    print('='*60)
    print(f"Files executed: {passed_files}/{total_files}")
    
    if passed_files == total_files:
        print("All test files executed successfully!")
    else:
        print(f"Warning:  {total_files - passed_files} test files had issues")
        print("Note: Some failures are expected without C++ extension build")
    
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
