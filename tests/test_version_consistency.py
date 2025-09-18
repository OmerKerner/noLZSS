"""
Test version consistency across the noLZSS package.

This test ensures that all version sources report the same version.
"""


def test_version_consistency():
    """Test that all version sources agree."""
    import noLZSS
    from importlib.metadata import version as pkg_version, PackageNotFoundError
    
    try:
        package_version = pkg_version("noLZSS")
        reported_version = noLZSS.__version__
        
        print(f"Package metadata version: {package_version}")
        print(f"noLZSS.__version__: {reported_version}")
        
        assert package_version == reported_version, f"Version mismatch: package={package_version}, module={reported_version}"
        print("✓ All versions are consistent")
    except PackageNotFoundError:
        print("Package not installed via pip, testing C++ extension version directly...")
        try:
            from noLZSS._noLZSS import __version__ as cpp_version
            assert noLZSS.__version__ == cpp_version, f"Version mismatch: python={noLZSS.__version__}, cpp={cpp_version}"
            print(f"✓ C++ extension version matches: {cpp_version}")
        except ImportError:
            print("✓ C++ extension not available, fallback version used")


def test_version_fallback():
    """Test that version fallback works when C++ extension is missing."""
    import sys
    import importlib
    
    # Temporarily hide the C++ extension to test fallback
    original_modules = sys.modules.copy()
    
    try:
        # Remove the C++ extension if it's loaded
        if 'noLZSS._noLZSS' in sys.modules:
            del sys.modules['noLZSS._noLZSS']
        if 'noLZSS' in sys.modules:
            del sys.modules['noLZSS']
        
        # Mock the C++ import to fail
        class MockImportError:
            def __init__(self):
                pass
            def __getattr__(self, name):
                raise ImportError("Mock C++ extension not available")
        
        # Test that the fallback mechanism works
        sys.modules['noLZSS._noLZSS'] = MockImportError()
        
        # Re-import noLZSS, should use fallback
        import noLZSS
        importlib.reload(noLZSS)
        
        # Check that version is still available via fallback
        assert hasattr(noLZSS, '__version__')
        
        # Determine expected fallback version based on package installation
        from importlib.metadata import version as pkg_version, PackageNotFoundError
        try:
            expected_version = pkg_version("noLZSS")
            print(f"Fallback version: {noLZSS.__version__} (from package metadata)")
            assert noLZSS.__version__ == expected_version, f"Expected fallback to package metadata version {expected_version}, got {noLZSS.__version__}"
        except PackageNotFoundError:
            # Package not installed, should fallback to "0.0.0"
            print(f"Fallback version: {noLZSS.__version__} (package not installed)")
            assert noLZSS.__version__ == "0.0.0", f"Expected fallback version 0.0.0, got {noLZSS.__version__}"
        
        print("✓ Version fallback mechanism works")
        
    finally:
        # Restore original modules
        sys.modules.clear()
        sys.modules.update(original_modules)


if __name__ == "__main__":
    print("Running version consistency tests...")
    
    try:
        test_version_consistency()
        test_version_fallback()
        print("\n✓ All version tests passed!")
    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        raise