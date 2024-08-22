import subprocess
from importlib.metadata import version as get_version
from importlib.metadata import PackageNotFoundError

def get_version_info():
    package_name = "satute"  # Your package name
    default_version = "1.0.0"  # Default fallback version if all else fails

    try:
        # First, try to get the version from installed package metadata
        return get_version(package_name)
    except PackageNotFoundError:
        # If the package is not found, perhaps it's being run directly from the source
        try:
            # Try to fetch the current Git commit hash
            commit_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'], stderr=subprocess.STDOUT).decode().strip()
            return f"dev-{commit_hash[:7]}"  # Return a development version based on commit hash
        except subprocess.CalledProcessError:
            # If Git is not available or fails, return the default version
            return default_version


SATUTE_VERSION = f"""\n\nSatuTe Version {get_version_info()} \n\nTo cite SatuTe please use:\nC. Manuel, C. Elgert, E. Sakalli, H. A. Schmidt, C. Viñas and A. von Haeseler\nWhen the past fades: Detecting phylogenetic signal with SatuTe\nNature Methods\nDOI: \n\n"""