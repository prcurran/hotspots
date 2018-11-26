from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

print find_packages()

setup(
    name="hotspots",
    authors="Chris Radoux, Peter Curran",
    version="0.1",
    packages=find_packages(),
    include_package_data=True,
    # scripts=['src/run_hotspot.py'],
    long_description=long_description,
    long_description_content_type="text/markdown",
    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires=['numpy>=1.9', 'csd-python-api>=2.0', 'matplotlib', 'scipy', 'sklearn', 'scikit-image', 'pandas',
                      'futures', 'cython', 'tqdm', 'nglview'],
    package_data={
        # If any package contains *.txt or *.rst files, include them:
        '': ['*.txt', '*.rst', '*.pkl'],
        # And include any *.mol2 files found in the 'probes' package, too:
        'hotspots': ['probes/*.mol2'],
    }
)
