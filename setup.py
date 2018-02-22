from setuptools import setup, find_packages

setup(
    name="Fragment Hotspot Maps",
    version="0.1",
    packages=find_packages('src'),
    package_dir={'': 'src'},
    include_package_data=True,
    scripts=['src/run_hotspot.py'],

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires=['numpy', 'csd-python-api', 'ccdc_internal', 'matplotlib', 'scipy', 'sklearn'],

    package_data={
        # If any package contains *.txt or *.rst files, include them:
        '': ['*.txt', '*.rst','*.pkl'],
        # And include any *.mol2 files found in the 'probes' package, too:
        'fhm': ['probes/*.mol2'],
    }
)
