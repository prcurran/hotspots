
from __future__ import print_function
from setuptools import setup, find_packages

print(find_packages())

setup(
    name="hotspots",
    author="Chris Radoux, Peter Curran, Mihaela Smilova",
    author_email="pcurran@ccdc.cam.ac.uk",
	license="MIT",
    version="1.0.0-rc-1",
    url="https://github.com/prcurran/hotspots",
    packages=find_packages(),
    include_package_data=True,
    # scripts=['src/run_hotspot.py'],
    long_description=read_file('README.md'),
    long_description_content_type='text/markdown',
    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires=['numpy==1.15.4',
                      'csd-python-api>=2.0.0',
                      'rdkit==2018.09.1',
                      'matplotlib==2.2.3',
                      'scipy',
                      'sklearn',
                      'scikit-image==0.14.2',
                      'pandas',
                      'futures',
                      'cython==0.29.5',
                      'tqdm==4.31.1',
                      'xmltodict==0.12.0'],
    package_data={
        # If any package contains *.txt or *.rst files, include them:
        '': ['*.txt', '*.rst', '*.pkl'],
        # And include any *.mol2 files found in the 'probes' package, too:
        'hotspots': ['probes/*.mol2'],
    }
)
