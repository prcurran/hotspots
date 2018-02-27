#!/usr/bin/env python

"package fragment_hotspot_maps 0.1"

from __future__ import division, absolute_import, print_function

import glob
import itertools
import os
import re
import sys
import zipfile
try:
    import zlib
    compression = zipfile.ZIP_DEFLATED
except:
    compression = zipfile.ZIP_STORED

modes = {zipfile.ZIP_DEFLATED: 'deflated', zipfile.ZIP_STORED: 'stored'}

# Parse, don't import to get the version because importing may fail as its
# dependencies may not be in place at the time this is run.
with open(os.path.join('fragment_hotspot_maps', 'src', '__init__.py')) as init_py:
    version = re.search(r'''\b__version__\s*=\s*['"]([^'"]+)''', init_py.read()).group(1)


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)  
    parser.add_argument('--outdir', default='.')
    args = parser.parse_args(argv)

    zf_name = os.path.join(args.outdir, 'fragment-hotspot-maps' + version + '.zip')
    zf = zipfile.ZipFile(zf_name, mode='w')
    try:
        for example_file in itertools.chain(
                glob.glob(os.path.join('fragment_hotspot_maps', '*')),
                glob.glob(os.path.join('fragment_hotspot_maps', 'src', '*')),
                glob.glob(os.path.join('fragment_hotspot_maps', 'src' 'probes'))):
            zf.write(example_file, compress_type=compression)
    finally:
        print('Done')
        zf.close()

    return 0

if __name__ == '__main__':
    sys.exit(main())
