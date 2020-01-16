from __future__ import print_function
import os
from concurrent import futures


def run(cmd):
    print(cmd)
    os.system(cmd)

parent = "/home/pcurran/github_packages/hotspots/examples/2_docking/virtual_screening/"

data = os.path.join(parent, 'akt1')
h = os.path.join(parent,'hotspot')
script = "virtual_screening.py"
weights = [0, 10, 100]

cmds = ["""python virtual_screening.py {} {} {}""".format(script, data, h, weight) for weight in weights]

with futures.ProcessPoolExecutor(max_workers=3) as executor:
    executor.map(run, cmds)
