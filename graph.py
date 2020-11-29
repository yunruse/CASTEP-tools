'''
Display a Matplotlib Axes file as dumped.
'''

from pickle import load
from argparse import ArgumentParser

import matplotlib

# the default for mac, 'MacOSX', is glitchy
matplotlib.use('TkAgg')

from matplotlib import pyplot

parser = ArgumentParser(description=__doc__)
parser.add_argument('path', type=str, help='''
    the file path of the .md file OR a .analysis file
''')
argv = ['/Users/yunruse/Documents/University Work/SH Project/lab/graphs/CH5_naumova_P_1-20gpa-300K-msd.ax']
path = parser.parse_args(argv).path

# allow running from applescript
if ':' in path:
    _, path = path.replace(':', '/').split('/', 1)

with open(path, 'rb') as f:
    ax = load(f)

pyplot.show()
