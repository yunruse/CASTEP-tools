'''
Analysis script by Mia Dobson.
Outputs graph of desired statistic.
'''

from argparse import ArgumentParser
from itertools import product
from os.path import isfile

from matplotlib import pyplot
from numpy import array, arange, pi
from numpy.linalg import pinv

from parse_cell import CellFile
from parse_md import MDFile

D = [-1, 0, 1]
DELTAS = array(list(product(D, D, D)))

TIMELABEL = 'Time / ps'

# TODO: some state object so you don't need
# to keep asking for cell and cellinv dammit!

# TODO: multiple graphs! take in `out path/to/yee-$type.png` maybe?
# use a wrapper on each func I guess
    
def cell_square_dist(a, b, cell, cellinv):
    '''
    Return |b - a|**2 considering that
    a and b may have "wrapped around" a cell
    (i.e. pick the smallest from all combinations).
    '''
    if cell is None:
        return ((b - a)**2).sum()
    # form all 3x3x3 images of b and find
    # minimum distance to a
    b_imgs = ((b @ cellinv) + DELTAS) @ cell
    img_dists = ((b_imgs - a)**2).sum(1)
    return min(img_dists)

OPT_TEMP = 'temp'
def temperature(md):
    pyplot.title(f'Temperature of `{md.path}`')
    pyplot.xlabel(TIMELABEL)
    pyplot.ylabel('Temperature / K')
    pyplot.plot(
        [step.t / 1000 for step in md.steps],
        [step.T for step in md.steps]
    )

OPT_PRESSURE = 'pressure'
def pressure(md):
    pyplot.title(f'Pressure of `{md.path}`')
    pyplot.xlabel(TIMELABEL)
    pyplot.ylabel('Pressure / GPa')
    pyplot.plot(
        [step.t / 1000 for step in md.steps],
        [step.P for step in md.steps]
    )

def msd(step, ref, cell, cellinv):
    '''Mean square displacement from reference cell'''        
    species = {}
    for a, b in zip(step.ions, ref.ions):
        assert a.species == b.species and a.number == b.number
        if a.species not in species:
            species[a.species] = []
        dist = cell_square_dist(a.pos, b.pos, cell, cellinv)
        species[a.species].append(dist)
    return {k: sum(v)/len(v) for k, v in species.items()}

OPT_MSD = 'msd'
def msds(md, cell, cellinv):
    pyplot.title(f'Mean square displacement of `{md.path}`')
    
    pyplot.xlabel(TIMELABEL)
    t = [step.t / 1000 for step in md.steps[1:]]
    
    pyplot.ylabel('Mean square displacement / Å²')
    msds = [
        msd(step, md.steps[0], cell, cellinv)
        for step in md.steps[1:]]
    for k in msds[0]:
        pyplot.plot(t, [i[k] for i in msds], label=k)

    pyplot.legend()


def bonds(step, cell, cellinv):
    C = [i.pos for i in step.ions if i.species == 'C']
    H = [i.pos for i in step.ions if i.species == 'H']
    CH = array([
        cell_square_dist(c, h, cell, cellinv) ** 0.5
        for c in C for h in H
    ])
    HH = array([
        cell_square_dist(H[a], H[b], cell, cellinv) ** 0.5
        for a in range(len(H)) for b in range(a)
    ])

    return CH, HH

def bonds_gen(md, cell, cellinv):
    CH, HH = bonds(md.steps[0], cell, cellinv)
    N = 0
    for i in arange(dstep, len(md.steps), dstep):
        N += 1
        ch, hh = bonds(md.steps[i], cell, cellinv)
        CH += ch
        HH += hh
        
    CH = [i / N for i in CH]
    HH = [i / N for i in HH]
    return CH, HH

OPT_RDF = 'rdf'
def rdf(md, cell, cellinv):
    BOND_CH = 1.08595
    BOND_HH = 0.7414
    BOND_CH = BOND_HH = 1
    DR = 0.08
    
    pyplot.title(f'Radial distribution function of `{md.path}`')
    pyplot.xlabel('Bond length $r$')
    pyplot.ylabel('Radial distribution $g(r)$')

    CH, HH = bonds(md.steps[-1], cell, cellinv)
    #CH, HH = bonds_gen(md, cell, cellinv)

    CH = list(CH)
    HH = list(HH)
    CH.sort(reverse=True)
    HH.sort(reverse=True)
    
    maxbond = max(HH[0], CH[0])
    lengths = arange(DR, maxbond, DR)
    
    ch_hist = []
    hh_hist = []
    for r in lengths:
        n_ch = 0
        while CH and CH[-1] <= r:
            CH.pop()
            n_ch += 1
        ch_hist.append(n_ch / (4*pi*r**2) / DR)
        
        n_hh = 0
        while HH and HH[-1] <= r:
            HH.pop()
            n_hh += 1
        hh_hist.append(n_hh / (4*pi*r**2) / DR)
    
    pyplot.plot(lengths / BOND_CH, ch_hist, label='CH')
    pyplot.plot(lengths / BOND_HH, hh_hist, label='HH')
    pyplot.plot([0, maxbond], [1, 1], color='gray', linestyle='-')
    pyplot.plot([1, 1], [0, 150], color='gray', linestyle='-')
    pyplot.legend()
    

OPTS_CELLS = [OPT_MSD, OPT_RDF]
OPTS = OPTS_CELLS + [OPT_TEMP, OPT_PRESSURE]

parser = ArgumentParser(description=__doc__)
parser.add_argument('path', type=str, help='the file path of the .md file')
parser.add_argument('out',  type=str, help='the file to output the graph')
parser.add_argument('option', choices=OPTS)
parser.add_argument('-c', '--cell', metavar='path', type=str, default=None,
                    help='''
manually specify .cell file path
(defaults to same name as .md file)
''')
parser.add_argument('--cellstat', metavar='path', type=str, default=None,
                    help='''
path for hydrogen-tagging file for .cell files
''')
parser.add_argument('-t', '--terminal', action='store_true',
                    help='configure graph for better terminal display')

def find(path):
    if not isfile(path):
        raise FileNotFoundError(f'Could not find `{path}`')

def main(args):
    args = parser.parse_args(args)
    
    if args.terminal:
        pyplot.rcParams.update({
            'font.size': 12,
            'figure.figsize': (4.4, 2.2),
            'figure.dpi': 50
        })

    find(args.path)
    cell = None
    cellinv = None
    if args.option in OPTS_CELLS:
        cellpath = args.cell
        if cellpath is None:
            cellpath = args.path.replace('.md', '.cell')
        find(cellpath)
        cell = CellFile(cellpath).cell_vectors
        cellinv = pinv(cell)
    
    md = MDFile(args.path, args.cellstat)
    
    if args.option == OPT_MSD:
        msds(md, cell, cellinv)
    elif args.option == OPT_RDF:
        rdf(md, cell, cellinv)
    elif args.option == OPT_TEMP:
        temperature(md)
    elif args.option == OPT_PRESSURE:
        pressure(md)
    else:
        print(f'Unknown option {args.option}')
    pyplot.savefig(args.out)
    
    return md, cell, cellinv

if __name__ == '__main__':
    from os import chdir
    chdir('/Users/yunruse/Downloads/')
    argv = [
        "castep/CH5_naumova_P_1-20gpa-300K.md",
        "rdf-ye.png",
        "rdf"
    ]
    md, cell, cellinv = main(argv)
