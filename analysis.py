'''
Analysis script by Mia Dobson.
Outputs graph(s) of desired statistics.
'''

from argparse import ArgumentParser
from functools import wraps
from itertools import product
from os.path import isfile, isdir, split, splitext

from matplotlib import pyplot
from numpy import array, arange, pi, nan, dot, cross
from numpy.linalg import pinv, norm

from parse_cell import CellFile
from parse_md import MDFile

def find(path):
    if not isfile(path):
        raise FileNotFoundError(f'Could not find `{path}`')

FUNCS = {}
def method(command):
    # small wrapper for methods, binding them to the terminal
    # and adding niceties
    def wrapper(func):
        @wraps(func)
        def newfunc(self, *args, **kwargs):
            # ensure pyplot is clear
            pyplot.clf()
            result = func(self, *args, **kwargs)
            if result is NotImplemented:
                return
            # save with helpful formatting
            path = self.outpath.replace(
                '$name', self.name).replace(
                '$method', command)
            print(path)
            pyplot.savefig(path)
        FUNCS[command] = func.__name__
        return newfunc
    return wrapper

D = [-1, 0, 1]
DELTAS = array(list(product(D, D, D)))

TIMELABEL = 'Time / ps'

class Analysis:
    __slots__ = ['name', 'outpath', 'record_every', 'steps', 'cell', 'cellinv']
    
    def __init__(
        self,
        mdpath: str,
        outpath: str,
        cellpath=None,
        hydro_path=None,
        record_every: int = 1,
    ):
        find(mdpath)
        _, name = split(mdpath)
        self.name, _ = splitext(name)
        self.outpath = outpath
        
        self.record_every = record_every
        
        md = MDFile(
            mdpath,
            hydro_path,
            step_is_valid=lambda i: i % record_every == 0,
        )
        self.steps = [s for s in md.steps if s.ions]
        
        if cellpath is None:
            cellpath = mdpath.replace('.md', '.cell')

        self.has_hydro_tags = hydro_path is not None
        
        find(cellpath)
        self.cell = CellFile(cellpath).cell_vectors
        self.cellinv = pinv(self.cell)

    def plot_variable(self, name, unit, func):
        pyplot.title(f'{name} of `{self.name}`')
        pyplot.xlabel(TIMELABEL)
        pyplot.ylabel(f'{name} / {unit}')
        
        time = [step.t / 1000 for step in self.steps]
        var = [func(step) for step in self.steps]
        valid = [i for i in var if i is not nan]
        pyplot.plot(time, var)
        pyplot.plot([0, max(time)], [valid[0], valid[0]], color='gray', linestyle='-')

    @method('temp')
    def temperature(self):
        self.plot_variable('Temperature', 'K', lambda step: step.T)

    @method('pressure')
    def pressure(self):
        self.plot_variable('Pressure', 'as raw', lambda step: step.P)
    
    def cell_square_dist(self, a, b):
        '''
        Return |b - a|**2 considering that
        a and b may have "wrapped around" a cell
        (i.e. pick the smallest from all combinations).
        '''
        # form all 3x3x3 images of b and find
        # minimum distance to a
        b_imgs = ((b @ self.cellinv) + DELTAS) @ self.cell
        img_dists = ((b_imgs - a)**2).sum(1)
        return min(img_dists)
    
    def msd(self, step):
        '''Mean square displacement from reference cell'''
        species = {}
        for a, b in zip(step.ions, self.steps[0].ions):
            if not (a.species == b.species and a.number == b.number):
                print(f'unmatched: {a.species} {a.number}'
                      f' / {b.species} {b.number}')
            if a.species not in species:
                species[a.species] = []
            dist = self.cell_square_dist(a.pos, b.pos)
            species[a.species].append(dist)
        
        return {k: sum(v)/len(v) for k, v in species.items()}

    @method('msd')
    def msds(self):
        pyplot.title(f'Mean square displacement of `{self.name}`')
        steps = self.steps[1:]
        
        pyplot.xlabel(TIMELABEL)
        t = [step.t / 1000 for step in steps]
        
        pyplot.ylabel('Mean square displacement / Å²')
        msds = [self.msd(step) for step in steps]
        for k in msds[0]:
            pyplot.plot(t, [i[k] for i in msds], label=k)

        pyplot.legend()

    def bonds(self, step):
        C = [i.pos for i in step.ions if i.species == 'C']
        H_H2 = [i.pos for i in step.ions if i.species == 'H(H2)']
        H_CH4 = [i.pos for i in step.ions if i.species == 'H(CH4)']
        CH = array([
            self.cell_square_dist(c, h) ** 0.5
            for c in C for h in H_CH4
        ])
        HH = array([
            self.cell_square_dist(H_H2[a], H_H2[b]) ** 0.5
            for a in range(len(H_H2)) for b in range(a)
        ])

        return CH, HH

    def bonds_gen(self, dstep=500):
        steps = self.steps
        CH, HH = self.bonds(steps[0])
        N = 0
        for i in arange(dstep, len(steps), dstep):
            N += 1
            ch, hh = self.bonds(steps[i])
            CH += ch
            HH += hh
            
        CH = [i / N for i in CH]
        HH = [i / N for i in HH]
        return CH, HH

    @method('rdf')
    def rdf(self, dstep=250):
        BOND_CH = 1.08595
        BOND_HH = 0.7414
        BOND_CH = BOND_HH = 1
        DR = 0.04

        if not self.has_hydro_tags:
            print('Error! rdf requires hydrogen tagging')
            return NotImplemented
        
        pyplot.title(f'Radial distribution function of `{self.name}`')
        pyplot.xlabel(f'Bond length $r$ (dstep={dstep})')
        pyplot.ylabel('Radial distribution $g(r)$')

        #CH, HH = self.bonds(self.steps[-1])
        CH, HH = self.bonds_gen(dstep)
        
        maxbond = max(HH[0], CH[0])
        lengths = arange(DR, maxbond, DR)

        def hist(bonds, rho):
            bonds = sorted(bonds, reverse=True)
            hist = []
            for r in lengths:
                n = 0
                while bonds and bonds[-1] <= r:
                    bonds.pop()
                    n += 1
                hist.append(rho * n / (4*pi*r**2) / DR)
            return hist

        ions = self.steps[0].ions
        
        n_H2 = len([i.pos for i in ions if i.species == 'H(H2)'])
        n_CH4 = len([i.pos for i in ions if i.species == 'H(CH4)'])
        
        a, b, c = self.cell
        V = norm(dot(cross(a, b), c))
        ch_hist = hist(CH, n_CH4 / V)
        hh_hist = hist(HH, n_H2 / V)
        
        pyplot.plot(lengths / BOND_CH, ch_hist, label='CH (CH4)')
        pyplot.plot(lengths / BOND_HH, hh_hist, label='HH (H2)')
        pyplot.legend()

        # the top two values will be the first spike, the
        # intramolecular covalent bonds
        # (one for each histogram)
        yvals = sorted(ch_hist + hh_hist)
        maxy = yvals[-3]*1.1
        #maxy = 6
        pyplot.ylim(-0.1, maxy)
        
        pyplot.plot([-1, maxbond+0.1], [1, 1], color='gray', linestyle='-')
        pyplot.xlim(-0.1, maxbond)


parser = ArgumentParser(description=__doc__)
parser.add_argument('-t', '--terminal', action='store_true', help='''
    configure graph for better terminal display
''')

parser.add_argument('path', type=str, help='''
    the file path of the .md file OR a .analysis file
''')
parser.add_argument('out', type=str, help='''
    the file to output the graph, replacing
    `$name` (filename without ext) and
    `$method` (name of graph method as input) as applicable.
    (if these aren't provided, the graphs will all overwrite each other!)
''')
parser.add_argument('options', choices=FUNCS.keys(), nargs='+', help='''
    the graph(s) to produce.
''')
parser.add_argument('--every', metavar='N', type=int, default=1, help='''
    only record data at every N timesteps
''')
parser.add_argument('--cell', metavar='path', type=str, default=None, help='''
    manually specify .cell file path
    (defaults to same name as .md file; must be provided)
''')
parser.add_argument('--hydropath', metavar='path', type=str, default=None, help='''
    path for .hydrogens.txt for .cell files
''')

# TODO: actually do stuff!

def main(argv=None):
    args = parser.parse_args(argv)
    
    if args.terminal:
        pyplot.rcParams.update({
            'font.size': 12,
            'figure.figsize': (4.4, 2.2),
            'figure.dpi': 50
        })
        
    self = Analysis(args.path, args.out, args.cell, args.hydropath, args.every)
    print(f'read {len(self.steps)} timesteps)')
        
    graphs = [FUNCS[i] for i in args.options]
    graphs = list(dict.fromkeys(graphs))  # remove duplicates
    for name in graphs:
        getattr(self, name)()
    return self

if __name__ == '__main__':
    argv = None
    argv = []
    self = main(argv)

    # todo: pickle
        