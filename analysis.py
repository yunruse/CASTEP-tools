'''
Analysis script by Mia Dobson.
Outputs graph(s) of desired statistics.
'''

from pickle import dump
from argparse import ArgumentParser
from functools import wraps
from itertools import product
from os.path import isdir, split, splitext

from matplotlib import pyplot
import numpy as np
from numpy import dot, cross
from numpy.linalg import pinv, norm

from helpers import find
from parse_cell import CellFile
from parse_md import MDFile


FUNCS = {}


def method(command, anim=False):
    '''
    Meta-wrapper for functions in Analysis.

    Registers an identifier, used both in commands and as a key
    in Analysis.graphs. Additionally, saves output to files.
    '''
    # and adding niceties
    def wrapper(func):
        if anim:
            def newfunc(self):
                for i, step in enumerate(self.steps):
                    pyplot.clf()
                    ax = pyplot.axes()
                    data = func(self, ax, step) or dict()
                    data['t'] = round(step.t, 1)
                    self.output(command, ax, data)
        else:
            def newfunc(self):
                pyplot.clf()
                ax = pyplot.axes()
                data = func(self, ax)
                self.output(command, ax, data)

        FUNCS[command] = func.__name__
        return wraps(func)(newfunc)
    return wrapper


D = [-1, 0, 1]
DELTAS = np.array(list(product(D, D, D)))

TIMELABEL = 'Time / ps'


class Analysis:
    def __init__(
        self,
        mdpath: str,
        outpath: str,
        cellpath=None,
        hydro_path=None,
        record_every: int = 1,
        save_state=False,
    ):
        find(mdpath)
        _, name = split(mdpath)
        self.name, _ = splitext(name)
        self.outpath = outpath
        self.save_state = save_state

        self.record_every = record_every

        md = MDFile(
            mdpath,
            hydro_path,
            step_is_valid=lambda i: i % record_every == 0,
        )
        self.steps = [s for s in md.steps if s.ions]

        if cellpath is None:
            cellpath = mdpath.replace('.md', '.cell')
        elif cellpath.lower() == 'none':
            cellpath = None

        if cellpath is not None:
            self.has_hydro_tags = hydro_path is not None

            find(cellpath)
            self.cell = CellFile(cellpath).cell_vectors
            self.cellinv = pinv(self.cell)

    def output(self, command, ax, data=None):
        '''Handler for methods. May be called in own graphs.'''
        if data is NotImplemented:
            return
        if data is None:
            data = dict()

        data['name'] = self.name
        data['method'] = command

        if self.outpath:
            path = self.outpath
            for k, v in data.items():
                path = path.replace('$'+k, str(v))
            print(path)
            pyplot.savefig(path)
            if self.save_state:
                path, _ = splitext(path)
                path += '.ax'
                print(path)
                with open(path, 'wb') as f:
                    dump(ax, f)

    def plot_variable(self, ax, name, unit, func):
        ax.set_title(f'{name} of `{self.name}`')
        ax.set_xlabel(TIMELABEL)
        ax.set_ylabel(f'{name} / {unit}')

        time = [step.t / 1000 for step in self.steps]
        var = [func(step) for step in self.steps]
        valid = [i for i in var if i is not np.nan]
        ax.plot(time, var)
        ax.plot([0, max(time)], [valid[0], valid[0]],
                color='gray', linestyle='-')

    @method('temp')
    def temperature(self, ax):
        self.plot_variable(ax, 'Temperature', 'K', lambda step: step.T)

    @method('pressure')
    def pressure(self, ax):
        self.plot_variable(ax, 'Pressure', 'as raw', lambda step: step.P)

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
    def msds(self, ax):
        ax.set_title(f'Mean square displacement of `{self.name}`')
        steps = self.steps[1:]

        ax.set_xlabel(TIMELABEL)
        t = [step.t / 1000 for step in steps]

        ax.set_ylabel('Mean square displacement / Å²')
        msds = [self.msd(step) for step in steps]
        for k in msds[0]:
            ax.plot(t, [i[k] for i in msds], label=k)

        ax.legend()

    def bonds_step(self, step):
        C = [i.pos for i in step.ions if i.species == 'C']
        H_H2 = [i.pos for i in step.ions if i.species == 'H(H2)']
        H_CH4 = [i.pos for i in step.ions if i.species == 'H(CH4)']

        CH = np.array([
            self.cell_square_dist(c, h) ** 0.5
            for c in C for h in H_CH4
        ])
        HH = np.array([
            self.cell_square_dist(H_H2[a], H_H2[b]) ** 0.5
            for a in range(len(H_H2)) for b in range(a)
        ])
        CC = np.array([
            self.cell_square_dist(C[a], C[b]) ** 0.5
            for a in range(len(C)) for b in range(a)
        ])

        return CC, CH, HH

    def bonds_average(self, dstep):
        steps = self.steps
        # TODO: look after 1ps
        CC, CH, HH = self.bonds_step(steps[0])
        N = 0
        for i in np.arange(dstep, len(steps), dstep):
            N += 1
            cc, ch, hh = self.bonds_step(steps[i])
            CC += cc
            CH += ch
            HH += hh

        CC = [i / N for i in CC]
        CH = [i / N for i in CH]
        HH = [i / N for i in HH]
        return CC, CH, HH

    @method('bond')
    def bonds_graph(self, ax):
        ax.set_title(f'Bond lengths of `{self.name}`')
        CC, CH, HH = self.bonds_average(dstep=100)

        ax.set_xlabel('Bond length $r$ (Angstroms)')
        ax.set_xlim(-0.1, 10)
        ax.set_ylabel('Cumulative number of bonds')
        #ax.set_ylim(-100, len(self.steps[0].ions)**2)

        ax.plot(sorted(CC), range(len(CC)), label='C-C (CH4)')
        ax.plot(sorted(CH), range(len(CH)), label='C-H (CH4)')
        ax.plot(sorted(HH), range(len(HH)), label='H-H (H2)')
        ax.legend()

    @method('bondanim', anim=True)
    def bonds_anim(self, ax, step):
        ax.set_title(f'Bond lengths of `{self.name}` (t={step.t}')
        CC, CH, HH = self.bonds_step(step)

        ax.set_xlabel('Bond length $r$ (Angstroms)')
        ax.set_xlim(-0.1, 10)
        ax.set_ylabel('Cumulative number of bonds')
        #ax.set_ylim(-100, len(self.steps[0].ions)**2)

        ax.plot(sorted(CC), range(len(CC)), label='C-C (CH4)')
        ax.plot(sorted(CH), range(len(CH)), label='C-H (CH4)')
        ax.plot(sorted(HH), range(len(HH)), label='H-H (H2)')
        ax.legend()

    @method('rdf')
    def rdf(self, ax):
        DR = 0.01
        DSTEP = 100
        EPSILON = 1e-15

        a, b, c = self.cell
        V = norm(np.dot(np.cross(a, b), c))
        maxbond = min(norm(a), norm(b), norm(c))
        lengths = np.arange(DR, maxbond, DR)

        if not self.has_hydro_tags:
            print('Error! rdf requires hydrogen tagging')
            return NotImplemented

        ax.set_title(f'Radial distribution function of `{self.name}`')
        ax.set_xlabel(f'Bond length $r$ (dstep={DSTEP})')
        ax.set_ylabel('Radial distribution $g(r)$')

        CC, CH, HH = self.bonds_average(DSTEP)

        def hist(bonds, rho):
            bonds = sorted(bonds, reverse=True)
            hist = []
            for r in lengths:
                n = 0
                while bonds and bonds[-1] <= r + EPSILON:
                    bonds.pop()
                    n += 1
                hist.append(n / (4/3*np.pi*((r+DR)**3 - r**3)) / rho)
            return hist

        ions = self.steps[0].ions

        n_H2 = len([i.pos for i in ions if i.species == 'H(H2)'])
        n_CH4 = len([i.pos for i in ions if i.species == 'H(CH4)'])
        n_C = len([i.pos for i in ions if i.species == 'C'])

        cc_hist = hist(CC, n_C**2 / V)
        ch_hist = hist(CH, n_C * n_CH4 / V)
        hh_hist = hist(HH, n_H2**2 / V)

        ax.plot(lengths, cc_hist, label='C-C (CH4)')
        ax.plot(lengths, ch_hist, label='C-H (CH4)')
        ax.plot(lengths, hh_hist, label='H-H (H2)')
        ax.legend()

        # the top two values will be the first spike, the
        # intramolecular covalent bonds
        # (one for each histogram)
        yvals = sorted(ch_hist + hh_hist)
        maxy = yvals[-3]*1.1
        # maxy = 6
        ax.set_ylim(-0.1, maxy)

        ax.plot([-1, maxbond+0.1], [1, 1], color='gray', linestyle='-')
        ax.set_xlim(-0.1, maxbond)


parser = ArgumentParser(description=__doc__)
parser.add_argument('-t', '--terminal', action='store_true', help='''
    configure graph for better terminal display
''')
parser.add_argument('-s', '--savestate', action='store_true', help='''
    save state to an .analysis file (with same name scheme as `out`)
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
    path for .cell file. Defaults to same name as .md file.
    Use `--cell none` to provide no cell
    (which some graphs will explicitly fail on!)
''')
parser.add_argument('--hydropath', metavar='path', type=str, default=None, help='''
    path for .hydrogens.txt for .cell files
''')


def main(argv=None):
    args = parser.parse_args(argv)

    if args.terminal:
        pyplot.rcParams.update({
            'font.size': 12,
            'figure.figsize': (4.4, 2.2),
            'figure.dpi': 50
        })

    self = Analysis(
        args.path, args.out, args.cell,
        args.hydropath, args.every, args.savestate)
    print(f'read {len(self.steps)} timesteps')

    graphs = [FUNCS[i] for i in args.options]
    graphs = list(dict.fromkeys(graphs))  # remove duplicates
    for name in graphs:
        method = getattr(self, name)
        method()
    return self


if __name__ == '__main__':
    argv = None
    self = main(argv)

    # todo: pickle
