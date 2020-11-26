'''
Animation tool for hydrocarbons.

Produces snapshots at various frames, and may color-code hydrogens.
'''

from os.path import isfile, isdir, split, splitext

import numpy as np
from matplotlib import pyplot
import matplotlib.animation as animation

from parse_cell import CellFile
from parse_md import MDFile

ATOMS = {
    'H': (1, 'white'),
    'C': (12, 'gray'),
}

H_COLORS = ['b', 'g', 'r', 'c', 'm', 'y']

# Sphere parameters
RADIUS = 0.2
N_SPHERE = 15
psi = np.linspace(0, 2 * np.pi, N_SPHERE*2)
theta = np.linspace(0, np.pi, N_SPHERE)
SPHERE = np.array([
    np.outer(np.cos(psi), np.sin(theta)),
    np.outer(np.sin(psi), np.sin(theta)),
    np.outer(np.ones(np.size(psi)), np.cos(theta))
]) * RADIUS


class Animator(MDFile):
    def __init__(
        self, path, savepath,
        start, stop, step,
        steps_per_picosecond=2000,
    ):
        self.times = np.arange(start, stop, step)
        steps = (self.times * steps_per_picosecond).astype(int)
        max_step = max(steps)
        MDFile.__init__(
            self, path,
            step_is_valid=lambda i: i in steps,
            step_should_stop=lambda i: i > max_step
        )
        self.savepath = savepath
        self.steps_per_picosecond = steps_per_picosecond
        
        _, name = split(path)
        self.name, _ = splitext(name)
        self.fig = pyplot.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        
        maxs, mins = np.array([0, 0, 0]), np.array([np.inf, np.inf, np.inf])
        self.allowed = set()
        first_step = steps[0]
        for ion in self.steps[first_step].ions:
            mins = np.min((mins, ion.pos), axis=0)
            maxs = np.max((maxs, ion.pos), axis=0)
            x, y, z = ion.pos
            #if min(x, y, z) > 4 and max(x, y, z) < 8:
            #    self.allowed.add((ion.number, ion.species))
        self.minmax = (mins, maxs)


    def clear(self):
        self.ax.clear()
        self.ax.set_title(f'MD of `{self.name}` (Angstroms)')
        mins, maxs = self.minmax
        self.ax.set_xlim3d([mins[0], maxs[0]])
        self.ax.set_ylim3d([mins[1], maxs[1]])
        self.ax.set_zlim3d([mins[2], maxs[2]])
        # todo: make spheres actually look spherical?
        # todo: display axes
    
    def plot(self, t=0):
        '''Basic plot of atoms.'''
        i = int(t * self.steps_per_picosecond // 1)
        
        self.clear()
        artists = []
        ions = self.steps[i].ions
        # todo: something to make it more... visible?
        # ions = ions[:len(ions)//8]
        for ion in ions:
            x, y, z = ion.pos
            #if (ion.number, ion.species) not in self.allowed:
            #    continue
            n_e, col = ATOMS[ion.species]
            if ion.species == 'H':
                col = H_COLORS[int(ion.number) % len(H_COLORS)]
            sx, sy, sz = SPHERE * n_e ** (1/3)
            artist = self.ax.plot_surface(x+sx, y+sy, z+sz, color=col)
            artists.append(artist)

    def run(self, start=None, end=None, step=None):
        if start is not None:
            times = np.arange(start, stop, step)
        else:
            times = self.times
        for t in times:
            self.plot(t)
            path = self.savepath
            path = path.replace('$t', f'{t:06.2f}')
            path = path.replace('$n', self.name)
            pyplot.savefig(path)
        

if __name__ == '__main__':
    # TODO: make this a cli!
    DIR = '/System/Volumes/Data/Volumes/Media/castep/cells/'
    PATH = DIR + 'CH6_P212121-100gpa-600K.md'

    PATH = DIR + 'CH6_P212121-20gpa-300K.md'
    SAVEPATH = '/Users/yunruse/Downloads/test/anim-$n-$t.png'
    start, stop, step = 7.9, 8.65, 0.01

    self = Animator(
        PATH, SAVEPATH,
        start=start, stop=stop, step=step,
        steps_per_picosecond=2000,
    )
    self.run()
