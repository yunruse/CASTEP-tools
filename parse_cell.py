from math import log, tau
from itertools import product

import numpy as np

from parse_block import ParseError, BlockFile

POS_REL = 'positions_frac'
POS_ABS = 'positions_abs'
CELL_VEC = 'lattice_cart'
K_POINTS = 'kpoint_list'

COORD = "{:14f} {:14f} {:14f}"
RIGHT_ANGLE = np.pi / 2


def angle_between(a, b):
    '''
    The angle in multiples of 90deg between two vectors.
    '''
    return np.arccos(np.clip(np.dot(
        a / np.linalg.norm(a),
        b / np.linalg.norm(b)
    ), -1.0, 1.0)) / RIGHT_ANGLE


def cube_difference(length, cell_vectors):
    '''
    How un-cubelike a cell is given its vectors.

    Mathematically, the sum of the square difference between
    its angle parameters and ninety degrees.
    '''
    a, b, c = cell_vectors
    A = np.linalg.norm(a)
    B = np.linalg.norm(b)
    C = np.linalg.norm(c)

    return sum((x-1)**2 for x in (
        angle_between(b, c),
        angle_between(a, c),
        angle_between(a, b),
        B/C,
        A/C,
        A/B
    ))


def shear_matrices(n):
    '''Enumerate candidate shear matrices.'''
    return list(product(
        range(1, n*2),
        range(1, n*2),
        range(1, n*2),
        range(n),
        range(n),
        range(n)
    ))


class CellFile(BlockFile):
    __slots__ = BlockFile.__slots__ + ('_ions', 'force_abs')
    '''
    A file corresponding to a (CASTEP) cell file.

    This stores its ions internally in absolute coordinates.
    '''

    def __init__(self, path, force_abs=False, lattice_save=True):
        '''
        Generate a new cell file.

        If force_abs is True, this will forcefully delete positions_rel
        and generate a positions_abs block. If False, it will
        use whichever coordinate system is given.

        If lattice_save is not set to False, the lattice_cart block
        is modified so that it ignores units in header. This means
        that the lattice_cart block will execute undefined behaviour
        if it is deleted (which you probably shouldn't do anyway)
        '''
        BlockFile.__init__(self, path)

        self.force_abs = force_abs

        # If ANG or other unit present, sneakily modify pointer
        # to point to starting from next line.
        if False and CELL_VEC in self:
            headerline = self[CELL_VEC].split('\n')[0]
            try:
                [float(x) for x in headerline.split()]
            except ValueError:
                h = len(headerline + '\n')
                s, l = self._blocks[CELL_VEC]
                self._blocks[CELL_VEC] = s+h, l-h

        # pre-convert positions to absolute if not found
        if POS_ABS not in self:
            cell = self.cell_vectors
            self._ions = [
                (sym, pos @ cell)
                for (sym, pos) in self._ions_get(POS_REL)]
        else:
            self._ions = self._ions_get(POS_ABS)

    @property
    def cell_vectors(self):
        lines = self[CELL_VEC].split('\n')
        # sometimes lines are given with units,
        # so we attempt to skip invalid headers
        # (and always assume angstroms!)
        while lines:
            try:
                [float(x) for x in lines[0].split()]
                break
            except ValueError:
                lines.pop(0)

        return np.array([
            [float(x) for x in line.split()]
            for line in lines if line
        ])

    @cell_vectors.setter
    def cell_vectors(self, cell):
        self[CELL_VEC] = '\n'.join(
            COORD.format(*vec)
            for vec in cell
        )

    def reciprocal_vectors(self):
        return tau * np.linalg.pinv(self.cell_vectors).T

    def baldereschi(self):
        '''replace k points with (1/4, 1/4, 1/4) in reciprocal space'''
        k_point = self.reciprocal_vectors().sum(1) / 4
        # x, y, z, weight
        self[K_POINTS] = COORD.format(*k_point) + ' 1\n'

    def _ions_get(self, block):
        ions = []
        for line in self[block].split('\n'):
            if line:
                sym, *pos = line.split()
                pos = np.array(list(map(float, pos)))
                if len(pos):
                    ions.append((sym, pos))
        return ions

    @property
    def ions(self):
        return self._ions

    @ions.setter
    def ions(self, value):
        self._ions = value

        if self.force_abs and POS_REL in self:
            del self[POS_REL]
        use_rel = POS_REL in self and not self.force_abs
        cellinv = np.linalg.pinv(self.cell_vectors)

        block = ''
        seen = set()
        for sym, pos in value:
            p = tuple(pos.round(3))
            if p in seen and sym == 'Q':
                print('aaaah!', sym, p)
            seen.add(p)

            if use_rel:
                pos = pos @ cellinv
            block += f'{sym:<4} {COORD.format(*pos)}\n'

        self[POS_REL if use_rel else POS_ABS] = block

    def species(self):
        species = {}
        for (s, r) in self.ions:
            if s not in species:
                species[s] = 0
            species[s] += 1
        return species

    def supercell(self, nx=1, ny=None, nz=None):
        a, b, c = self.cell_vectors

        if ny is None:
            ny = nx
            nz = nx

        # embiggen cell vectors
        self.cell_vectors = cell = np.array([
            a*nx, b*ny, c*nz
        ])

        # generate grid of ions
        ions = self.ions
        newions = []
        for x, y, z in product(range(nx), range(ny), range(nz)):
            #print(x, y, z)
            r = np.array([x/nx, y/ny, z/nz]) @ cell
            for sym, p in ions:
                newions.append((sym, r + p))
        # input('--')

        self.ions = newions

    def cubify_old(self, n=1):
        # outdated assumption of method:
        # map unit cube to relative coordinates and round
        # this generates a supercell with shear, which
        # maintains translational symmetry

        cell = self.cell_vectors
        vec = (n**2 * np.linalg.pinv(cell)).round()
        self.supercell(int(vec.max()))
        self.cell_vectors = vec @ cell
        self.wrap_around()

    def cubify(self, n=1, ion_limit=None):
        '''
        Produce a supercell sheared so as to be closest to right angles.
        As it stands, this is not necessarily closest to a cube;
        the cell will be repeated n times in all three axes.
        '''
        # iterate over 8n**3 possible lower diagonal shears
        # to find that with the minimal angle

        cell = self.cell_vectors

        def shear(a, b, c, d, e, f):
            return np.array([
                [a, 0, 0],
                [d, b, 0],
                [e, f, c]
            ])

        shears = shear_matrices(n)
        if ion_limit is not None:
            max_supersize = ion_limit // len(self.ions)
            shears = [
                (a, b, c, *d) for (a, b, c, *d) in shears
                if a*b*c <= max_supersize
            ]

        shears.sort(key=lambda sh: cube_difference(
            n, shear(*sh) @ cell
        ))

        sh = shear(*shears[0])
        self.supercell(*np.diagonal(sh))
        self.cell_vectors = sh @ cell
        self.wrap_around()

    def wrap_around(self):
        # For each cell, move to relative coordinates, then 'wrap about'.
        cell = self.cell_vectors
        cellinv = np.linalg.pinv(cell)
        self.ions = [
            (sym, ((pos @ cellinv) % 1) @ cell)
            for (sym, pos) in self.ions]
