'''
MD parser.

Hardcoded to H2+CH4 systems: if a `hydro_path` text file is used,
H become H(mol) with the molecule as indicated by the file.
The cell will likely be a supercell, so this loops around the file.
'''

from numpy import array, NaN

# Atomic units are expressed in the Hartree system, where
# hbar = e = m_e = a_0 = 1

# hartree / hbar^2 / (m_e * a_0^2) to joules
ATOMIC_ENERGY   = 4.359744722e-18
# hbar / hartree to femtoseconds
ATOMIC_TIME     = 2.418884326e-2
# a_0 (or bohr) to angstrom
ATOMIC_LENGTH   = 0.529177210903
# hartree / k_B to kelvin
ATOMIC_TEMP     = 3.157750248e5
# hartree / a_0**3 to GPa
ATOMIC_PRESSURE = 2.942101427e4

class Ion:
    __slots__ = 'species number pos'.split()
    def __init__(self, data):
        self.species, self.number, *vec = data
        self.pos = array(vec, dtype=float) * ATOMIC_LENGTH

class Step:
    '''
    A single molecular dynamics timestep.

    Contains ion
    '''
    __slots__ = 't E T P h ions species'.split()
    def __init__(self):
        self.t = None
        self.T = NaN
        self.P = NaN
        self.E = array((NaN, NaN, NaN))
        self.h = []
        self.ions = []
        self.species = None 

    def genspecies(self):
        ss = {}
        for ion in self.ions:
            s = ion.species
            if s not in ss:
                ss[s] = []
            ss[s].append(vec)
        self.species = ss

class MDFile:
    '''
    A CASTEP .md file.
    Automatically converts to SI-like units:
    angstrom, seconds, kelvin, GPa, joule
    '''
    __slots__ = 'path file steps'.split()
    def __init__(self, path, hydro_path=None):
        self.path = path
        self.file = open(path)
        self.steps = []
        
        lines = self.file.readlines()
        for n, line in enumerate(lines):
            if line.strip() == 'END header':
                break
        else:
            return

        num_hydrogen = 0
        hydro_tags = None
        if hydro_path is not None:
            hydro_tags = []
            with open(hydro_path) as f:
                for line in f.readlines():
                    _, mol = line.strip().split()
                    mol = 'H('+mol+')'
                    hydro_tags.append(mol)
        
        step = Step()
        for line in lines[n+1:]:
            chunks = line.split()
            
            if not chunks:
                if step.t is not None:
                    if step.t % 10 == 0:
                       print(step.t)
                    self.steps.append(step)
                step = Step()
            elif len(chunks) == 1:
                step.t = float(chunks[0]) * ATOMIC_TIME
            elif len(chunks) < 3 or chunks[-2] != '<--':
                print(chunks)
                break
            else:
                *chunks, _, tag = chunks
                if tag == 'E':
                    step.E = array(chunks, dtype=float) * ATOMIC_ENERGY
                elif tag == 'T':
                    step.T = float(chunks[0])  * ATOMIC_TEMP
                elif tag == 'P':
                    step.P = float(chunks[0])
                elif tag == 'h':
                    step.h += list(map(float, chunks))
                elif tag == 'R':
                    ion = Ion(chunks)
                    if hydro_tags and ion.species == 'H':
                        ion.species = hydro_tags[num_hydrogen % len(hydro_tags)]
                        num_hydrogen += 1
                    step.ions.append(ion)
                elif tag == 'V':
                    pass  # NotImplemented
                elif tag == 'F':
                    pass  # NotImplemented
                else:
                    print(f'Unknown tag `{tag}`')
                    break

    def __del__(self):
        self.file.close()
