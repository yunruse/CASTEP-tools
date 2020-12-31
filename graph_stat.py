from math import isnan
import os

import numpy as np

from parse_md import MDFile
from parse_cell import CellFile

record_every = 200

print('name,n,V,P,T')

for mdpath in os.listdir():
    if mdpath.endswith('.md'):
        name, _ = mdpath.split('.')
        md = MDFile(
            mdpath,
            step_is_valid=lambda i: i % record_every == 0
        )
        a, b, c = CellFile(name + '.cell').cell_vectors
        
        n = len(md.steps[0].ions)
        V = np.linalg.norm(np.dot(np.cross(a, b), c))
        
        Ps = [step.P for step in md.steps if not isnan(step.P)]
        Ts = [step.T for step in md.steps if not isnan(step.P)]
        P = sum(Ps) / len(Ps)
        T = sum(Ts) / len(Ts)
        print(f'{name},{n},{V},{P},{T}')
