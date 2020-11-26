

def old_cell_dist(a, b, cell, cellinv=None):
    # Assumes orthorhombic, so incorrect
    if cell is None:
        return ((b - a)**2).sum()
    # change to relative coordinates (a 1x1x1 cube)
    ai = a @ cellinv
    bi = b @ cellinv
    # with a temp change of coordinates,
    # ensure b is inside a cube with a at the centre
    zero = ai - 0.5
    bi = ((bi - zero) % 1) + zero
    return square_dist(a, bi @ cell)
    
def old_cell_dist_2(a, b, cell, cellinv=None):
    '''
    Return square_dist(a, b) considering that
    a and b may have "wrapped around" a cell
    (i.e. pick the smallest from all combinations)
    '''
    if cell is None:
        return ((b - a)**2).sum()
    ai = a @ cellinv
    bi = b @ cellinv
    D = array([-1, 0, 1])
    return min(
        square_dist(a, ((bi + delta) % 1) @ cell)
        for delta in product(D, D, D)
    )

OPT_BOND_H = 'bond_h'
def bond_h(md, cell, cellinv):
    pyplot.title(f'Hydrogen bond length change of `{md.path}`')
    pyplot.xlabel('$n$th smallest hydrogen bond')
    pyplot.ylabel('Bond length / Å')
    
    ions = [x for s, n, x in md.steps[0].ions if s == 'H']
    bonds = [
        cell_square_dist(ions[a], ions[b], cell, cellinv)
        for a in range(len(ions)) for b in range(a)
    ]
    bonds.sort()
    pyplot.plot(bonds)

OPT_BOND_CH = 'bond_ch'
def bond_okd_ch(md, cell, cellinv):
    pyplot.title(f'C-H bond length change of `{md.path}`')
    pyplot.xlabel('$n$th smallest C-H bond')
    pyplot.ylabel('Bond length / Å')
    
    C = [i.pos for i in md.steps[0].ions if i.species == 'C']
    H = [i.pos for i in md.steps[0].ions if i.species == 'H']
    bonds = [
        cell_square_dist(c, h, cell, cellinv)
        for c in C for h in H
    ]
    bonds.sort()

    N = len(C) * 4
    bonds = bonds[:N*10]
    pyplot.plot(bonds)
    pyplot.plot([N, N], [0, max(bonds)])

def bond_ch(md, cell, cellinv):
    pyplot.title(f'C-H bond length change of `{md.path}`')
    pyplot.xlabel('$n$th smallest C-H bond')
    pyplot.ylabel('Bond length / Å')
    
    C = [i.pos for i in md.steps[0].ions if i.species == 'C']
    H = [i.pos for i in md.steps[0].ions if i.species == 'H']
    for c in C:
        bonds = [
            cell_square_dist(c, h, cell, cellinv) for
            h in H
        ]
        bonds.sort()
        pyplot.plot(bonds[:20], alpha=0.05)
