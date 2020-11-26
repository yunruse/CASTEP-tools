import pandas

def parse_md(text):
    
    timesteps = []
    t = None
    data = {}
    i = 0  # position in current timechunk
    for line in text.split('\n'):
        i += 1
        line = line.strip()
        if line == 'END header':
            t = 0
            continue
        elif t is None:
            continue
        
        if not line:
            i = 0
        elif i == 1:
            data['t'] = t
            timesteps.append(data)
            data = dict()
            t = float(line)
        else:
            *rest, arrow, tag = line.split()
            assert arrow == '<--'
            if tag in ('T', 'P'):
                data[tag] = float(rest[0])
            
            elif tag == 'E':
                data['E'], data['E_h'], data['E_k'] = map(float, rest)
            
            elif tag in ('R', 'V', 'F'):
                chem, num, *rest = rest
                rest = map(float, rest)
                
                # store details about ions
            
            elif tag in ('h', 'hv'):
                chem, num, *rest = rest
                rest = map(float, rest)
                # store details about cells
                
    data = pandas.DataFrame(timesteps[1:])
    data = data.set_index('t')

    return data

if __name__ == '__main__':
    from sys import argv
    _, infile, outfile = argv
    with open(infile) as f:
        data = f.read()
    data = parse_md(data)
    from matplotlib import pyplot
    pyplot.rcParams.update({
        'font.size': 12,
        'figure.figsize': (4.4, 2.2),
        'figure.dpi': 50
    })
    data['T'].plot()
    pyplot.savefig(outfile)
