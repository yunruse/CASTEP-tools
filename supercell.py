'''
CASTEP supercell script for molecular dynamics.
By Mia Dobson.

Will modify the cell file given; use `cp src.cell dest.cell`
before running to avoid overwriting.

Best paired with a shell script, such as:

rm -r $1
cp -r source_directory $1
python3.8 supercell.py -kfcv -s $2 -p "$1/all.param" -r "$1/all.sh" $1/*.cell

where $1, a destination folder, is given, with cell size $2,
and the source directory contains an all.param and all.sh to copy through.

'''

from argparse import ArgumentParser
from os.path import isfile, split, join
from shutil import copyfile
from sys import stderr

from numpy import gcd

from parse_block import ParseError
from parse_cell import CellFile

def error(code, msg):
    print(msg, file=stderr)
    exit(code)

parser = ArgumentParser(description=__doc__)
parser.add_argument('paths', type=str, nargs='+', metavar='path',
                    help='the file path of the cell file')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='output information of operations')
parser.add_argument('-k', '--kpoint', action='store_true', help='''
replace k points with singular Baldereschi point (1/4 reciprocal lattice vectors)
''')
parser.add_argument('-f', '--forceabs', action='store_true',
                    help='force creation of positions_abs block')
parser.add_argument('-c', '--cube', action='store_true', help='''
make cell as cube-like as possible, maintaining
''')
parser.add_argument('-s', '--supercell', metavar='N', type=int, default=1,
                    help='resize the cell to be N times larger')
parser.add_argument('-N', '--number', metavar='N', type=int, default=None,
                    dest='ion_limit',
                    help='maximum number of ions permitted (requires -c)')
parser.add_argument('-p', '--param', metavar='path', type=str, default=None,
                    help='a parameter file to copy and rename for all cells')
parser.add_argument('-r', '--runscript', metavar='path', type=str, default=None,
                    help='''
a run script to copy and rename for all cells,
in which the following variables will be transformed in text:
- $NAME the name of the .cell file, sans extension
''')


def parse(args, path):
    if args.verbose:
        print('\n# Processing', path)
    if not isfile(path):
        error(2, 'fatal: path provided is not a file')       
    try:
        crystal = CellFile(path, args.forceabs)
    except ParseError as f:
        msg, line_no, line = f.args
        error(3, f'''\
> {line}
fatal error parsing cell file on line {line_no}:
{msg}''')
    except NotImplementedError as f:
        error(101, f'not yet implemented: {f.args[0]}')
        
    if args.kpoint:
        crystal.baldereschi()
    
    if args.verbose:
        N_ions = len(crystal.ions)
        print(f'Found: {N_ions} ions')
        species = crystal.species()
        if species:
            GCD = gcd.reduce(list(species.values()))
            print(f'Equivalent: {GCD} molecules of', ''.join(
                f"{s}{'' if v == GCD else v//GCD}"
                for s, v in sorted(species.items())))
    if args.cube:
        crystal.cubify(args.supercell, args.ion_limit)
    elif args.supercell > 1:
        crystal.supercell(args.supercell)
    if args.verbose:
        print(f'Result: {len(crystal.ions)} ions')

    return len(crystal.ions)

def is_file(path):
    return path is not None and path and isfile(path)

def copy_if_dest(source, dest, extension):
    if is_file(source):
        folder, name = split(dest)
        param_dest = join(folder, name.replace('.cell', extension))
        copyfile(source, param_dest)
        return param_dest

def __main__():
    args = parser.parse_args()
    num_ions = []
    for path in args.paths:
        n = parse(args, path)
        num_ions.append(n)
        
        copy_if_dest(args.param, path, '.param')
        script_path = copy_if_dest(args.runscript, path, '.sh')
        if is_file(script_path):
            with open(script_path) as f:
                content = f.read(-1)
            transforms = {
                '$NAME': split(path)[1].replace('.cell', '')
            }
            for k, v in transforms.items():
                content = content.replace(k, v)
            with open(script_path, 'w') as f:
                try:
                    f.write(content)
                except Exception as e:
                    print(type(e), e)
    if args.verbose:
        print('# Result')
        print(f'Smallest cell: {min(num_ions)} ions')
        print(f'Largest cell: {max(num_ions)} ions')
        rms = sum(i**2 / len(num_ions) for i in num_ions) ** 0.5
        print(f'Effective cell (rms): {rms:.2f} ions')
            

if __name__ == '__main__':
    __main__()
    
