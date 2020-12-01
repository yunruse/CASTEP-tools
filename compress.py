'''
CASTEP compression script.
Typically shaves about 30% of disk usage from .md and .castep
files – only by removing unnecessary space characters!

Run as `python compress.py *.{md,castep}`.
'''

from argparse import ArgumentParser

from shutil import move
from os import remove, stat


def compress(path, out):
    '''Compress a space-laden file.'''
    with open(path) as f, open(out, 'w') as fout:
        while line := f.readline():
            line = ' '.join(line.strip().split())
            fout.write(line + '\n')


def compress_inplace(path):
    tmp = path + '.compressed'
    compress(path, tmp)
    remove(path)
    move(tmp, path)


parser = ArgumentParser(description=__doc__)
parser.add_argument('paths', type=str, nargs='+', help='''
    file(s) to compress in place
''')
paths = parser.parse_args().paths

start_bytes = sum(stat(path).st_size for path in paths)

for path in paths:
    print(path)
    compress_inplace(path)

end_bytes = sum(stat(path).st_size for path in paths)

print(f'# Before: {start_bytes / 1024**2 :>5.2f} GiB')
print(f'# After:  {end_bytes / 1024**2 :>5.2f} GiB')
print(f'# Ratio:  {100*end_bytes/start_bytes :>.2f}%')
