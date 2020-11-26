import os

from common import files_ext

def get_max_iter(path):
    n = 0
    with open(path) as f:
        for l in f.readlines():
            if l.startswith(' Starting MD iteration'):
                try:
                    n = max(n, int(l.split()[3]))
                except ValueError:
                    pass
    return n

def main(path):
    for f in files_ext(path, '.castep'):
        iters = get_max_iter(f)
        print(f'{f:20}: {iters}')

if __name__ == '__main__':
    PATH = os.path.expanduser('~/Scratch/')
    PATH = os.path.expanduser('~/Documents/University Work/SH Project/lab')
    main(PATH)
