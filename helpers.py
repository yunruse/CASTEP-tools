from os.path import expanduser


def find(path):
    path = expanduser(path)
    if not isfile(path):
        raise FileNotFoundError(f'Could not find `{path}`')


class ParseError(ValueError):
    pass
