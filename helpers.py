from os.path import expanduser, isfile


def find(path):
    path = expanduser(path)
    if not isfile(path):
        raise FileNotFoundError(f'Could not find `{path}`')
    return path


def rolling(vals, n):
    '''
    Return a rolling average. For each value x, this is the average
    of (2n+1) points, centred on the value x.
    '''
    rolling = []
    for i, v in enumerate(vals):
        vals_to_left = i
        vals_to_right = (len(vals) - 1) - i
        size = min(vals_to_left, vals_to_right, int(n/2))
        roll = vals[i-size:i+size+1]
        rolling.append(sum(roll) / len(roll))
    return rolling


class ParseError(ValueError):
    pass
