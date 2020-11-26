class ParseError(ValueError):
    pass


def char_line_content(file):
    file.seek(0)
    chars = 0
    lines = 0
    line = ''
    while c := file.read(1):
        if c == '\n':
            yield chars, lines, line
            lines += 1
            line = ''
        else:
            line += c
        chars += 1


_START = '%BLOCK {}\n'
_END = '%ENDBLOCK {}\n'


def find_blocks(file):
    '''
    given a cell file return a map of blocks.

    mapping is {blockname: (start, length)}
    such that file.seek(start) and file.read(length)
    produce the block's exact contents, sans headers'''
    blocks = {}
    block = None
    blockline = 0
    blockchar = 0

    for stream, line_no, line in char_line_content(file):
        line = line.strip()
        if line.startswith('%'):
            cmd, name, *rest = line[1:].split()
            cmd = cmd.lower()

            if rest:
                raise ParseError(
                    'cannot parse this block structure',
                    line_no, line)

            elif cmd == 'block':
                if block is not None:
                    raise ParseError(
                        f'block {name!r} started before {block!r} could end',
                        line_no, line)
                if block in blocks:
                    raise ParseError(
                        f'block {name} started but was already declared')
                block = name
                blockstart = stream + 1

            elif cmd == 'endblock':
                if block != name:
                    raise ParseError(
                        f'block {name!r} ended but never started',
                        line_no, line)
                blockend = stream - len(line) - 1
                # We assume there are two newlines we are ignoring,
                # separating the header start and end lines from contents.
                # If there is only one newline, as the block is empty,
                # the length is -1. This is intentional, and is caught.
                blocks[block] = (blockstart, blockend - blockstart)
                block = None

            else:
                raise ParseError(
                    f"cannot understand {name!r} instruction",
                    line_no, line)

    if block is not None:
        raise ParseError(
            f'block {block!r} started on line {blockstart} but never ended',
            line_no, line)

    return blocks


class BlockFile:
    '''
    A file, which stores the coordinates of all of its blocks,
    allowing blocks to be modified in-place without disturbing
    outside content.

    These blocks may be accessed by indexing.
    '''

    __slots__ = ('_file', '_blocks')

    def __init__(self, path, lattice_save=True):
        self._file = open(path, 'r+')
        self._blocks = find_blocks(self._file)

    def __getitem__(self, name):
        start, length = self._blocks.get(name, (0, 0))
        if length < 0:
            return ''
        f = self._file

        f.seek(start)
        return f.read(length)

    def __delitem__(self, name):
        if name not in self._blocks:
            return
        start, length = self._blocks.pop(name)
        f = self._file

        # this assumes there is no extr. whitespace in block declarations!

        f.seek(start + length + len(_END.format(name)))
        pos_end = f.tell()
        rest_of_file = f.read()
        f.seek(start - len(_START.format(name)))
        pos_start = f.tell()
        f.truncate()
        f.write(rest_of_file)

        # pull all indices to correct values
        for n, (s, l) in self._blocks.items():
            if s > start:
                self._blocks[n] = (s - (pos_end - pos_start), l)

    def __setitem__(self, name, contents):
        f = self._file

        if name in self:
            start, length = self._blocks[name]
            length2 = len(contents)

            was_empty = length < 0
            if was_empty:
                length = 0
                contents = contents + '\n'

            # r+ mode can only write to end, so do some cutting
            f.seek(start + length)
            rest_of_file = f.read()
            f.seek(start)
            f.truncate()
            f.write(contents)
            f.write(rest_of_file)

            self._blocks[name] = (start, length2)

            # push all indices to correct values
            offset = length2 - length
            for n, (s, l) in self._blocks.items():
                if s > start:
                    self._blocks[n] = (s + offset, l)
        else:
            # not found; append to end
            f.write('\n\n' + _START.format(name))
            start = f.tell()
            f.write(contents)
            length = f.tell() - start
            f.write(_END.format(name) + '\n')
            self._blocks[name] = (start, length)

    def __contains__(self, name):
        return name in self._blocks

    def __del__(self):
        self._file.close()
