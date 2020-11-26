import os

def files_ext(path, ext):
    for root, _, files in os.walk(path):
        for f in files:
            name, e = os.path.splitext(f)
            if e == ext:
                yield os.path.relpath(os.path.join(root, f), path)
