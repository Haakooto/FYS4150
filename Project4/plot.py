import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time


def burntime(fname, *args):
    pass


def main():
    if len(sys.argv) == 1:
        return
    elif len(sys.argv) > 2:
        args = [a for a in sys.argv[2:]]
    else:
        args = None

    func = sys.argv[1]
    if args is not None:
        exec(func + f"(*{args})")
    else:
        exec(func + "()")

if __name__ == "__main__":
    main()