import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import subprocess
import sys
import pandas as pd

datapath = "./data/"

def burntime(fname, M, new_run=False):
    M = int(M)
    if fname in glob(datapath):
        pass


def main():
    if "help" in sys.argv:
        print("First argument after 'python plot.py' must be a function in this file")
        print("If this function takes any arguments, the are to be passed as well")
        print("Example: 'python plot.py burntime burn 10")

    try:
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
    except:
        print("Bad usage! Pass 'help' from commandline to see guide")

if __name__ == "__main__":
    main()