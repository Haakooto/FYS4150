import numpy as np
import plotly.express as px
import matplotlib.pyplot as plt
from glob import glob
import subprocess
import sys
import pandas as pd

datapath = "./data/"


def burntime(fname, M=0, new_run=False):
    runname = fname
    if fname[-4:] != ".csv":
        fname += ".csv"
    file = datapath + fname
    if file not in glob(datapath + "*") or new_run:
        subprocess.Popen(f"./burn {runname} {M}".split(" ")).wait()

    data = pd.read_csv(file, header=0, sep="  ")
    print(data["er1"])


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
        print(args)
        func = sys.argv[1]
        if args is not None:
            exec(func + f"(*{args})")
        else:
            exec(func + "()")
    except:
        print("Bad usage! Pass 'help' from commandline to see guide")
        raise
if __name__ == "__main__":
    main()