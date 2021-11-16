import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from glob import glob
import subprocess
import sys
import pandas as pd

datapath = "./data/"


def burntime(fname, T=1, M=0, new_run=False):
    """
    Plot the energy and magnetization for the first M cycles
    Determine burntime

    Arguments:
        fname: str
            filename of datafile. If file already in data-directory,
            read from that. Else, call ./burn to make one with the given
            configuration. If new_run, overwrite the file with new data to use
        T: int
            temperature
        M: int
            MC-cycles
        new_run: bool
            force new run
    Returns:
        lines and dots with shiny colours
    """
    runname = fname
    if fname[-4:] != ".csv":
        fname += ".csv"
    M = int(M)
    file = datapath + fname
    if file not in glob(datapath + "*") or new_run:
        subprocess.Popen(f"./burn {runname} {T} {M}".split(" ")).wait()

    data = pd.read_csv(file, header=0, sep=",")
    data = data.head(M)

    fig = go.Figure()
    colors = px.colors.qualitative.Plotly
    c = -1
    for init in ["rnd", "low", "hig"]:
        c += 1
        for quant in ["e", "m"]:
            q = "energy" if quant == "e" else "magnetization"
            raw = f"{quant}_{init}"
            avg = "avg_" + raw
            name = q + " with " + init
            fig.add_trace(go.Scatter(y=data[avg], mode="lines", line=dict(width=4, color=colors[c]), name=name))
            fig.add_trace(go.Scatter(y=data[raw], mode="markers", marker=dict(size=4, color=colors[c]), showlegend=False))
    fig.update_layout()
    fig.show()

def main():
    if "help" in sys.argv:
        print("First argument after 'python plot.py' must be a function in this file")
        print("If this function takes any arguments, the are to be passed as well")
        print("Example: 'python plot.py burntime burn 2.4 10")

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
        raise

if __name__ == "__main__":
    main()