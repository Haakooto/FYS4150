import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from glob import glob
import time
import subprocess
import sys
import pandas as pd
from uncertainties import ufloat


datapath = "./data/"


def analytic(T=1, M=[1,], R=1, paraRell=False):
    """
    Make table comparing numerical values against analytical ones

    Arguments:
        T: int
            Temperature
        M: list
            log10 of cycles
        R: int
            Number of cocurrent runs
    Returns:
        nicely formatted table (hopefully)
    """
    M = [int(i) for i in M]
    if paraRell:
        compiled = "paralytic.out"
        print("Running parallelized")
    else:
        compiled = "analytic.out"
        print("Running serialized")
    table = []
    cols = ["M", "$\langle \varepsilon \rangle$", "$\langle m \rangle$", "$C_v$", "$\chi$"]
    for i, m in enumerate(M):
        table.append([f"$10^{m}$"])
        t1 = time.time()
        do = subprocess.run(f"./{compiled} {T} {10**m} {R} ugly".split(" "), stdout=subprocess.PIPE)
        t2 = time.time()
        print(f"M: 10^{m}, duration: {t2 - t1}")
        result = do.stdout.decode().strip().split(" ")
        analy = result[:4]
        computed = result[4:]
        for v in range(4):
            val = ufloat(float(computed[v]), float(computed[v + 4]))
            table[-1].append(f"${val:.1u}$".replace("+/-", "\pm"))

    table.append(["Analytic"] + [f"${float(a):.4f}$" for a in analy])
    table = pd.DataFrame(table, columns=cols).transpose()
    latex = table.to_latex(escape=False, index=cols, column_format="c|" + "c"*len(M) + "|c")
    latex = latex.splitlines()
    toprule = latex.pop(4)
    latex[2] = toprule
    print("\n".join(latex))




def burntime(fname, T=1, M=1, new_run=False):
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
    M = int(M)
    fname += f"_{T}".replace(".", "")  # save temp in filename
    runname = fname
    fname += ".csv"
    file = datapath + fname
    if file not in glob(datapath + "*") or new_run:
        subprocess.Popen(f"./burn.out {runname} {T} {M}".split(" ")).wait()

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
