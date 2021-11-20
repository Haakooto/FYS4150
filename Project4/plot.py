import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from glob import glob
import time
import subprocess
import sys
import pandas as pd
from uncertainties import ufloat
import matplotlib.pyplot as plt


datapath = "./cdata2/"


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
    for init in ["rnd", "low"]:
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


def run_temps(fname, Tmin=1, Tmax=2, Ts=2, M=1, R=1, new_runs=False, L=[40,60,80,100]):
    """
    Make datafile for each L in temperature analysis.
    Run this func as if it plots the result, while
    in fact is only makes sure the data files are there.
    It then call on plot_temps() to do the actual plotting.
    This is just to separate the tasks-
    It is not nessisary to understand how this function works.
    It works.

    Arguments:
        fname: str
            filename for run. datafiles will get added L and csv
        Tmin: float
            start temp
        Tmax: float
            end temp, inclusive
        Ts: int
            Number of temps to run
        M: int
            log10 of cycles
        R: int
            Number of cocurrent runs
        new_run: bool
            wether to force new run
        L: list of ints
            Lattice sizes to run for
    Returns:
        lines and dots with shiny colours
    """
    if type(L) == str:
        L = eval(L)
    if type(new_runs) == str:
        new_runs = eval(new_runs)
    Ls = [int(i) for i in L]
    Lruns = {}  # make dict with all L
    for L in Ls:
        Lruns[L] = {}  # each element is a dict
        runname = fname + f"_{L}"
        name = runname + ".csv"
        file = datapath + name
        Lruns[L]["data"] = file
        if file not in glob(datapath + "*") or new_runs:
            print(f"Starting run with L = {L}")
            Lruns[L]["start"] = time.time()
            Lruns[L]["process"] = subprocess.Popen(f"./tempting.out {runname} {int(M)} {int(R)} {L} {float(Tmin)} {float(Tmax)} {Ts}".split(" "))#.wait()
            Lruns[L]["done"] = None
            # print(f"{L} now done. Time spent: {time.time() - Lruns[L]['start']}")

        else:
            print(f"There was already data for L = {L}")
            Lruns[L]["done"] = True  # mark as done
            Lruns[L]["start"] = time.time()

    done = np.ones(len(Ls))
    while sum(done):
        for L in Ls:
            if Lruns[L]["done"] is None:
                Lruns[L]["done"] = Lruns[L]["process"].poll()
            else:
                if done[Ls.index(L)] != 0:
                    print(f"{L} now done. Time spent: {time.time() - Lruns[L]['start']}")
                    done[Ls.index(L)] = 0
    print(f"Data now ready for all L")
    Ts = np.linspace(float(Tmin), float(Tmax), int(Ts))
    plot_temps(Lruns, Ts)

def plot_temps(Ls, Ts):
    """
    Does the actual plotting of temperature plots.
    Do not run this func directly, rather call upon run_temps()
    See that function for docstring

    The errors in the values have not yet been accounted for, and are assumed to be 0
    """
    e = np.zeros((len(Ls), len(Ts)))
    m = np.zeros((len(Ls), len(Ts)))
    Cv = np.zeros((len(Ls), len(Ts)))
    chi = np.zeros((len(Ls), len(Ts)))
    l = []
    for i, (L, info) in enumerate(Ls.items()):
        l.append(L)
        data = pd.read_csv(info["data"], header=1, sep=",")
        e[i] = data["e_avg"]
        m[i] = data["m_avg"]
        Cv[i] = data["Cv"]
        chi[i] = data["chi"]
    print(e.shape)
    print(l)
    e = pd.DataFrame(e.T, columns=l)
    m = pd.DataFrame(m.T, columns=l)
    Cv = pd.DataFrame(Cv.T, columns=l)
    chi = pd.DataFrame(chi.T, columns=l)
    print(data)
    print(e)
    print(m)
    print(Cv)
    print(chi)
    plt.plot(Ts, e)
    plt.show()
    plt.plot(Ts, m)
    plt.show()
    plt.plot(Ts, Cv)
    plt.show()
    plt.plot(Ts, chi)
    plt.show()



def plot_pdf():

    data_T_low = pd.read_csv("data/pdf_T1.csv", header = 0, sep = ",")
    data_T_high = pd.read_csv("data/pdf_T2.4.csv", header = 0, sep = ",")

    fig_low = px.histogram(
        x=data_T_low["e_avg"],
        y=data_T_low["prob"],
        nbins=250)

    fig_low.update_layout(
        #xaxis_range=[-2, 0.2],   #all p=0 after e=0.2
        #yaxis_range=[0,0.9],
        font_family="Open sans",
        font_size=30,
        title="Histogram of measured probability distribution of the energy for T = 1",
        xaxis_title=r"$ \huge \text{Energy  }  \epsilon$",
        yaxis_title="Probability",
        legend=dict(yanchor="top", xanchor="left", x=0.01, y=0.99))




    fig_high = px.histogram(
        x=data_T_high["e_avg"],
        y=data_T_high["prob"],
        nbins=250)

    fig_high.update_layout(
        #xaxis_range=[-2, 0.2],   #all p=0 after e=0.2
        #yaxis_range=[0, 0.9],
        font_family="Open sans",
        font_size=30,
        title="Histogram of measured probability distribution of the energy for T = 2.4",
        xaxis_title=r"$ \huge \text{Energy  }  \epsilon$",
        yaxis_title="Probability",
        legend=dict(yanchor="top", xanchor="left", x=0.01, y=0.99))

    fig_low.show()
    fig_high.show()


def paralympics():
    """
    MC STOP. Hammer time!

    Arguments:
        None. All shit is hard coded
    Returns:
        Speed-up factor
    """

    M = 2000
    R = 100
    L = 20
    T = 2
    A = 10
    for p in ["nopara", "para"]:
        for o in ["time", "optime"]:
            times = np.zeros(A)
            for i in range(A):
                run = subprocess.run(f"./{o}.out {M} {R} {L} {T} {p}".split(" "), stdout=subprocess.PIPE)#.wait()
                times[i] = float(run.stdout.decode().strip())
            print(p, o, np.mean(times))


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
    # plot_pdf()
