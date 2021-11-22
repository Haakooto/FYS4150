import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from glob import glob
import time
import subprocess
import sys
import pandas as pd
from uncertainties import ufloat
from scipy.interpolate import UnivariateSpline
from scipy.stats import linregress


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
    fname += f"_{T}"
    runname = fname
    fname += ".csv"
    file = datapath + fname
    if file not in glob(datapath + "*") or new_run:
        subprocess.Popen(f"./burn.out {runname} {T} {M}".split(" ")).wait()

    data = pd.read_csv(file, header=0, sep=",")
    data = data.head(M)

    fig1 = go.Figure()
    colors = px.colors.qualitative.Plotly
    c = -1
    for init, method in zip(["rnd", "low"], ["Random ", "Low energy "]):
        raw = f"e_{init}"
        c += 1
        name = method + " spin initialisation"
        fig1.add_trace(go.Scatter(y=data[raw], mode="markers", marker=dict(size=7, color=colors[c]), name=name))

    fig1.update_layout(
        font_family="Open sans",
        font_size=30,
        title=f"Mean energy as function of MC cycles for a 20 x 20 lattice at T = {T}",
        xaxis_title="MC cycles",
        yaxis_title="Mean energy [J]",
        legend=dict(yanchor="top", xanchor="right", x=0.99, y=0.99))


    fig2 = go.Figure()
    c = -1
    for init, method in zip(["rnd", "low"], ["Random ", "Low energy "]):
        raw = f"m_{init}"
        c += 1
        name = method + "spin initialisation"
        fig2.add_trace(go.Scatter(y=data[raw], mode="markers", marker=dict(size=7, color=colors[c]), name=name))


    fig2.update_layout(
        font_family="Open sans",
        font_size=30,
        title=f"Mean magnetisation as function of MC cycles for a 20 x 20 lattice at T = {T}",
        xaxis_title="MC cycles",
        yaxis_title="Mean  magnetisation",
        legend=dict(yanchor="bottom", xanchor="right", x=0.99, y=0.01, font_size=30))

    fig1.show()
    fig2.show()



def run_temps(fname, Tmin=2.1, Tmax=2.5, Ts=30, M=1, R=1, new_runs=False, L=[40,60,80,100]):
    """
    Make datafile for each L in temperature analysis.
    Run this func as if it plots the result, while
    in fact is only makes sure the data files are there.
    It then call on plot_temps() to do the actual plotting.
    This is just to separate the tasks-
    It is not necessary to understand how this function works.
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
        time.sleep(2)  # only check for finished programs every 2 seconds
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

    e = pd.DataFrame(e.T, columns=l)
    m = pd.DataFrame(m.T, columns=l)
    Cv = pd.DataFrame(Cv.T, columns=l)
    chi = pd.DataFrame(chi.T, columns=l)

    fig = go.Figure()
    colors = px.colors.qualitative.Plotly

    for variable, title, unit in zip([e, m, Cv, chi], ["Mean energy ", "Mean magnetisation ", "Heat capacity ", "Susceptibility "], ["[J]", "", "", "[1/J]"]):
        c = 0
        fig = go.Figure()

        for size in l:
            name = f"Lattice size: {size} x {size}"
            fig.add_trace(go.Scatter(x = Ts, y = variable[size], mode="lines", line=dict(width = 5, color=colors[c]), name=name))
            c += 1

        Title = title + "as function of temperature for different lattice sizes"

        fig.update_layout(
            font_family="Open sans",
            font_size=30,
            title = Title,
            xaxis_title="Temperature [J]",
            yaxis_title= title + unit,
            legend=dict(yanchor="top", xanchor="right", x=0.99, y=0.99))

        if title == "Mean energy ":
            fig.update_layout(
                legend=dict(yanchor="top", xanchor="left", x=0.01, y=0.99))


        fig.show()


def critical_temp(fname, Tmin = 2.25, Tmax = 2.325, Ts = 15, L = "[40, 60, 80, 100]"):
    L = np.asarray(eval(L))
    T = np.linspace(float(Tmin), float(Tmax), int(Ts))
    m = np.linspace(float(Tmin), float(Tmax), 10001)
    for q, qname, ylab in zip(["Cv", "chi"], ["Heat capacity", "Susceptibility"], ["Heat capacity", "Susceptibility [1/J]"]):
        Tc = np.zeros(len(L))
        Qs = np.zeros(len(L))

        splines = go.Figure()
        c = 0
        colors = px.colors.qualitative.Plotly
        for i, l in enumerate(L):
            file = datapath + fname + f"_{l}.csv"
            data = pd.read_csv(file, header=1, sep=",")
            spline = UnivariateSpline(T, data[q], k=5, s=4)
            Tc[i] = m[np.argmax(spline(m))]
            print(qname, "Size: ", l, "Tc: ", Tc[i])   #print the critical temperatures
            Qs[i] = spline(Tc[i])

            name = f"Size: {l} x {l}"
            splines.add_trace(go.Scatter(x=T, y=data[q], mode="markers", marker=dict(size=10, color=colors[c]), name=name))
            splines.add_trace(go.Scatter(x=m, y=spline(m), mode="lines", line=dict(width=4, color=colors[c]), name="Fitted line"))
            c += 1

        title = f"{qname} with fitted lines for different lattice sizes"
        splines.update_layout(
                font_family="Open sans",
                font_size=30,
                title = title,
                xaxis_title=r"$\LARGE \text{Temperature  } [J]$",
                yaxis_title= ylab,
                legend=dict(yanchor="top", xanchor="right", x=0.99, y=0.99))

        splines.show()

        res = linregress(1 / L, Tc)
        slope = ufloat(res.slope, res.stderr)
        Tinfty = ufloat(res.intercept, res.intercept_stderr)
        print(f"Tc(L=∞) = {Tinfty:.1u}")
        print(f"slope a = {slope:.1u}")



def pdf():

    data_T_low = pd.read_csv("data/pdf_T1.csv", header = 0, sep = ",")
    data_T_high = pd.read_csv("data/pdf_T2.4.csv", header = 0, sep = ",")

    binwidth = 0.01


    fig_low = px.histogram(
        x=data_T_low["e_avg"],
        y=data_T_low["prob"],
        nbins=int(2.4/binwidth))

    fig_low.update_layout(
        xaxis_range=[-2.2, 0.2],   #all p=0 after e=0.2
        #yaxis_range=[0,0.9],
        font_family="Open sans",
        font_size=30,
        title="Histogram of measured probability distribution of the energy for T = 1",
        xaxis_title="Mean energy [J]",
        yaxis_title="Frequency",
        legend=dict(yanchor="top", xanchor="left", x=0.01, y=0.99))


    fig_high = px.histogram(
        x=data_T_high["e_avg"],
        y=data_T_high["prob"],
        nbins=300)

    fig_high.update_layout(
        xaxis_range=[-2.2, 0.2],   #all p=0 after e=0.2
        #yaxis_range=[0, 0.9],
        font_family="Open sans",
        font_size=30,
        title="Histogram of measured probability distribution of the energy for T = 2.4",
        xaxis_title="Mean energy [J]",
        yaxis_title="Frequency",
        legend=dict(yanchor="top", xanchor="left", x=0.01, y=0.99))

    fig_low.show()
    fig_high.show()


def paralympics():
    """
    STOP. MC Hammer time!

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
                run = subprocess.run(f"./{o}.out {M} {R} {L} {T} {p}".split(" "), stdout=subprocess.PIPE)
                times[i] = float(run.stdout.decode().strip())
            print(p, o, np.mean(times))


def main():
    if "help" in sys.argv:
        print("First argument after 'python plot.py' must be a function in this file")
        print("If this function takes any arguments, the are to be passed as well")
        print("Example: 'python3 plot.py burntime burn 2.4 10")

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
