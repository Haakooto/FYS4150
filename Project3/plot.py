import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import glob
from PIL import Image


img = Image.open("giffel.jpg")
giffel = lambda x, y, s=1: dict(source=img, x=int(x), y=int(y), sizex=s, sizey=s, xref="x", yref="y", xanchor="center", yanchor="middle", layer="above", opacity=1)

def plot_z():  # Ex9p1
    df = pd.read_csv("outputs/oneP_endurance.txt", header=0, sep=" ")

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=df["time"], y=df["z"], mode="lines"))

    N = df.shape[0]

    dt = df["time"][1]
    yf = np.fft.rfft(df["z"])
    xf = np.fft.rfftfreq(N, d=dt)

    wz = np.sqrt(2 * 9.65 / 40.078)
    freq = 2 * np.pi * xf[np.argmax(np.abs(yf))]
    print(freq)
    print(1 - freq / wz)  # relative error in freq

    anal = 0.5 * (np.exp(1j * wz * df["time"]) + np.exp(-1j * wz * df["time"]))

    trace = go.Scatter(x=df["time"], y=df["z"] / df["z"][0], mode="lines", line=dict(width=5))
    trace2 = go.Scatter(x=df["time"], y=np.real(anal), mode="lines", line=dict(dash="dot"))
    fig = go.Figure(data=[trace, trace2])

    fig.show()


def plot_xy_plane():
    max_time = 20
    fig = go.Figure()
    for b in ["", "_no"]:
        df = pd.read_csv(f"outputs/twoP{b}_ppi.txt", header=0, sep=" ")
        df = df.loc[lambda x: x["time"] < max_time, :]
        for i in range(1, 3):
            dP = df.loc[lambda x: x["particle"] == i, :]
            fig.add_trace(go.Scatter(x=dP["x"], y=dP["y"], mode="lines", line=dict(width=4), name=f"particle {i} {b} ppi"))
            start = dP.loc[lambda x: x["time"] == 0, :]
            fig.add_layout_image(giffel(start["x"], start["y"]))

    fig.update_layout(title="Position of 2 particles in xy-plane, with and without interaction",
                      xaxis_title="$x [\mu m]",
                      yaxis_title="$y [\mu m]",
                      font_family="Open sans",
                      font_size=45,
                      )
    fig.show()

def plot_phase_diagrams():
    max_time = 10
    for x in "xyz":
        fig = go.Figure()
        for b in ["", "_no"]:
            df = pd.read_csv(f"outputs/twoP{b}_ppi.txt", header=0, sep=" ")
            df = df.loc[lambda x: x["time"] < max_time, :]
            for i in range(1, 3):
                dP = df.loc[lambda x: x["particle"] == i, :]
                fig.add_trace(go.Scatter(x=dP[x], y=dP["v" + x], mode="lines", line=dict(width=4), name=f"particle {i} {b} ppi"))
                start = dP.loc[lambda x: x["time"] == 0, :]
                fig.add_layout_image(giffel(start[x], start["v" + x], 0.5 if x == "z" else 0.7))
        fig.update_layout(title=f"Phase space plot for {x}, with and without interaction",
                    xaxis_title=f"{x} [\mu m]",
                    yaxis_title=f"v{x} [\mu m / \mu s]",
                    font_family="Open sans",
                    font_size=45,
                    )
        fig.show()



def plot_rel_errors(method="RK4"):  # Ex9p5
    files = sorted(glob.glob(f"outputs/rel_errors_{method}*"))
    traces = []
    for file in files:
        data = pd.read_csv(file, header=0, sep=" ")
        data["err"][0] = 0
        s = file.rfind("_")
        st = file.find(".")
        dt = file[s + 1: st]

        trace = go.Scatter(x=data["t"], y=np.log10(data["err"]), mode="lines", line=dict(width=5), name=f"log10(h) = - {dt}")
        traces.append(trace)
    fig = go.Figure(data=traces)
    fig.update_layout(title=f"Relative error as function of time using {method} for different timesteps",
                      xaxis_title=r"$\Huge \text{Time} [\mu s]$",
                      yaxis_title="log10(relative error)",
                      font_size=45,
                      font_family="Open sans",
                      )
    fig.show()

def error_conv_rate(method="RK4"):  # Ex9p6
    files = sorted(glob.glob(f"outputs/rel_errors_{method}*"))  # load files
    datas = [pd.read_csv(file, header=0, sep=" ") for file in files]  # read files
    errs = np.asarray([data["err"][0] for data in datas])  # extract max abs_err
    dts = np.asarray([data["t"][1] for data in datas])  # extract dt
    derrs = np.log(errs[1:] / errs[:-1])  # find log of ratios
    ddts = np.log(dts[1:] / dts[:-1])  # find log of ratios
    conv = sum(derrs / ddts) / (len(files) - 1)  # do sum
    print(f"Error convergence rate for {method} is {conv}")  # print



    N = df.shape[0]

    dt = df["time"][1]
    yf = np.fft.rfft(df["z"])
    xf = np.fft.rfftfreq(N, d=dt)

    data = pd.read_csv("outputs/ex10_TimeTrap_particle_track_f0.4_w0.44.txt", header = 0, sep = " ")

    fig.add_trace(go.Scatter(x=data["time"], y=data["z"], mode="lines"))

    fig.update_layout(
    xaxis_range=[0, 200],
    yaxis_range=[-500,500],
    font_family="Open sans",
    font_size=45,
    xaxis_title=r"$\Huge \text{Time} [\mu s] $",
    yaxis_title=r"$\Huge \text{z} [\mu m] $",
    title=r"$\Huge{\text{Position  along  z-axis  with  } \textit{f} = 0.4, \omega_V = 0.44}$")
    fig.show()


def ex_10_plot_both_tracks_z():
    f=0.4
    w=0.44
    fig = go.Figure()
    dfTime = pd.read_csv(f"outputs/ex10_TimeTrap_particle_track_f{f}_w{w}.txt", header = 0, sep = " ")
    dfRegular = pd.read_csv(f"outputs/ex10_RegularTrap_particle_track_f{f}_w{w}.txt", header = 0, sep = " ")

    fig.add_trace(go.Scatter(
        x=dfTime["time"],
        y=dfTime["z"],
        mode="lines",
        line=dict(width=4),
        name = "Time-dependent potential"))

    fig.add_trace(go.Scatter(
        x=dfRegular["time"],
        y=dfRegular["z"],
        mode="lines",
        line=dict(width=4),
        name = "Constant potential"))


    fig.update_layout(
    xaxis_range=[0, 150],
    yaxis_range=[-600,600],
    font_family="Open sans",
    font_size=45,
    title=r"$\Huge{\text{Position  along  z-axis  with  } \textit{ f = 0.4}, \omega_V \textit{= 0.44}}$",
    xaxis_title=r"$\Huge \text{Time  } [\mu s] $",
    yaxis_title=r"$\Huge \text{z  } [\mu m] $",
    legend=dict(yanchor="top", xanchor="left", x=0.01, y=0.99))
    fig.show()

def plot_freqs_z():
    f = 0.4
    w = 0.44
    fig = go.Figure()
    dfTime = pd.read_csv(f"outputs/ex10_TimeTrap_particle_track_f{f}_w{w}.txt", header=0, sep=" ")
    dfRegular = pd.read_csv(f"outputs/ex10_TimeRegular_particle_track_f{f}_w{w}.txt", header=0, sep=" ")

    for df in (dfTime, dfRegular):
        N = df.shape[0] 

        dt = df["time"][1]
        yf = np.fft.rfft(df["z"])
        xf = np.fft.rfftfreq(N, d=dt)

        fig.add_trace(go.Scatter(x=xf, y=np.abs(yf), mode="lines"))

    fig.show()

def ex_10_track_xy():
    fig = go.Figure()
    f = 0.4
    w = 0.44

    dfTime = pd.read_csv(f"outputs/ex10_TimeTrap_particle_track_f{f}_w{w}.txt", header = 0, sep = " ")
    dfRegular = pd.read_csv(f"outputs/ex10_RegularTrap_particle_track_f{f}_w{w}.txt", header = 0, sep = " ")

    dfTime = dfTime.loc[lambda df: df["time"] < 130, :]
    dfRegular = dfRegular.loc[lambda df: df["time"] < 130, :]

    fig.add_trace(go.Scatter(
        x=dfTime["x"],
        y=dfTime["y"],
        mode="lines",
        line=dict(width=3),
        name = "Time-dependent potential"))

    fig.add_trace(go.Scatter(
        x=dfRegular["x"],
        y=dfRegular["y"],
        mode="lines",
        line=dict(width=3),
        name = "Constant potential"))


    fig.update_layout(
    xaxis_range=[-200,300],
    yaxis_range=[-300,200],
    font_family="Open sans",
    font_size=45,
    title = r"$\Huge{\text{Position  in the xy-plane  with  } \textit{ f = 0.4}, \omega_V \textit{= 0.44}}$",
    xaxis_title=r"$\Huge \text{x  } [\mu m] $",
    yaxis_title=r"$ \Huge \text{y  } [\mu m] $",
    legend=dict(yanchor="top", xanchor="left", x=0.01, y=0.99))
    fig.show()




def ex10_broad_plot_fraction_remaining():
    fig = go.Figure()
    file = pd.read_csv("outputs/broad_freq_search_new.txt", header = 2, sep = " " )
    print(file)
    #timestep er 0.0025
    for f in [0.1, 0.4, 0.7]:

        data = file.loc[lambda x: x["ampl"] == f, :]
        fig.add_trace(go.Scatter(x = data["wV"],
            y = data["fracRem"],
            mode = "lines",
            line=dict(width=4),
            name = f"Amplitude: {f}"))

        fig.update_layout(
            font_family="Open sans",
            font_size=45,
            title = r"$\Huge{\text{Fraction  of  particles  remaining  in  the  Penning  trap  after  0.5  } \mu s }$",
            xaxis_title=r"$\Huge\omega_V$",
            yaxis_title="Fraction",
            legend=dict(yanchor="bottom", xanchor="right", x=0.99, y=0.01)
            )
    fig.show()


def ex10_narrow_plot_no_ppi_fraction_remaining():
    fig = go.Figure()
    file = pd.read_csv("outputs/narrow_freq_search_no_ppi.txt", header = 0, sep = " " )
    #timestep er 0.0025
    for f in [0.1, 0.4, 0.7]:
        data = file.loc[lambda x: x["ampl"] == f, :]
        fig.add_trace(go.Scatter(x = data["wV"],
        y = data["fracRem"],
        mode = "lines",
        line=dict(width=4),
        name = f"Amplitude: {f}"))

    fig.update_layout(
        font_family="Open sans",
        font_size=45,
        title=r"$\uge{\text{Fraction  of  particles  remaining  in  the  Penning  trap  after  0.5  }  \mu s}$",
        xaxis_title=r"$\Huge\omega_V$",
        yaxis_title="Fraction",
        legend=dict(yanchor="top", xanchor="left", x=0.01, y=0.99)
        )
    fig.show()


if __name__ == "__main__":
    # plot_z()
    # plot_rel_errors()
    # plot_rel_errors("Euler")
    #error_conv_rate()
    #error_conv_rate("Euler")
    # plot_xy_plane()
    # plot_phase_diagrams()
    # ex10_broad_plot_fraction_remaining()
    plot_freqs_z()
    # ex_10_plot_both_tracks_z()
    # ex_10_track_xy()
    #ex10_narrow_plot_no_ppi_fraction_remaining()
