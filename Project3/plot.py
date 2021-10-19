import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd


def main():
    df = pd.read_csv("outputs/oneP_endurace.txt", header=0, sep=" ")
    print(df)
    time0 = df.loc[lambda df: df["time"] == 0, :]
    # fig = go.Figure(
    #     data=[go.Scatter(time0, x="x", y="y")]
    # )

    fig = px.scatter_3d(df, x="x", y="y", z="z", animation_frame="time", animation_group="particle",
                        range_x=[-1e4, 1e4], range_y=[-1e4, 1e4], range_z=[-1e4, 1e4]
                        )
    fig.show()

def plot_z():
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
    fig = go.Figure()
    for b in ["", "_no"]:
        df = pd.read_csv(f"outputs/twoP{b}_ppi.txt", header=0, sep=" ")
        for i in range(1, 2):
            dP = df.loc[lambda x: x["particle"] == i, :]
            fig.add_trace(go.Scatter(x=dP["x"], y=dP["y"], mode="lines", name=f"particle {i} {b} ppi"))

    fig.update_layout(xaxis_range=[-30, 30], yaxis_range=[-30,30])
    # fig = px.scatter(df, x="x", y="y", animation_frame="time", animation_group="particle", range_y=[-20, 20])
    fig.show()



def ex_10_plot_xy_plane():
    fig = go.Figure()
    data = pd.read_csv("outputs/ex10_particle_track.txt", header = 0, sep = " ")

    fig.add_trace(go.Scatter(x=data["x"], y=data["y"], mode="lines"))

    fig.update_layout(
    xaxis_range=[-500, 500],
    yaxis_range=[-500,500],
    font_family="Garamond",
    font_size=30,
    title="",
    xaxis_title="x",
    yaxis_title="y",)
    fig.show()



def ex10_broad_plot_fraction_remaining():
    fig = go.Figure()
    file = pd.read_csv("outputs/broad_freq_search_test.txt", header = 0, sep = " " )
    #timestep er 0.0025
    for f in [0.1, 0.4, 0.7]:

        data = file.loc[lambda x: x["ampl"] == f, :]
        fig.add_trace(go.Scatter(x = data["wV"],
            y = data["fracRem"],
            mode = "lines",
            line=dict(width=3),
            name = f"Amplitude: {f}"))

        fig.update_layout(
            font_family="Garamond",
            font_size=30,
            title="Fraction of particles remaining in the Penning trap after 0.5 milliseconds testing",
            xaxis_title="Frequency",
            yaxis_title="Fraction of particles remaining",
            legend=dict(yanchor="bottom", xanchor="left", x=0.01, y=0.01)
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
        line=dict(width=3),
        name = f"Amplitude: {f}"))

    fig.update_layout(
        font_family="Garamond",
        font_size=30,
        title="Fraction of particles remaining in the Penning trap after 0.5 milliseconds",
        xaxis_title="Frequency",
        yaxis_title="Fraction of particles remaining",
        legend=dict(yanchor="top", xanchor="left", x=0.01, y=0.99)
        )
    fig.show()


if __name__ == "__main__":
    # main()
    plot_z()
    #plot_xy_plane()
    #ex10_broad_plot_fraction_remaining()
    #ex_10_plot_xy_plane()
    # ex10_broad_plot_fraction_remaining()
    #ex10_narrow_plot_no_ppi_fraction_remaining()
