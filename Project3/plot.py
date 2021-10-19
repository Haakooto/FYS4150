import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd


def make_nice_plot(fig, xlabel="", ylabel="", title="", legend=None):
    layout = go.Layout(
        font_family="Garamond",
        font_size=35,
        xaxis=dict(title=xlabel),
        yaxis=dict(title=ylabel),
        title=title,
    )
    # fig.update_layout(font_family="Garamond",
    #                   font_size=35,
    #                   title=title,
    #                   xaxis_title=xlabel,
    #                   yaxis_title=ylabel,
    #                   )
    # if legend is not None:
        # fig.update_layout(legend=legend)
    return layout

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


def ex9_plot_z():
    df = pd.read_csv("outputs/oneP_endurance.txt", header=0, sep=" ")
    N = df.shape[0]

    dt = df["time"][1]
    yf = np.fft.rfft(df["z"])
    xf = np.fft.rfftfreq(N, d=dt)

    wz = np.sqrt(2 * 9.65 / 40.078)
    freq = 2 * np.pi * xf[np.argmax(np.abs(yf))]
    print(1 - freq / wz)  # relative error in freq

    trace = go.Scatter(x=df["time"], y=df["z"], mode="lines")
    lo = make_nice_plot(r"$Time [\mu s]$", r"$z [\mu m]$",
                   r"$\text{Movement in z-direction for single particle}$")
    fig = go.Figure(data=[trace,], layout=lo)
    fig.show()


def plot_xy_plane():
    fig = go.Figure()
    for b in ["", "_no"]:
        df = pd.read_csv(f"outputs/twoP{b}_ppi.txt", header=0, sep=" ")
        for i in range(1, 2):
            dP = df.loc[lambda x: x["particle"] == i, :]
            fig.add_trace(go.Scatter(
                x=dP["x"], y=dP["y"], mode="lines", name=f"particle {i} {b} ppi"))

    fig.update_layout(xaxis_range=[-30, 30], yaxis_range=[-30, 30])
    # fig = px.scatter(df, x="x", y="y", animation_frame="time", animation_group="particle", range_y=[-20, 20])
    fig.show()


if __name__ == "__main__":
    # main()
    ex9_plot_z()
