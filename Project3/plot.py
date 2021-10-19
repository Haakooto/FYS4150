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

if __name__ == "__main__":
    # main()
    plot_xy_plane()