import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd


def main():
    df = pd.read_csv("test.txt", header=0, sep=" ")
    print(df)
    time0 = df.loc[lambda df: df["time"] == 0, :]
    # fig = go.Figure(
    #     data=[go.Scatter(time0, x="x", y="y")]
    # )

    fig = px.scatter_3d(df, x="x", y="y", z="z", animation_frame="time", animation_group="particle",
                        range_x=[-1e4, 1e4], range_y=[-1e4, 1e4], range_z=[-1e4, 1e4]
                        )
    fig.show()


if __name__ == "__main__":
    main()
