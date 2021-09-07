import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression as LR
import sys

n = int(sys.argv[1])

f = lambda x: f"datas/data_{x}.txt"


def num_sol():
    data = pd.read_csv(f(n), header=0, sep=" ")

    fig = px.line(data, x="x", y=["u", "v"])
    fig.update_layout(
        font_family="Garamond",
        font_size=30,
        title=rf"Numerical solution of N = 10 \U+207{n}",
        xaxis_title="x",
        yaxis_title="v(x)",
    )

    fig.write_image(f"figures/num_sol_n{n}.pdf")
    fig.show()


def first_few_num_sols():
    fig = go.Figure()

    colors = ["seagreen", "aqua", "gold"]
    for i in range(1, n + 1):
        data = pd.read_csv(f(i), header=0, sep=" ")

        fig.add_trace(go.Scatter(x=data["x"], y=data["v"], mode="lines", line=dict(width=7, color=colors[i - 1]), name=rf"N = 10\U+207{i}"))

    fig.add_trace(go.Scatter(x=data["x"], y=data["u"], line=dict(width=7, dash="dot", color="firebrick"), name=f"Analytic solution"))

    fig.update_layout(
        font_family="Garamond",
        font_size=40,
        title="Numerical solution for different values of N",
        xaxis_title="x",
        yaxis_title="v(x)",

    )
    fig.write_image("figures/num_sol_N123.pdf")
    fig.show()


def absolute_error():
    fig = go.Figure()
    for i in range(1, n + 1):
        data = pd.read_csv(f(i), header=0, sep=" ")

        fig.add_trace(go.Scatter(x=data["x"][1:-1],
                                 y=data["abs_err"][1:-1],
                                 mode="lines",
                                 line=dict(width=7),
                                 name=f"N = 10^{i}"))
    fig.update_layout(
        font_family="Garamond",
        font_size=40,
        title="Absolute error for different values of N",
        xaxis_title="x",
        yaxis_title="Absolute error",
    )
    fig.show()


def relative_error():
    fig = go.Figure()
    for i in range(1, n + 1):
        data = pd.read_csv(f(i), header=0, sep=" ")

        fig.add_trace(go.Scatter(x=data["x"][1:-1],
                                 y=data["rel_err"][1:-1],
                                 mode="lines",
                                 line=dict(width=7),
                                 name=f"N = 10^{i}"))
    fig.update_layout(
        font_family="Garamond",
        font_size=40,
        title="Relative error for different values of N",
        xaxis_title="x",
        yaxis_title="Relative error",
    )
    fig.show()


def max_error():
    fig = go.Figure()
    x = np.arange(n) + 1
    y = np.zeros(n)
    for i in range(n):
        y[i] = np.min(pd.read_csv(f(i + 1), header=0, sep=" ")["rel_err"][1:])

    model = LR()
    model.fit(x[np.nonzero(y)].reshape(-1, 1), y[np.nonzero(y)].reshape(-1, 1))
    y_predict = model.predict(x.reshape(-1, 1))

    fig.add_trace(go.Scatter(x=x, y=y,
                             mode="markers",
                             marker=dict(size=20,
                                         color="maroon"),
                             name="Errors"
                             )
                  )

    fig.add_trace(go.Scatter(x=x, y=y_predict.T[0],
                             mode="lines",
                             line=dict(dash="dot",
                                       width=7,
                                       color="firebrick"),
                             name=f"Fitted line. Slope = {model.coef_[0][0]:.5}"
                             )
                  )
    fig.update_layout(
        font_family="Garamond",
        font_size=40,
        title="Max relative error as function of N",
        xaxis_title="Log10(N)",
        yaxis_title="Max relative error",
        legend=dict(yanchor="bottom", xanchor="left", x=0.01, y=0.01),
    )
    fig.show()


def main():
    if len(sys.argv) > 2:
        method = sys.argv[2]
    else:
        method = "num_sol"
    exec(method + "()")


if __name__ == "__main__":
    main()