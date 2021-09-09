import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression as LR
import sys

n = int(sys.argv[1])
if len(sys.argv) == 4:
    s = "_optim_"
else:
    s = ""
f = lambda type, algo, n: f"datas/{type}{algo}" + (f"_{n}" if type == "full" else "") + ".txt"
data = lambda type, algo, n="full": pd.read_csv(f(type, algo, n), header=0, sep=" ")

def num_sol():
    fig = go.Figure()
    dat = data("full", s, n)
    #fig = px.line(data("full", s, n), x="x", y=["u"], line=dict(width=7))

    fig.add_trace(go.Scatter(x=dat["x"], y=dat["u"], line=dict(width=7, color="firebrick")))

    fig.update_layout(
        font_family="Garamond",
        font_size=35,
        title=f"Analytical solution for N = 10^{n}",
        xaxis_title="x",
        yaxis_title="u(x)",
    )

    #fig.write_image(f"figures/num_sol_n{n}.pdf")
    fig.show()


def first_few_num_sols():
    fig = go.Figure()

    colors = ["seagreen", "aqua", "gold"]
    for i in range(1, n + 1):
        dat = data("full", s, i)

        fig.add_trace(go.Scatter(x=dat["x"], y=dat["v"], mode="lines", line=dict(width=7, color=colors[i - 1]), name=rf"N = 10^{i}"))

    fig.add_trace(go.Scatter(x=dat["x"], y=dat["u"], line=dict(width=7, dash="dot", color="firebrick"), name=f"Analytic solution"))

    fig.update_layout(
        font_family="Garamond",
        font_size=40,
        title="Numerical solution for different values of N",
        xaxis_title="x",
        yaxis_title="v(x)",

    )
    #fig.write_image("figures/num_sol_N123.pdf")
    fig.show()


def absolute_error():
    fig = go.Figure()
    for i in range(1, n + 1):
        dat = data("full", s, i)

        fig.add_trace(go.Scatter(x=dat["x"][1:-1],
                                 y=dat["abs_err"][1:-1],
                                 mode="lines",
                                 line=dict(width=7),
                                 name=f"N = 10^{i}"))
    fig.update_layout(
        font_family="Garamond",
        font_size=40,
        title="Absolute error for different values of N",
        xaxis_title="x",
        yaxis_title="log10(absolute error)",
    )
    fig.show()


def relative_error():
    fig = go.Figure()
    for i in range(1, n + 1):
        dat = data("full", s, i)

        fig.add_trace(go.Scatter(x=dat["x"][1:-1],
                                 y=dat["rel_err"][1:-1],
                                 mode="lines",
                                 line=dict(width=7),
                                 name=f"N = 10^{i}"))
    fig.update_layout(
        font_family="Garamond",
        font_size=40,
        title="Relative error for different values of N",
        xaxis_title="x",
        yaxis_title="log10(relative error)",
    )
    fig.show()


def max_error():
    fig = go.Figure()
    points_to_exclude = 2
    x = np.asarray(data("limited", s, "This agrument is useless")["n"])
    y = np.asarray(data("limited", s, "This agrument is useless")["max_error"])

    model = LR()
    model.fit(x[:-points_to_exclude].reshape(-1, 1), y[:-points_to_exclude].reshape(-1, 1))
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
        font_size=35,
        title="Max relative error for " + ("optimized " if s == "_optim_" else "general ") + "Thomas algorithm as function of N",
        xaxis_title="Log10(N)",
        yaxis_title="Max relative error",
        legend=dict(yanchor="bottom", xanchor="left", x=0.01, y=0.01),
    )
    fig.show()


def both_max_error():
    fig = go.Figure()
    points_to_exclude = (0, 2)
    x = np.asarray(data("limited", "", "np.log10(this is a usless argument)")["n"])
    y_norma = np.asarray(data("limited", "", "np.log10(this is a usless argument)")["max_error"])
    y_optim = np.asarray(data("limited", "_optim_", "np.log10(this is a usless argument)")["max_error"])

    model = LR()
    model.fit(x[:-points_to_exclude[1]].reshape(-1, 1), y_norma[:-points_to_exclude[1]].reshape(-1, 1))
    y_predict_norma = model.predict(x.reshape(-1, 1))

    fig.add_trace(go.Scatter(x=x, y=y_norma,
                             mode="markers",
                             marker=dict(size=20,
                                         color="maroon"),
                             name="General Thomas algorithm"
                             )
                  )

    fig.add_trace(go.Scatter(x=x, y=y_predict_norma.T[0],
                             mode="lines",
                             line=dict(dash="dot",
                                       width=7,
                                       color="firebrick"),
                             name=f"Fitted line. Slope = {model.coef_[0][0]:.5}"
                             )
                  )

    model = LR()
    model.fit(x.reshape(-1, 1), y_optim.reshape(-1, 1))
    y_predict_optim = model.predict(x.reshape(-1, 1))

    fig.add_trace(go.Scatter(x=x, y=y_optim,
                             mode="markers",
                             marker=dict(size=20,
                                         color="darkgreen"),
                             name="Optimized Thomas algorithm"
                             )
                  )

    fig.add_trace(go.Scatter(x=x, y=y_predict_optim.T[0],
                             mode="lines",
                             line=dict(dash="dot",
                                       width=7,
                                       color="mediumseagreen"),
                             name=f"Fitted line. Slope = {model.coef_[0][0]:.5}"
                             )
                  )

    fig.update_layout(
        font_family="Garamond",
        font_size=40,
        title="Max relative error for algoritms as function of N",
        xaxis_title="log10(N)",
        yaxis_title="log10(Max relative error)",
        legend=dict(yanchor="bottom", xanchor="left", x=0.01, y=1-0.99),
    )
    fig.show()

def timing():
    fig = go.Figure()

    x = np.arange(1, n + 1, 1)
    x_raw = np.asarray(data("limited", "", "np.log10(this is a usless argument)")["n"])
    y_g_raw = np.log10(np.asarray(data("limited", "", "np.log10(this is a usless argument)")["time"]))
    y_o_raw = np.log10(np.asarray(data("limited", "_optim_", "np.log10(this is a usless argument)")["time"]))
    y_g_mean = np.zeros(n)
    y_g_std = np.zeros(n)
    y_o_mean = np.zeros(n)
    y_o_std = np.zeros(n)

    for i in range(n):
        y_g = y_g_raw[np.where(x_raw == i + 1)]
        y_g_mean[i] = np.mean(y_g)
        y_g_std[i] = np.std(y_g)
        y_o = y_o_raw[np.where(x_raw == i + 1)]
        y_o_mean[i] = np.mean(y_o)
        y_o_std[i] = np.std(y_o)

    model = LR()
    model.fit(x.reshape(-1, 1), y_g_mean.reshape(-1, 1))
    y_predict_norma = model.predict(x.reshape(-1, 1))

    fig.add_trace(go.Scatter(x=x, y=y_g_mean,
                             error_y=dict(type='data', array=y_g_std, visible=True),
                             mode="markers",
                             marker=dict(size=6,
                                         color="maroon"),
                             name="General Thomas algorithm"
                             )
                  )

    fig.add_trace(go.Scatter(x=x, y=y_predict_norma.T[0],
                             mode="lines",
                             line=dict(dash="dot",
                                       width=3,
                                       color="firebrick"),
                             name=f"Fitted line. Slope = {model.coef_[0][0]:.5}"
                             )
                  )

    model = LR()
    model.fit(x.reshape(-1, 1), y_o_mean.reshape(-1, 1))
    y_predict_optim = model.predict(x.reshape(-1, 1))

    fig.add_trace(go.Scatter(x=x, y=y_o_mean,
                             error_y=dict(type='data', array=y_o_std, visible=True),
                             mode="markers",
                             marker=dict(size=6,
                                         color="darkgreen"),
                             name="Optimized Thomas algorithm"
                             )
                  )

    fig.add_trace(go.Scatter(x=x, y=y_predict_optim.T[0],
                             mode="lines",
                             line=dict(dash="dot",
                                       width=3,
                                       color="mediumseagreen"),
                             name=f"Fitted line. Slope = {model.coef_[0][0]:.5}"
                             )                  )

    fig.update_layout(
        font_family="Garamond",
        font_size=35,
        title="Running time of algorithm as function of N",
        xaxis_title="Log10(N)",
        yaxis_title="Log10(Time)",
        legend=dict(yanchor="top", xanchor="left", x=0.01, y=0.99),
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
