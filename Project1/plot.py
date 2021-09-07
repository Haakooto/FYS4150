import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import sys

n = int(sys.argv[1])

def num_sol():
    file = f"data_{n}.txt"
    data = pd.read_csv(file, header=0, sep=" ")

    fig = px.line(data, x="x", y=["u", "v"], title=f"Numerical solution for N=10^{n}")
    fig.update_layout(
        font_family="Garamond",
        font_size=40,
        title="Numeric solution for different values of N",
        xaxis_title="x",
        yaxis_title="v(x)",
    )

    fig.write_image(f"plots/num_sol_n{n}.pdf")
    fig.show()


def first_few_num_sols():
    fig = go.Figure()

    colors = ["seagreen", "aqua", "gold"]
    for i in range(1, n + 1):
        file = f"data_{i}.txt"
        data = pd.read_csv(file, header=0, sep=" ")

        fig.add_trace(go.Scatter(x=data["x"], y=data["v"], mode="lines", line=dict(width=7, color=colors[i - 1]), name=f"N = 10^{i}"))

    fig.add_trace(go.Scatter(x=data["x"], y=data["u"], line=dict(width=7, dash="dot", color="firebrick"), name=f"Analytic solution"))

    fig.update_layout(
        font_family="Garamond",
        font_size=40,
        title="Numerical solution for different values of N",
        xaxis_title="x",
        yaxis_title="v(x)",

    )
    fig.write_image("plots/num_sol_N123.pdf")
    fig.show()


def main():
    if len(sys.argv) > 2:
        method = sys.argv[2]
    else:
        method = "num_sol"
    exec(method + "()")

if __name__ == "__main__":
    main()