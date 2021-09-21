import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
import sys


def iterations(N):
    Iterations = []
    Dim = []

    for i in range (5, N+1):
        data = pd.read_csv("results/data_" + str(i) + ".txt", header=0, sep=" ")

        Iterations.append(data.iloc[0, 0])
        Dim.append(i)


    poly = PolynomialFeatures(degree=2)
    Dim_transform = poly.fit_transform(np.reshape(Dim, (-1, 1)))

    model = LinearRegression()
    model.fit(Dim_transform, np.reshape(Iterations, (-1, 1)))

    predict = model.predict(Dim_transform)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=Dim, y=Iterations, mode="markers", marker=dict(size=9, color="firebrick"), name="Data"))
    fig.add_trace(go.Scatter(x=Dim, y=predict.ravel(),
        line=dict(width=4,
        color="darkgreen"),
        name=f"Fitted line. Iterations = {model.coef_[0][0]:.5} {model.coef_[0][1]:.3}N + {model.coef_[0][2]:.3}N^2" ))

    fig.update_layout(
        font_family="Garamond",
        font_size=35,
        title="Number of iterations as a function of N",
        xaxis_title="N",
        yaxis_title="Iterations",
        legend=dict(yanchor="top", xanchor="left", x=0.01, y=0.99)
        )

    fig.show()


N = 100

iterations(N)
