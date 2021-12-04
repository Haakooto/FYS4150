import numpy as np
import plotly.graph_objects as go


def animate():
    data = np.load("npz/prob.npz")
    m = int(np.sqrt(data.shape[1] - 1))
    t = data[:, 0]
    Z = data[:, 1:].reshape(len(t), m, m)
    x = np.linspace(0, 1, m)

    data = go.Contour(x=x, y=x,
                      z=Z[0],
                      #   contours=dict(start=0, end=0.1, size=0.01, showlines=False),
                      )

    frames = [go.Frame(data=[go.Contour(x=x, y=x, z=Z[i].reshape(m, m))])
              for i in range(1, len(t))]

    layout = go.Layout(updatemenus=[dict(type="buttons", buttons=[
                       dict(label="Play", method="animate", args=[None])])])

    f = go.Figure(data=data, layout=layout, frames=frames)
    f.show()


def plot_states():
    prob = np.load("npz/prob.npz")
    m = int(np.sqrt(prob.shape[1] - 1))
    t = prob[:, 0]
    p = prob[:, 1:].reshape(len(t), m, m)
    x = np.linspace(0, 1, m)
    fig = go.Figure(data=go.Contour(x=x, y=x, z=p[np.argmin(abs(t - 0))]))
    fig.show()
    fig = go.Figure(data=go.Contour(x=x, y=x, z=p[np.argmin(abs(t - 0.004))]))
    fig.show()
    fig = go.Figure(data=go.Contour(x=x, y=x, z=p[np.argmin(abs(t - 0.008))]))
    fig.show()


def plot_initial_u():
    data = np.load("initial_u.npz")
    data = np.conj(data) * data
    print(data)

    fig = go.Figure(data=go.Contour(z=np.real(data)))
    fig.show()


if __name__ == "__main__":
    animate()
    # plot_states()
    # plot_initial_u()
