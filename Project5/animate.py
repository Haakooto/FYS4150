import numpy as np
import plotly.graph_objects as go

def main():
    data = np.load("test.npz")
    m = int(np.sqrt(data.shape[1] - 1))
    t = data[:, 0]
    Z = data[:, 1:]
    x = np.linspace(0, 1, m)

    data = go.Contour(x=x, y=x,
                      z=Z[0].reshape(m,m),
                    #   contours=dict(start=0, end=0.1, size=0.01, showlines=False),
                      )

    frames = [go.Frame(data=[go.Contour(x=x, y=x, z=Z[i].reshape(m,m))]) for i in range(1, len(t))]

    layout = go.Layout(updatemenus=[dict(type="buttons", buttons=[dict(label="Play", method="animate", args=[None])])])

    f = go.Figure(data=data, layout=layout, frames=frames)
    f.show()


if __name__ == "__main__":
    main()
