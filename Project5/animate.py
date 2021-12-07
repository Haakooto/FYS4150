import numpy as np
import plotly.graph_objects as go


def animate():
    data = np.load("npz/prob.npz")
    m = int(np.sqrt(data.shape[1] - 1))
    t = data[:, 0]
    Z = data[:, 2:].reshape(len(t), m, m)
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
    p = prob[:, 2:].reshape(len(t), m, m)
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



def p7_probs():
    data = np.load("npz/p7_no_slit.npz")
    m = int(np.sqrt(data.shape[1] - 1))
    t = data[:, 0]
    probs = (data[:, 1])

    fig = go.Figure()

    fig.add_trace(go.Scatter(x = t, y = probs, mode="lines", line=dict(width = 5)))

    fig.update_layout(
        font_family="Garamond",
        font_size=30,
        title = "Deviance of total probability from 1",
        xaxis_title= "Time",
        yaxis_title= "Total probability deviance")

    fig.show()



def p8_three_steps():
    data = np.load("npz/p8.npz")
    m = int(np.sqrt(data.shape[1] - 1))
    t = data[:, 0]
    p = data[:, 2:].reshape(len(t), m, m)

    real = np.load("npz/p8_real.npz")
    imag = np.load("npz/p8_imag.npz")

    axes = np.linspace(0, 1, m)

    times = [0.0, 0.001, 0.002]

    for time in times:
        fig = go.Figure(data=go.Contour(x=axes, y=axes, z=p[np.argmin(abs(t - time))]))
        fig.update_layout(
            font_family="Garamond",
            font_size=30,
            title = f"Probability distribution at time t = {time}",
            xaxis_title= "x",
            yaxis_title= "y")

        fig.show()


    for time in times:
        for name, part in zip(["Imaginary", "Real"], [imag, real]):
            fig = go.Figure(data=go.Contour(x=axes, y=axes, z=part[np.argmin(abs(t - time))]))
            fig.update_layout(
                font_family="Garamond",
                font_size=30,
                title = f"{name} part of probability distribution at time t = {time}",
                xaxis_title= "x",
                yaxis_title= "y")

            fig.show()




if __name__ == "__main__":
    #p7_no_slit()
    p7_probs()
    #animate()
    # plot_states()
    # plot_initial_u()
