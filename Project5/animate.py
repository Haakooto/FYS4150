import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio


def position_to_index(position, h=0.005):
    """ Convert position to index """
    index = position/h - 1
    return int(index)


def time_to_index(time, step=2.5e-5):
    """ Convert time to index """
    index = time/step
    return int(index)


def animate(file="p9_3_slit.npz", h=0.005):
    """ Animate time evolution of a simulation """
    data = np.load("npz/" + file)
    m = int(np.sqrt(data.shape[1] - 1))
    t = data[:, 0]
    Z = data[:, 2:].reshape(len(t), m, m)
    axes = np.linspace(h, 1-h, m)

    data = go.Contour(x=axes, y=axes,
                      z=Z[0],
                      )

    frames = [go.Frame(data=[go.Contour(x=axes, y=axes, z=Z[i].reshape(m, m))])
              for i in range(1, len(t))]

    layout = go.Layout(updatemenus=[dict(type="buttons", buttons=[
                       dict(label="Play", method="animate", args=[None])])])

    fig = go.Figure(data=data, layout=layout, frames=frames)
    fig.show()


def p7_no_slit():
    """ Plot deviation from 1 of total probability as function of time """
    data = np.load("npz/p7_no_slit.npz")
    m = int(np.sqrt(data.shape[1] - 1))
    t = data[:, 0]
    probs = (data[:, 1])

    fig = go.Figure()

    fig.add_trace(go.Scatter(x=t, y=probs, mode="lines", line=dict(width=5)))

    fig.update_layout(
        font_family="Garamond",
        font_size=30,
        title="Deviance of total probability from 1 without slit",
        xaxis_title="Time",
        yaxis_title="Total probability deviance")

    fig.show()


def p7_double_slit():
    """ Plot deviation from 1 of total probability as function of time """
    print("Somehow I dont exist and must be copied from p7_no_slit()")


def animate_real_imag():
    """ Animate the time evolution of real and imag part of the wave function """
    data = np.load("npz/p8.npz")
    m = int(np.sqrt(data.shape[1] - 1))
    t = data[:, 0]

    real = np.load("npz/p8_real.npz").reshape(len(t), m, m)
    imag = np.load("npz/p8_imag.npz").reshape(len(t), m, m)

    axes = np.linspace(0, 1, m + 2)[1:-1]

    real_ = go.Contour(x=axes, y=axes, z=real[0])
    imag_ = go.Contour(x=axes, y=axes, z=imag[0])

    real_frames = [go.Frame(data=[go.Contour(x=axes, y=axes, z=real[i])])
                   for i in range(1, len(t))]
    imag_frames = [go.Frame(data=[go.Contour(x=axes, y=axes, z=imag[i])])
                   for i in range(1, len(t))]

    layout = go.Layout(updatemenus=[dict(type="buttons", buttons=[
                       dict(label="Play", method="animate", args=[None])])])

    rf = go.Figure(data=real_, layout=layout, frames=real_frames)
    imf = go.Figure(data=imag_, layout=layout, frames=imag_frames)
    rf.show()
    imf.show()


def p8_three_steps(h=0.005):
    """ Plot probability function as well as real and imag part of Ψ at 3 time points """
    data = np.load("npz/p8.npz")
    m = int(np.sqrt(data.shape[1] - 1))
    t = data[:, 0]
    p = data[:, 2:].reshape(len(t), m, m)

    real = np.load("npz/p8_real.npz").reshape(len(t), m, m)
    imag = np.load("npz/p8_imag.npz").reshape(len(t), m, m)

    axes = np.linspace(0, 1, m + 2)[1:-1]
    times = [0.0, 0.001, 0.002]

    for k, time in enumerate(times):
        if time == 0.001:
            x1 = np.argmin(abs(axes - 0.2))
            y1 = np.argmin(abs(axes - 0.15))
            x2 = np.argmin(abs(axes - 0.8))
            y2 = np.argmin(abs(axes - 0.85))
        else:
            x1, y1, x2, y2 = 0, 0, -1, -1

        fig = make_subplots(rows=3, cols=1,
                            column_widths=[1],
                            row_heights=[0.33, 0.33, 0.33],
                            specs=[[{"type": "contour"}] for _ in range(3)],
                            subplot_titles=[f"Probability distribution at time t = {time}.",
                                            "Real part of quantum state.", "Imaginary part of quantum state."],
                            horizontal_spacing=0.1,
                            vertical_spacing=0.1,
                            )
        fig.add_trace(go.Contour(x=axes[x1:x2], y=axes[y1:y2],
                      z=p[np.argmin(abs(time - t))][y1:y2, x1:x2], colorbar=dict(len=0.3, x=1, y=0.87)), row=1, col=1)
        fig.add_trace(go.Contour(x=axes, y=axes,
                      z=real[np.argmin(abs(time - t))], colorbar=dict(len=0.3, x=1, y=0.5)), row=2, col=1)
        fig.add_trace(go.Contour(x=axes, y=axes,
                      z=imag[np.argmin(abs(time - t))], colorbar=dict(len=0.3, x=1, y=0.13)), row=3, col=1)
        fig.for_each_annotation(lambda a: a.update(
            font_size=80, font_family="Garamond"))

        for a in ["", "2", "3"]:
            for xy in ["x", "y"]:
                fig["layout"][f"{xy}axis{a}"]["title"] = {"text": xy}

        fig.update_layout(
            font_family="Garamond",
            font_size=70,
            coloraxis=dict(colorbar=dict(len=0.5, x=1, y=0.22)),
            width=2000,
            height=3000,
        )
        pio.write_image(fig, f"p8_{k}.pdf")
        # fig.show()


def p9(h=0.005):
    """ Plot probability distribution at x=0.8 after 2ms in the case of 1, 2, 3 slits """
    nice = {1: "one slit", 2: "two slits", 3: "three slits"}
    for slits in [1, 2, 3]:
        data = np.load(f"npz/p9_{slits}_slit.npz")
        m = int(np.sqrt(data.shape[1] - 1))
        t = data[:, 0]

        t_index = time_to_index(0.002)

        p = data[:, 2:].reshape(len(t), m, m)
        x_index = position_to_index(0.8)

        screen = p[t_index, :, x_index]
        screen /= np.sum(screen)  # Normalise

        axis = np.linspace(h, 1-h, m)

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=axis, y=screen,
                                 mode="lines", line=dict(width=5))
                      )
        fig.update_layout(
            font_family="Garamond",
            font_size=30,
            title=f"Probability distribution along y-axis at x = 0.8, t = 0.002, for {nice(slits)}.",
            xaxis_title="Position along y-axis",
            yaxis_title="Probability")
        fig.show()


def plot_potential(slits, h=0.005):
    """ Plot the potential in unit space """
    V = np.load(f"npz/potential_{slits}_slits.npz")
    m = V.shape[0]
    x = np.linspace(h, 1-h, m)
    fig = go.Figure(go.Contour(x=x, y=x, z=V))
    fig.show()


if __name__ == "__main__":
    # p7_no_slit()
    # p7_double_slit()
    # animate_real_imag()
    p8_three_steps()
    # animate()
    # p9(1)
    # p9(2)
    # p9(3)
