import plotly.express as px
import pandas as pd

file = "data.txt"
data = pd.read_csv(file, header=0, sep=" ")

fig = px.line(data, x="x", y=data.keys()[1:], title="Analytical solutions")
fig.update_layout(
    font_family="Garamond",
    font_size=25,
)
fig.show()
