import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

file = "data.txt"
data = pd.read_csv(file, header=0, sep=" ")
x = data["x"]
u = data["u"]

plt.plot(x, u, lw=5)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Analytical solution")
plt.show()