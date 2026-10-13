
import numpy as np
import matplotlib.pyplot as plt

print("Computational Nuclear Physics")
print("NumPy version:", np.__version__)

time = np.linspace(0, 10, 100)
N0 = 1000
decay_constant = 0.1

N = N0 * np.exp(-decay_constant * time)

plt.plot(time, N)
plt.xlabel("Time (s)")
plt.ylabel("Remaining nuclei")
plt.title("Radioactive Decay")
plt.grid(True)
plt.show()
