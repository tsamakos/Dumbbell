import numpy as np
import matplotlib.pyplot as plt


N = np.array([5, 10, 20, 40])
Dpar = np.array([0.183481778, 0.184029609, 0.178087061, 0.165786456])
Dver = np.array([0.16302855, 0.106654952, 0.134356669, 0.148948004])



plt.figure()
plt.plot(N, Dpar, 'o-', label="D$_{\\parallel}$")
plt.plot(N, Dver, 'o-', label="D$_{\\perp}$")
plt.xlabel("N")
plt.ylabel("D")
plt.legend()
plt.savefig("D.png")

# plot ratio
plt.figure()
plt.plot(N, Dver/Dpar, 'o-')
plt.xlabel("N")
plt.ylabel("D$_{\\perp}$/D$_{\\parallel}$")
plt.savefig("D_ratio.png")

