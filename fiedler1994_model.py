import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt


def N(a, b):
    """
    (A6) from Fiedler+1994.
    """
    return 0.5*(erf(b/np.sqrt(2)) - erf(a/np.sqrt(2)))


def Fiedler(d_array, theta_l_, theta_b_, mu, core_dominance):
    """
    (A5) from Fiedler+1994.
    """
    c = np.sqrt(8*np.log(2))
    d_array = d_array*mu
    int1 = N(c*(d_array - 0.5*theta_l_)/(1 + theta_b_**2)**0.5,
             c*(d_array + 0.5*theta_l_)/(1 + theta_b_**2)**0.5)
    int2 = N(c*(0.5*theta_l_ - d_array), 100.)
    int3 = N(-100., c*(-0.5*theta_l_ - d_array))
    return core_dominance*(int1 + int2 + int3) + (1. - core_dominance)


epochs = np.arange(-0.5, 0.5, 0.01)
epoch_of_min = (np.max(epochs) + np.min(epochs))/2
epochs_centered = epochs - epoch_of_min

# Model parameters: theta_l, theta_b, mu, core_dominance
trues = [1, 2, 6, 1.0]
flux = Fiedler(epochs, *trues)


fig, axes = plt.subplots(1, 1)
axes.plot(epochs, flux)
axes.set_title("Fiedler 1994")
axes.set_xlabel("Epoch (normalized and centered)")
axes.set_ylabel("Normalized flux density")
plt.grid()
fig.savefig("Fiedler1994.png", bbox_inches="tight", dpi=300)
plt.show()

