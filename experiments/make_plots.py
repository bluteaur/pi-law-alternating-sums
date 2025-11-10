import matplotlib
matplotlib.use("Agg")  # headless; no Tk

import numpy as np
import matplotlib.pyplot as plt
from mpmath import mp
from altpower.core import S_trunc, slope_via_fit

mp.dps = 80

def main():
    a_vals = np.linspace(6.5, 8.0, 16)
    N, p, q = 2000, 9, 4

    S = [S_trunc(a, N, p, q) for a in a_vals]
    B, a_bar = slope_via_fit(lambda a, N_, p_, q_: S_trunc(a, N_, p_, q_),
                             a_vals[0], a_vals[-1], len(a_vals), N, p, q)
    # y = A + B (a - a_bar)
    A = (mp.fsum(S) / len(S)) - B * (mp.mpf(np.mean(a_vals)) - a_bar)
    S_fit = [A + B * (mp.mpf(a) - a_bar) for a in a_vals]

    # Cast to float for matplotlib
    Sf = np.array([float(x/1e16) for x in S])
    Sfitf = np.array([float(x/1e16) for x in S_fit])

    plt.figure(figsize=(7, 4))
    plt.plot(a_vals, Sf, "o", label="S_N(a)")
    plt.plot(a_vals, Sfitf, "-", label="linear fit")
    plt.xlabel("a"); plt.ylabel("S_N(a) / 1e16")
    plt.grid(True); plt.legend(); plt.tight_layout()
    plt.savefig("figures/S_vs_a_linear_fit.png", dpi=200, bbox_inches="tight")
    print("Saved: figures/S_vs_a_linear_fit.png")

if __name__ == "__main__":
    main()
