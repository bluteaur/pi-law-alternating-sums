import os
from typing import Callable, Iterable, Tuple
from mpmath import mp

# Precision: allow override via env, default 80 dp
mp.dps = int(os.getenv("ALTPWR_DPS", "80"))

def S_trunc(a: mp.mpf, N: int, p: int, q: int) -> mp.mpf:
    """Truncated alternating sum: S_N(a) = sum_{k=1}^N (-1)^k k^p / (k+a)^q."""
    a = mp.mpf(a)
    tot = mp.mpf("0")
    sign = -1
    for k in range(1, N + 1):
        tot += sign * (mp.power(k, p)) / mp.power(k + a, q)
        sign = -sign
    return tot

def S_abel(a: mp.mpf, p: int, q: int, r: str | float = "0.999", K: int = 500_000) -> mp.mpf:
    """Abel-regularized infinite sum:
       S_Abel(a; r) = sum_{k>=1} (-1)^k r^k k^p / (k+a)^q, 0<r<1.
       Converges absolutely; use to verify small-slope limit and sign.
    """
    a = mp.mpf(a)
    r = mp.mpf(r)
    if not (mp.mpf("0") < r < mp.mpf("1")):
        raise ValueError("r must be in (0,1) for Abel regularization")
    tot = mp.mpf("0")
    rk = -r  # (-1)*r^1
    for k in range(1, K + 1):
        tot += rk * mp.power(k, p) / mp.power(k + a, q)
        rk *= -r
        # tiny early break (optional)
        if mp.almosteq(rk, 0, rel_eps=0, abs_eps=mp.mpf("1e-80")):
            break
    return tot

def _linreg_through_window(a_lo: mp.mpf, a_hi: mp.mpf, m: int,
                           f: Callable[[mp.mpf], mp.mpf]) -> Tuple[mp.mpf, mp.mpf]:
    """Return (slope, a_bar) for a least-squares line fit over m evenly-spaced points in [a_lo, a_hi]."""
    a_lo, a_hi = mp.mpf(a_lo), mp.mpf(a_hi)
    if m < 2:
        raise ValueError("m must be >= 2 for a linear fit")
    h = (a_hi - a_lo) / (m - 1)
    xs = [a_lo + i*h for i in range(m)]
    ys = [f(x) for x in xs]
    xbar = mp.fsum(xs) / m
    ybar = mp.fsum(ys) / m
    Sxx = mp.fsum((x - xbar)*(x - xbar) for x in xs)
    if Sxx == 0:
        raise ZeroDivisionError("degenerate abscissae in linear fit")
    Sxy = mp.fsum((xs[i] - xbar)*(ys[i] - ybar) for i in range(m))
    slope = Sxy / Sxx
    return slope, xbar  # a_bar

def slope_via_fit(S_fn: Callable[..., mp.mpf],
                  a_lo: mp.mpf, a_hi: mp.mpf, m: int,
                  *S_args) -> Tuple[mp.mpf, mp.mpf]:
    """Estimate ∂S/∂a by linear fit on a ∈ [a_lo, a_hi] using m points.
       S_fn(a, *S_args) is evaluated at each point.
       Returns (slope, a_bar).
    """
    return _linreg_through_window(a_lo, a_hi, m, lambda a: S_fn(a, *S_args))

def leading_slope(N: int, a_bar: mp.mpf, p: int, q: int) -> mp.mpf:
    """Heuristic leading term for B_N(a) ≈ - (q/2) * (N/(N+a_bar))^{q+1} * N^{p-q-1}.
       Sign is negative (since derivative wrt a lowers denominator power).
    """
    N = mp.mpf(N); a_bar = mp.mpf(a_bar)
    kernel = mp.power(N / (N + a_bar), q + 1)
    return - (mp.mpf(q) / 2) * kernel * mp.power(N, p - q - 1)

def estimate_C1(N_list: Iterable[int],
                B_list: Iterable[mp.mpf],
                a_bar: mp.mpf, p: int, q: int) -> Tuple[mp.mpf, mp.mpf]:
    """Fit correction B_N - L(N) ≈ C1 * N^{p-q-2}. Returns (C1, RMS_relative_error)."""
    N_list = list(N_list)
    B_list = list(B_list)
    xs = [mp.power(mp.mpf(N), p - q - 2) for N in N_list]
    ys = [B_list[i] - leading_slope(N_list[i], a_bar, p, q) for i in range(len(N_list))]
    num = mp.fsum(xs[i] * ys[i] for i in range(len(xs)))
    den = mp.fsum(xs[i] * xs[i] for i in range(len(xs)))
    C1 = num / den if den != 0 else mp.nan
    errs = []
    for i, N in enumerate(N_list):
        model = leading_slope(N, a_bar, p, q) + C1 * mp.power(mp.mpf(N), p - q - 2)
        denom = max(mp.mpf("1e-30"), abs(B_list[i]))
        errs.append(mp.power((B_list[i] - model) / denom, 2))
    rms_rel = mp.sqrt(mp.fsum(errs) / len(errs)) if errs else mp.nan
    return C1, rms_rel
