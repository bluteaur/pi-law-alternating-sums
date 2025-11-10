from mpmath import mp
from .core import S_trunc, S_abel, slope_via_fit, leading_slope, estimate_C1

mp.dps = 80

def run_suite(pq_pairs=((9,4), (7,3), (11,5)),
              Ns=(500, 1000, 1500, 2000, 3000),
              a_interval=("6.5", "8.0"),
              m_points=16):
    a_lo, a_hi = mp.mpf(a_interval[0]), mp.mpf(a_interval[1])
    print("\n=== Universal limit check (corrected): |B_N| / (π N^{p-q-1}) → (q/2π) ===\n")
    for (p, q) in pq_pairs:
        print(f"(p,q)=({p},{q}),  a ∈ [{a_lo}, {a_hi}],  m={m_points}")
        target = mp.mpf(q) / (2 * mp.pi)
        print(f"  Target q/(2π) = {mp.nstr(target, 12)}")
        B_vals, N_vals = [], []
        for N in Ns:
            B, a_bar = slope_via_fit(lambda a, N_, p_, q_: S_trunc(a, N_, p_, q_),
                                     a_lo, a_hi, m_points, N, p, q)
            N_vals.append(N); B_vals.append(B)
            norm = abs(B) / (mp.pi * (mp.mpf(N)**(p - q - 1)))
            pred = leading_slope(N, a_bar, p, q)
            kernel_factor = (mp.mpf(N)/(mp.mpf(N)+a_bar))**(q+1)
            print(f"  N={N:4d}  slope={mp.nstr(B, 12)}")
            print(f"        |B|/(π N^{p-q-1})={mp.nstr(norm, 12)}   "
                  f"(expected ≈ {mp.nstr(target*kernel_factor, 12)})   "
                  f"pred_ratio={mp.nstr(B/pred, 12)}")
        C1, rms = estimate_C1(N_vals, B_vals, a_bar, p, q)
        print(f"  → C1 estimate (B_N - L(N) ≈ C1·N^{p-q-2}): {mp.nstr(C1, 12)}   (RMS rel. err={mp.nstr(rms, 12)})\n")

def run_abel(a_interval=("6.5", "8.0"), m_points=16):
    a_lo, a_hi = mp.mpf(a_interval[0]), mp.mpf(a_interval[1])
    print("=== Abel-regularized slopes (small, convergent) for (p,q)=(9,4) ===")
    for r in ("0.95", "0.98", "0.99", "0.995", "0.999"):
        B_abel, _ = slope_via_fit(lambda a, p_, q_, r_: S_abel(a, p_, q_, r_),
                                  a_lo, a_hi, m_points, 9, 4, r)
        print(f"  r={r}  slope ≈ {mp.nstr(B_abel, 12)}")
