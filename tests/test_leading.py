from mpmath import mp
from altpower.core import S_trunc, slope_via_fit, leading_slope

mp.dps = 60

def test_ratio_near_q_over_2pi():
    a_lo, a_hi, m = mp.mpf("6.5"), mp.mpf("8.0"), 16
    for (p,q,N) in [(9,4,2000),(7,3,2000),(11,5,2000)]:
        B, a_bar = slope_via_fit(lambda a, N_, p_, q_: S_trunc(a, N_, p_, q_),
                                 a_lo, a_hi, m, N, p, q)
        norm = abs(B) / (mp.pi * (mp.mpf(N)**(p - q - 1)))
        target = mp.mpf(q)/(2*mp.pi) * (mp.mpf(N)/(mp.mpf(N)+a_bar))**(q+1)
        # within ~0.2% at N=2000
        assert abs(norm/target - 1) < mp.mpf("2e-3")
