from mpmath import mp
import csv
from altpower.core import S_trunc, slope_via_fit, leading_slope

mp.dps = 80

def one_case(p, q, Ns=(500,1000,1500,2000,3000), a_lo="6.5", a_hi="8.0", m=16):
    a_lo, a_hi = mp.mpf(a_lo), mp.mpf(a_hi)
    rows = [("N","B","norm","expected","pred_ratio")]  # header row

    for N in Ns:
        B, a_bar = slope_via_fit(lambda a, N_, p_, q_: S_trunc(a, N_, p_, q_),
                                 a_lo, a_hi, m, N, p, q)

        norm = abs(B) / (mp.pi * (mp.mpf(N)**(p - q - 1)))
        kernel_factor = (mp.mpf(N)/(mp.mpf(N)+a_bar))**(q+1)
        expected = (mp.mpf(q)/(2*mp.pi)) * kernel_factor
        pred = leading_slope(N, a_bar, p, q)

        rows.append((
            int(N),
            mp.nstr(B, 12),
            mp.nstr(norm, 12),
            mp.nstr(expected, 12),
            mp.nstr(B/pred, 12),
        ))
    return rows

def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        for r in rows:
            w.writerow(r)

def write_latex_table(path, rows, caption=""):
    header, *body = rows
    with open(path, "w") as f:
        f.write("\\begin{table}[H]\n\\centering\n")
        f.write("\\begin{tabular}{rllll}\n\\toprule\n")
        f.write(" & ".join(header) + " \\\\\n\\midrule\n")
        for r in body:
            f.write(f"{r[0]} & {r[1]} & {r[2]} & {r[3]} & {r[4]} \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
        if caption:
            f.write(f"\\caption{{{caption}}}\n")
        f.write("\\end{table}\n")

if __name__ == "__main__":
    cases = [(9,4), (7,3), (11,5)]
    for p,q in cases:
        rows = one_case(p,q)
        write_csv(f"figures/table_p{p}_q{q}.csv", rows)
        write_latex_table(
            f"figures/table_p{p}_q{q}.tex",
            rows,
            caption=f"Results for $(p,q)=({p},{q})$."
        )
    print("Saved CSV + LaTeX tables in `figures/`")
