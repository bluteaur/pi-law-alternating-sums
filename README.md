# README.md
## Alternating Power Sums — Asymptotics and a Universal q/(2π) Constant

This repo reproduces the asymptotic slope law for truncated alternating sums
S_N(a)=∑_{k=1}^N (-1)^k k^p/(k+a)^q, p>q≥1:
B_N(a) ~ -(q/2)·(N/(N+a))^{q+1}·N^{p-q-1}, hence |B_N|/(π N^{p-q-1})→q/(2π).

### Quick start
```bash
python -m experiments.run_asymptotics
python -m experiments.make_plots
