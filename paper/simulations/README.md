# Paper Simulations

Scripts that generate figures for the aggregate formation manuscript.

## Setup
Each script loads `common.jl` for shared parameter sets and plot utilities.

## Workflow
1. Run a simulation script: `julia fig01_aggregate_lifecycle.jl`
2. Output goes to `paper/figures/`
3. Figures sync to Overleaf via Dropbox (see below)

## Overleaf integration
The `paper/figures/` directory should be symlinked to or synced with your
Dropbox folder that Overleaf reads from. Set this up once:

```bash
# Option A: Symlink (if Dropbox folder exists)
ln -s /path/to/SoilAggregateModel/paper/figures ~/Dropbox/Apps/Overleaf/YOUR_PROJECT/figures

# Option B: Copy script (add to your workflow)
# cp paper/figures/*.pdf ~/Dropbox/Apps/Overleaf/YOUR_PROJECT/figures/
```

In your LaTeX manuscript, reference figures as:
```latex
\includegraphics[width=\linewidth]{figures/fig01_aggregate_lifecycle.pdf}
```

## Repair log

**2026-07-28 — `single_aggregate_physics/run_simulations.jl` was broken.** Its
fungal diagnostic block read `trans.trans_i`, `trans.trans_n` and
`trans.insulation` from the value returned by `fungal_transitions`. That
function has returned `(immobil_i, immobil_n, Resp_F_conv)` since the February
2026 rewrite, so the script threw on the first diagnostic time and could not
have produced output since. Two further defects in the same block: the printed
`net` formula had α and β swapped relative to `fungi.jl`, and `ζ_T` was not
clamped to 1 after the Arrhenius factor the way `reactions.jl:116` clamps it.

All three are fixed. **The script has not been re-run** — the `.log` and `data/`
files in that directory predate the repair and should be regenerated before any
figure is taken from them.

Note also that the workflow section above names `fig01_aggregate_lifecycle.jl`,
which does not exist in this directory.
