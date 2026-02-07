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
