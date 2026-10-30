# Lateral-drag publication notebooks

Two notebooks own the analysis for this publication, both maintained directly on `develop`:

| Notebook | Purpose |
|---|---|
| [LD-pub-workspace.ipynb](LD-pub-workspace.ipynb) | Experiment/namelist inventory, processing recipes, analytical form functions, FIA statistics and climatology, matched controls, coefficient/geometry/rheology comparisons, classification robustness, exports and manuscript review. |
| [LD-performance.ipynb](LD-performance.ipynb) | Common-support FIP skill and paired panels; thickness, deformation, strength and actual face stress; occurrence/retention; anchor-distance context; pack-ice checks and optional point-observation context. |

`paper2.ipynb` has been retired; its original content remains in Git history. Its relevant scientific questions are covered in these two notebooks. Unrelated wave-forcing and pilot work are outside this publication. Detailed manuscript tasks belong in the relevant notebook; `lateral-drag_manuscript_notes.md` is the short co-author-facing outline.

## Running on Gadi

1. Load the analysis environment and make the standalone `shuga` package importable (the experiment-group path changes merged in `caae0c0768414be30f196754e9d8b57621b7e154` or newer). Set `SHUGA_REPO` to its checkout if it is not installed.
2. Start Jupyter from this directory, or set `CICE_FREE_SLIP_REPO` to the checkout. Review the single shared configuration in `ld_pub_config.json`.
3. Run the first notebook's inventory before enabling expensive sections. Confirm archive coverage and select metadata explicitly if several `ice_in` records exist. The three submitted runs remain pending until outputs are present and validated.
4. Enable the required `RUN_*` sections after reviewing their inputs and definitions. Run expensive reductions in an allocated analysis/PBS job; the notebooks never submit jobs automatically. Processed inputs are read-only. Export tables, figures and small caches go to `GRAPHICAL/lateral-drag/publication`; override with `LD_PUB_EXPORT` for a local trial.
5. Run the supporting notebook for the selected comparisons. Supply the explicitly documented native-grid and common-grid products. A review flag records a scientific decision; it does not perform that validation for you.

Defaults execute the inventory and write recipes/review records, with scientific stages disabled. They do **not** imply that results have been extracted or verified. Restart and run all after changing configuration; do not mix stale interactive state with new settings.

The model-diagnostic window is January 2000–December 2005, after the 1995–1999 spin-up. The attached analysis notes specify **April 2000–April 2005** for observational comparison, interpreted here as 1 April 2000–30 April 2005 inclusive. Confirm the endpoint before final export. Older panels ending December 2005 need regeneration; partial endpoint years are labelled. Observation composites require a reviewed temporal-support treatment: independently averaged FIP fields are not accepted as a common-support skill comparison.

`ld_pub_workflow.py` contains shared validation and reductions, rather than a third notebook. Test it with `python3 -m pytest tests/test_ld_pub_workflow.py` from the repository root. Reductions reject unverified units, incomplete daily support and stale analysis signatures. Cache reuse assumes immutable source stores: a directory modification time is not a hash of every Zarr chunk. Leave recomputation enabled after any source rewrite. Gadi scientific data and the complete plotting environment are required to validate enabled scientific stages.
