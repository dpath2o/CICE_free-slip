# Lateral-drag manuscript: notes towards circulation

30 October 2026

## Aim and timing

I'll bring the manuscript to a complete draft for Alex to read first, with the aim of circulating it to the full co-author team in early December. The purpose of these notes is to set out the remaining scientific questions and the work I'll use to resolve them. Comments on gaps in the argument or essential comparisons would be very welcome.

Before embarkation, I'll check the experiment coverage and prepare the analysis material. During the two-week voyage, as time allows, I'll work on the methods, discussion and figure captions. After the voyage, I'll finish the essential diagnostics and quantitative results through late November, leaving time for Alex's review before wider circulation.

## Comparisons I'll bring together

The existing experiments cover drag-coefficient sensitivity, coastline and grounded-iceberg contributions, and rheology. I've submitted three additional runs: `free-slip-no-drag`, `no-slip-LFI-match` and `Cs-high-ktens-zero`. These are the targeted additions I'll use to separate effects that the existing configurations could not isolate cleanly.

| Question | Comparisons I'll use |
|---|---|
| What changes with the coastal boundary condition alone? | `no-slip-LFI-match` against `free-slip-no-drag`, with matched rheology and lateral drag disabled. |
| What does lateral drag add under free slip? | `free-slip-no-drag` against `Cs-high` and `Cq-high`. |
| How sensitive are the results to drag strength and speed dependence? | The retained low, mid and high static and quadratic coefficient sequences. |
| What do coastline and grounded-iceberg geometry contribute? | `cst-drag`, `GIB-drag` and combined `Cs-high`; their responses need not add linearly. |
| How much does tensile strength matter? | `Cs-high-ktens-zero`, `Cs-high-ktens-low`, `Cs-high-ktens-mid` and `Cs-high`, spanning `Ktens=0`, `0.05`, `0.1` and `0.2`. |
| How robust is the response to other rheological choices? | The ellipse-ratio and ice-strength alternatives, with combined changes identified explicitly. |

All three new runs are intended to match the grid, forcing and other physical settings of `Cs-high`, apart from the stated changes. They start without ice in January 1995 and are planned through December 2005. I'll check their completion, configuration and output coverage before using them. I'll retain 1995–1999 as spin-up and 2000–2005 for model diagnostics.

I'll keep the older `no-slip-def` and `no-slip-LFI` as reference configurations, since their rheology differs from the new matched controls. I'll use the new zero-tensile run to isolate tensile strength: `Cs-high-def-rhe` changes ellipse ratios and numerical settings as well. I'll distinguish the matched ice-strength comparison (`Cs-high` versus `Cs-high-roth-rhe-high`) from the combined rheology alternatives.

## Evidence I'll check

I'll reconcile the experiment settings, manuscript tables and figure labels against the archived configurations. I'll also pull together the existing implementation evidence for the free-slip boundary treatment, C-grid mapping, stable momentum updates and drag opposing motion. Any further implementation test should address a specific unresolved issue.

For observational comparisons, I'll use **April 2000–April 2005**, confirming the exact endpoint before final exports. I'll regenerate older panels that extend through December 2005 and check that all reported statistics use the stated dates, units and common coverage. I'll pay particular attention to the regional FIA signs and totals, and to the temporal support of the fortnightly observations when comparing fast-ice persistence.

I'll quantify the matched changes in fast-ice area, seasonal growth and retreat, and spatial persistence. I'll then use thickness, deformation, ice strength and the archived lateral-drag stresses to interpret the response. I'll check whether improvements near anchors come with changes in the surrounding pack, and narrow the claims where a diagnostic is unavailable. I'll use classification sensitivity and interannual variability to assess how much weight the apparent improvements can bear.

The detailed checks, processing recipes and analysis steps now sit in two publication notebooks:

- [Main manuscript comparisons](https://github.com/dpath2o/CICE_free-slip/blob/develop/notebooks/LD-pub-workspace.ipynb)
- [Supporting spatial and mechanical diagnostics](https://github.com/dpath2o/CICE_free-slip/blob/develop/notebooks/LD-performance.ipynb)

## Scope and manuscript completion

I'll keep the paper focused on the formulation, its CICE implementation and the retained Antarctic experiments. Coupled-ocean experiments, waves and tides, Arctic evaluation, resolution sensitivity, a new iceberg dataset and extension to 2018 remain outside this draft. I'll discuss these as limitations or future work where relevant.

I'll replace the remaining development figures and placeholders with checked exports, reconcile regional names and captions, and move repetitive sector panels to the supplement where useful. I'll complete the discussion around the contribution beyond land-masked icebergs, parameter trade-offs, seasonal and regional errors, forcing and geometry limitations, and implications for coupled models.

I'll finish the conclusions and revise the abstract around verified quantitative findings. The central claim I'll assess is whether the lateral-drag formulation represents an important component of Antarctic fast-ice anchoring without adding grounded icebergs to the model land mask, and under which parameter choices and limitations that claim holds.

Before sending the draft to Alex, I'll check the references, figure resources, authorship and dataset attribution, funding and code/data statements, and compile the full manuscript and supplement. I'll circulate after the main argument and conclusions are complete, the essential comparisons are reproducible, and the claims reflect the evidence. I'll identify any remaining scientific questions explicitly for the co-authors.
