# CICE free-slip and lateral-drag development branch

**Version:** v1.0 
**Status:** research-development branch for CICE6 free-slip boundary conditions and lateral-drag parameterisation
**Primary application:** Antarctic landfast sea-ice process experiments

This repository develops and evaluates a free-slip coastal boundary condition and explicit lateral-drag parameterisations for representing Antarctic landfast sea ice in CICE.

The central physical problem is that Antarctic landfast sea ice is not simply slow pack ice. It forms where coastal geometry, ice shelves, islands, grounded icebergs, internal ice stress, and local deformation interact to immobilise otherwise mobile sea ice. This branch tests whether high-resolution coastline and grounded-iceberg form factors can provide geometric anchoring in CICE without relying on a no-slip boundary condition that directly suppresses coastal ice motion.

---

## Contents

- [Version v1.0 scope](#version-v10-scope)
- [Scientific motivation](#scientific-motivation)
- [Key source files](#key-source-files)
- [Boundary-condition motivation](#boundary-condition-motivation)
- [C-grid free-slip implementation](#c-grid-free-slip-implementation)
- [Free-slip strain-rate expectation](#free-slip-strain-rate-expectation)
- [Lateral-drag form factors](#lateral-drag-form-factors)
- [Lateral-drag form functions](#lateral-drag-form-functions)
- [`blend_strain` formulation](#blend_strain-formulation)
- [Strain-rate normalisation](#strain-rate-normalisation)
- [Diagnostics](#diagnostics)
- [Namelist parameters](#namelist-parameters)
- [v1.0 test-configuration experiments](#v10-experiment-configuration)
- [v2.0 publication experiments](#v20-experiment-publication-experiments)
- [Ice-strength thickness-exponent tests](#ice-strength-thickness-exponent-tests)
- [Analysis workflow](#analysis-workflow)
- [Known limitations](#known-limitations)
- [Relationship to upstream CICE](#relationship-to-upstream-cice)
- [License](#license)

---

## Version v1.0 scope

Version v1.0 marks the first internally consistent version of this branch in which:

1. C-grid free-slip strain-rate diagnostics compute `divergU`, `tensionU`, `shearU`, and `DeltaU`.
2. `DeltaU` is available to the lateral-drag `blend_strain` form function as an area-normalised effective strain-rate scale.
3. Lateral-drag diagnostics are written only on the final EVP subcycle to avoid excessive memory traffic.
4. `blend_exp` is supported through an integer exponent pathway suitable for sensitivity experiments.
5. Diagnostic fields are available for process-based evaluation of static, quadratic, linear, and `blend_strain` form functions.
6. The branch is ready for multi-month and multi-year Antarctic fast-ice experiments using free-slip dynamics and static coastline/grounded-iceberg form factors.

This is a research branch. It is not an official CICE Consortium release.

---

## Scientific motivation

No-slip boundary conditions can artificially favour coastal immobilisation by forcing tangential ice velocity to vanish at coastlines and grounded obstacles. That can help produce landfast ice, but it conflates boundary-condition imposition with emergent landfast-ice physics.

This branch instead uses a free-slip boundary condition and introduces an explicit lateral-drag stress of the form

$$
\boldsymbol{\tau}_{\mathrm{LD}} = - \mathrm{K}_u \, \phi \, \mathbf{u},
$$

where:

- $K_u$ is an effective mass/form-factor term based on local ice/snow mass and static geometric form factors;
- $\phi$ is a lateral-drag form function;
- $\mathbf{u}$ is the relevant C-grid velocity component.

The intent is that coastline and grounded-iceberg geometry define **where anchoring is possible**, while the form function determines **when the local ice state is dynamically eligible for locking**.

---

## Key source files

The main implementation files are:

```text
cicecore/cicedyn/dynamics/ice_dyn_evp.F90
cicecore/cicedyn/dynamics/ice_dyn_shared.F90
cicecore/cicedyn/general/ice_init.F90
cicecore/cicedyn/analysis/ice_history.F90
cicecore/cicedyn/analysis/ice_history_shared.F90
```

The main analysis-side companion code is in the `shuga` toolbox (<https://github.com/dpath2o/mawsons-chest/tree/main/shuga>):

```text
shuga/grid/lateral_drag.py
```

This is used to generate and evaluate high-resolution coastline and grounded-iceberg form factors.

---

## Boundary-condition motivation

### No-slip

In an idealised no-slip configuration, ice velocity is constrained at solid boundaries. For a coastal boundary, this suppresses tangential motion and can promote artificial coastal immobilisation.

### Free-slip

A free-slip boundary condition removes the tangential no-slip constraint. It permits sea ice to slide along coastlines and grounded-iceberg boundaries while retaining the appropriate impermeability condition. This is a better framework when the goal is to let landfast ice emerge from the resolved momentum balance rather than impose it directly at the boundary.

However, free-slip also removes an implicit source of unresolved lateral resistance. This motivates an explicit lateral-drag term.

---

## C-grid free-slip implementation

The C-grid implementation uses staggered velocity components on E and N faces. The EVP dynamics compute strain rates, internal stresses, stress divergence, and then update the E- and N-face velocity components through:

```fortran
stepu_C
stepv_C
```

in:

```text
cicecore/cicedyn/dynamics/ice_dyn_shared.F90
```

The EVP driver in:

```text
cicecore/cicedyn/dynamics/ice_dyn_evp.F90
```

calls the relevant strain-rate routine depending on the selected boundary condition, then passes the U-grid deformation invariant `DeltaU` into the lateral-drag velocity update.

---

## Free-slip strain-rate expectation

The free-slip strain-rate routine should preserve a rigid-translation, zero-deformation state on a uniform grid.

For a two-dimensional velocity field

$$
\mathbf{u}=(u,v),
$$

the strain-rate components are

$$
e_{11}=\frac{\partial u}{\partial x},
\qquad
e_{22}=\frac{\partial v}{\partial y},
\qquad
2e_{12}=\frac{\partial u}{\partial y}+\frac{\partial v}{\partial x}.
$$

The CICE-style diagnostics are

$$
\mathrm{div} = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y},
$$

$$
\mathrm{tension} = \frac{\partial u}{\partial x} - \frac{\partial v}{\partial y},
$$

$$
\mathrm{shear} = \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x},
$$

and the EVP deformation invariant is

$$
\Delta U = \sqrt{\mathrm{div}^{2} + e_{\mathrm{fac}}\left(\mathrm{tension}^{2} + \mathrm{shear}^{2}\right)}.
$$

For rigid translation,

$$
u(x,y)=U_0,
\qquad
v(x,y)=V_0,
$$

all spatial derivatives vanish, so

$$
\mathrm{div} = \mathrm{tension} = \mathrm{shear} = \Delta U = 0.
$$

On a uniform orthogonal C-grid,

$$
dx_E=dx_U=\Delta x,
\qquad
dy_N=dy_U=\Delta y,
$$

and the metric-gradient correction terms vanish. The free-slip routine uses even reflection at masked neighbouring faces. For example,

$$
uN_{i+1,j} = uvelN(i+1,j)\,npm(i+1,j) + \left[npm(i,j)-npm(i+1,j)\right]\,npm(i,j)\,uvelN(i,j).
$$

If both neighbouring faces are active, this returns the neighbouring value. If the neighbouring face is masked and the interior face is active, the reflected value becomes the interior value. Therefore, for uniform velocity,

$$
uN_{i+1,j}=uN_{ij}=U_0,
\qquad
vE_{i,j+1}=vE_{ij}=V_0,
$$

and similarly for the shear terms. Thus,

$$
\mathrm{divergU} = \mathrm{tensionU} = \mathrm{shearU} = \Delta U = 0.
$$

This is a verification condition for idealised uniform-grid tests. It does **not** imply that free-slip always produces zero deformation. In realistic Antarctic simulations, nonzero `DeltaU` is expected because of spatially variable forcing, internal stress, coastlines, grounded icebergs, form factors, ice thickness gradients, concentration gradients, and ocean-current gradients.

---

## Lateral-drag form factors

Lateral drag is applied where high-resolution coastline and grounded-iceberg geometry indicate unresolved lateral contact or anchoring potential.

The static geometric factors are represented on the C-grid velocity faces as:

```text
F2E
F2N
```

The effective lateral-drag mass/form factor is

$$
\mathrm{K}_u = \mathrm{M}_u \mathrm{F}_2,
$$

where $\mathrm{M}_u$ is the velocity-grid ice/snow mass and $\mathrm{F}_2$ is the local geometric form factor. This makes the lateral stress scale with both available ice mass and local geometric exposure to coastline or grounded-iceberg anchoring.

Compared with a coastline-only implementation, this branch explicitly supports form-factor contributions from grounded icebergs. Grounded icebergs are treated as additional static anchoring geometry, with geometric properties that can be incorporated into the high-resolution form-factor generation workflow.

---

## Lateral-drag form functions

The lateral-drag stress is applied as

$$
\boldsymbol{\tau}_{\mathrm{LD}} = - \mathrm{K}_u \, \phi \, \mathbf{u}.
$$

The scalar form function $\phi$ is selected by `form_func`.

### `static`

$$
\phi_{\mathrm{static}} = \frac{\mathrm{C}_{\mathrm{s}}}{|\mathbf{u}| + u_0}.
$$

This is the Liu-style reference branch. It produces strong low-speed damping while remaining regularised by $u_0$.

### `quad`

$$
\phi_{\mathrm{quad}} = \mathrm{C}_\mathrm{q} |\mathbf{u}|.
$$

This is a velocity-scaled mobile branch. It is weak near zero speed and stronger for faster-moving ice.

### `linear`

$$
\phi_{\mathrm{linear}} = \mathrm{C}_\mathrm{L}.
$$

This is a Rayleigh-style damping branch.

### `blend_strain`

The corrected `blend_strain` form function blends the static and quadratic branches using both speed and strain-rate gates.

---

## `blend_strain` formulation

The corrected `blend_strain` formulation computes:

$$
w_\epsilon = \frac{1} {1+\left(\epsilon_{\mathrm{eff}}/\epsilon_{\mathrm{blend}}\right)^p},
$$

$$
w_u = \frac{1}{1+\left(|\mathbf{u}|/u_{\mathrm{blend}}\right)^p},
$$

$$
w_{\mathrm{L}} = w_\epsilon w_u.
$$

The realised form function is

$$
\phi_{\mathrm{blend}} = w_{\mathrm{L}}\phi_{\mathrm{static}} + \left(1-w_{\mathrm{L}}\right)\phi_{\mathrm{quad}}.
$$

Thus:

```text
low speed + low strain rate  -> static/locking branch
high speed or high strain    -> quadratic/mobile branch
```

This formulation is designed to permit landfast ice to emerge where the ice is slow, coherent, and geometrically anchored, while allowing mobile or deforming pack ice to remain dynamically mobile.

For computational efficiency, `blend_exp` is treated as an integer-valued exponent in the optimized path. The intended supported values are:

```text
blend_exp = 1, 2, 3, 4, 5, ..., 30
```

---

## Strain-rate normalisation

The C-grid diagnostic `DeltaU` is an EVP deformation invariant multiplied by U-cell area. Therefore, the effective strain rate used by `blend_strain` is area-normalised before comparison with `eps_blend`.

For E-face velocity points:

$$
\epsilon_{\mathrm{eff,E}}(i,j) = \frac{\Delta U(i,j)+\Delta U(i,j-1)}{A_U(i,j)+A_U(i,j-1)}.
$$

For N-face velocity points:

$$
\epsilon_{\mathrm{eff,N}}(i,j) = \frac{\Delta U(i,j)+\Delta U(i-1,j)}{A_U(i,j)+A_U(i-1,j)}.
$$

This produces a strain-rate scale in $s^{-1}$, suitable for comparison with `eps_blend`.

---

## Diagnostics

This branch includes lateral-drag diagnostics intended to support process-based evaluation rather than tuning only against aggregate fast-ice area.

Key diagnostics include:

| Diagnostic | Meaning |
|---|---|
| `F2E`, `F2N` | Static geometric form factors |
| `KuxE`, `KuyE`, `KuxN`, `KuyN` | Lateral-drag stress components |
| `ldphiE`, `ldphiN` | Realised form function $\phi$ |
| `ldwgtE`, `ldwgtN` | Static/locking branch weight $w_{\mathrm{lock}}$ |
| `ldepsE`, `ldepsN` | Effective strain rate used by `blend_strain` |
| `ldspdE`, `ldspdN` | Speed used by the form function |
| `ldpstatE`, `ldpstatN` | Static branch damping-rate diagnostic |
| `ldpquadE`, `ldpquadN` | Quadratic branch damping-rate diagnostic |
| `ldplinE`, `ldplinN` | Linear branch damping-rate diagnostic |

These diagnostics are written only on the final EVP subcycle. This avoids excessive memory traffic inside the dynamics loop while preserving daily or sub-daily diagnostic output for process analysis.

All new history-field flags must be included consistently in:

```text
ice_history_shared.F90
ice_history.F90
```

including declaration, namelist entry, broadcast, field definition, and accumulation. Missing broadcasts can cause MPI/PIO collective mismatches during history writing.

---

## Namelist parameters

The lateral-drag implementation uses the following dynamics namelist parameters.

| Parameter | Meaning |
|---|---|
| `boundary_condition` | Boundary condition, e.g. `'free_slip'` |
| `lateral_drag` | Enables/disables lateral drag |
| `form_func` | Selects `static`, `quad`, `linear`, or `blend_strain` |
| `Cs` | Static branch coefficient |
| `Cq` | Quadratic branch coefficient |
| `C_L` | Linear branch coefficient |
| `u0` | Low-speed regularisation |
| `u_blend` | Speed transition scale for `blend_strain` |
| `eps_blend` | Effective strain-rate transition scale for `blend_strain` |
| `blend_exp` | Blend transition exponent |
| `u_cap` | Optional speed cap, if enabled in the active implementation |

---

## v1.0 experiment configuration

The current v1.0 `LD-blend-base` candidate uses:

```fortran
ndte                 = 360
Pstar                = 2.75e4
Cstar                = 20
Ktens                = 0.2
e_yieldcurve         = 1.5
e_plasticpot         = 1.5

boundary_condition   = 'free_slip'
lateral_drag         = .true.
form_func            = 'blend_strain'

Cs                   = 2.5e-4
Cq                   = 350.0
C_L                  = 0.0
eps_blend            = 5.0e-8
blend_exp            = 3.0
u_blend              = 5.0e-4
```

The intended multi-year comparison set includes:

```text
LD-NIL
LD-static-Cs5e-4
LD-quad-Cq75
LD-linear-0p25
LD-blend-base
LD-blend-Cq75-exp20
```

These experiments are designed to test whether the corrected `blend_strain` formulation improves Antarctic fast-ice behaviour relative to no lateral drag, static drag, quadratic drag, and linear drag.

---

## v2.0 publication experiment suite

The v2.0 experiment suite records the main CICE6 free-slip and lateral-drag simulations used for the publication analysis. Experiment names refer to the archived experiment directories under:

```text
~/AFIM_archive/<experiment>/
boundary_condition = 'free_slip'
lateral_drag       = .true.
form factors       = coastline + grounded icebergs
history output     = daily averaged
dt                 = 1800 s
ndte               = 360
kdyn               = 1
```

The main reference `static` form function setting:
```text
form_func          = 'static'
Cs                 = 1.0e-3
Pstar              = 2.75e4
Cstar              = 20
Ktens              = 0.2
e_yieldcurve       = 1.5
e_plasticpot       = 1.5
```

| Order | Experiment              | Family                                       | Key perturbation                                            | Scientific role                                                                                            |
| ----: | ----------------------- | -------------------------------------------- | ----------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------- |
|     0 | `no-slip-def`           | Boundary reference                           | `no-slip`; no lateral drag; default rheology                | Traditional CICE lateral-boundary reference.                                                               |
|     1 | `no-slip-LFI`           | Boundary/rheology reference                  | `no-slip`; no lateral drag; LFI-favourable rheology         | Tests whether rheology alone can promote mechanically persistent ice under no-slip conditions.             |
|     2 | `Cs-low`                | Static lateral drag                          | `form_func = 'static'`; low `Cs`                            | Low-end static lateral-resistance case.                                                                    |
|     3 | `Cs-mid`                | Static lateral drag                          | `form_func = 'static'`; intermediate `Cs`                   | Intermediate static lateral-resistance case.                                                               |
|     4 | `Cs-high`               | Static lateral drag                          | `form_func = 'static'`; `Cs = 1.0e-3`                       | Main strong static-drag case; currently the central free-slip lateral-drag reference.                      |
|     5 | `Cs-high-def-rhe`       | Static lateral drag + rheology               | `Cs-high` with default rheology                             | Separates the effect of strong static drag from the higher-tensile/high-rheology settings.                 |
|     6 | `Cs-high-ktens-low`     | Static lateral drag + tensile sensitivity    | `Cs-high` with lower tensile-strength setting               | Tensile-strength sensitivity around the main static-drag case.                                             |
|     7 | `Cs-high-ktens-mid`     | Static lateral drag + tensile sensitivity    | `Cs-high` with intermediate tensile-strength setting        | Intermediate tensile-strength sensitivity around the main static-drag case.                                |
|     8 | `Cs-high-eDef`          | Static lateral drag + yield/plastic geometry | `Cs-high` with modified yield/plastic eccentricity settings | Tests sensitivity to the rheological ellipse/flow-rule geometry.                                           |
|     9 | `Cq-low`                | Quadratic lateral drag                       | `form_func = 'quad'`; low `Cq`                              | Low-end velocity-dependent lateral resistance.                                                             |
|    10 | `Cq-mid`                | Quadratic lateral drag                       | `form_func = 'quad'`; intermediate `Cq`                     | Intermediate quadratic response.                                                                           |
|    11 | `Cq-high`               | Quadratic lateral drag                       | `form_func = 'quad'`; high `Cq`                             | Strong velocity-dependent lateral resistance.                                                              |
|    12 | `Cl-low`                | Linear lateral drag                          | `form_func = 'linear'`; low `C_L`                           | Weak Rayleigh-style lateral-resistance case.                                                               |
|    13 | `Cl-mid`                | Linear lateral drag                          | `form_func = 'linear'`; intermediate `C_L`                  | Intermediate Rayleigh-style lateral-resistance case.                                                       |
|    14 | `blend-strain-low`      | `blend_strain` lateral drag                  | Low blend-strain setting                                    | Low-end blend-strain case.                                                                                 |
|    15 | `blend-strain-mid`      | `blend_strain` lateral drag                  | Intermediate blend-strain setting                           | Intermediate blend-strain case.                                                                            |
|    16 | `blend-strain-high`     | `blend_strain` lateral drag                  | High blend-strain setting                                   | Stronger blend-strain case.                                                                                |
|    17 | `cst-drag`              | Geometry sensitivity                         | Coastline-only form-factor contribution                     | Separates coastline anchoring from grounded-iceberg anchoring.                                             |
|    18 | `GIB-drag`              | Geometry sensitivity                         | Grounded-iceberg-only form-factor contribution              | Separates grounded-iceberg anchoring from coastline anchoring.                                             |
|    19 | `kstrength-test01`      | Ice-strength sensitivity                     | `kstrength = 2`; `Pstar = 27500 / sqrt(2)`                  | Tests an (h^{3/2}) ice-strength law normalised to match the standard Hibler strength at 2 m ice thickness. |
|    20 | `kstrength-test02`      | Ice-strength sensitivity                     | `kstrength = 2`; `Pstar = 27500 / sqrt(3)`                  | Tests an (h^{3/2}) ice-strength law normalised to match the standard Hibler strength at 3 m ice thickness. |
|    21 | `Cs-high-roth-def`      | Ice-strength/rheology sensitivity            | `Cs-high`; Rothrock strength; default rheology              | Separates Rothrock ice strength from the default-rheology static-drag case.                                |
|    22 | `Cs-high-roth-rhe-high` | Ice-strength/rheology sensitivity            | `Cs-high`; Rothrock strength; high-rheology settings        | Tests Rothrock strength under the high-rheology static-drag configuration.                                 |

The `static`, `quadratic`, `linear`, and `blend-strain` cases test the form-function dependence of the lateral-drag parameterisation. The `cst-drag` and `GIB-drag` cases test the geometric origin of the lateral-drag form factors. The `kstrength` and *Rothrock* cases test whether the simulated fast-ice strength/thickness behaviour is controlled by the ice-strength formulation rather than by lateral drag alone.

---

## Ice-strength thickness-exponent tests

The `kstrength-test01` and `kstrength-test02` experiments test a new Icepack/CICE ice-strength option:

$$
P = P^\ast h^{3/2} \exp[-C(1-A)].
$$

This is implemented as `kstrength = 2`. The tests are based on the `Cs-high` static lateral-drag configuration:

```text
boundary_condition = 'free_slip'
lateral_drag       = .true.
form factors       = coastline + grounded icebergs
form_func          = 'static'
Cs                 = 1.0e-3
Ktens              = 0.2
e_yieldcurve       = 1.5
e_plasticpot       = 1.5
Cstar              = 20
```

Two normalisations are:
| Experiment         | `kstrength` |  `Pstar` | Matching thickness | Role                                                                                                     |
|--------------------|------------:|---------:|-------------------:|----------------------------------------------------------------------------------------------------------|
| `kstrength-test01` |           2 | 19445.44 |                2 m | Conservative ($h^{3/2}$) test, close to standard Hibler strength near typical compact sea-ice thickness. |
| `kstrength-test02` |           2 | 15877.13 |                3 m | Weaker normalisation for 1--2 m ice; tests sensitivity to the chosen matching thickness.                 |

The ratio of the new ice strength to the standard linear-thickness Hibler strength is:

$$
\frac{\matrhm{P}_{\mathrm{new}}}{\matrhm{P}_{\mathrm{old}}} = \sqrt{\frac{\mathrm{hi}}{\mathrm{hi}_{\mathrm{match}}}}
$$

This makes the the 3 metre-normalised case systematically weaker than the 2 metre-normalised case for 1--2 m ice.

### September--November 2005 pilot comparison

The 90-day simulation starts from the 1 Sep 2005 restart and uses daily history output. Southern Ocean ice-mask-mean strength separates quickly between the two normalisations, while the bulk ice-area and extent diagnostics remain comparatively close over this short pilot window.

[Southern Ocean ice-mask-mean strength for the two kstrength pilot tests](figs/kstrength/metric_strength_mean_ice_mask.png)

A 30-day October map animation compares ice strength and sea-ice thickness over the Mawson Coast to Shackleton Ice Shelf sector:

[Download or view the October 2005 strength/thickness MP4](figs/kstrength/kstrength_strength_hi_200510.mp4)

---

## Analysis workflow

The recommended analysis workflow uses the `shuga` Python toolbox for:

1. CICE history conversion to Zarr;
2. binary-days fast-ice classification;
3. fast-ice area and persistence metrics;
4. comparison against the AF2020 Antarctic fast-ice dataset;
5. lateral-drag diagnostic analysis.

The preferred process diagnostics include:

- FIA time series;
- regional FIA;
- hit, miss, and false-alarm maps against AF2020;
- `ldwgt` phase-space plots;
- `ldphi` phase-space plots;
- `ldeps` and `ldspd` distributions;
- pack-ice collateral checks, including SIA, thickness, and speed.

The most important diagnostic distinction is:

```text
ldwgt -> which branch is active
ldphi -> what damping rate the momentum equation receives
```

A high `ldphi` value does not necessarily imply static locking, because the quadratic branch can also produce large damping at high speeds. `ldwgt` is the direct diagnostic of static/locking branch contribution.

---

## Known limitations

- This branch is a research-development branch, not an official CICE release.
- Form-factor generation depends on external high-resolution coastline and grounded-iceberg datasets.
- `eps_blend` is a resolved-grid strain-rate threshold, not a universal material constant.
- `blend_strain` sensitivity must be evaluated with diagnostics; circum-Antarctic FIA alone is insufficient.
- The current implementation focuses on C-grid dynamics.
- The final recommended parameter set remains experiment-dependent and should be evaluated regionally, not only circum-Antarctic.

## Relationship to upstream CICE

CICE is maintained by the CICE Consortium. This repository is a development branch/fork for testing free-slip and lateral-drag parameterisations.

Useful upstream resources:

- CICE Consortium repository: <https://github.com/CICE-Consortium/CICE>
- CICE documentation: <https://cice-consortium-cice.readthedocs.io/>
- CICE wiki: <https://github.com/CICE-Consortium/CICE/wiki>
- Icepack: <https://github.com/CICE-Consortium/Icepack>
- CICE Consortium forum: <https://xenforo.cgd.ucar.edu/cesm/forums/cice-consortium.146/>

For general CICE support, use the CICE Consortium forum rather than this development branch.

---

## License

This branch inherits the upstream CICE license and distribution policy. See:

```text
LICENSE.pdf
DistributionPolicy.pdf
```

where present in the repository.
