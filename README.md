
---

## 3. Suggested top-level `README.md` for `CICE_free-slip`

I would keep the repo README structured around the physical progression rather than code details first.

```markdown
# CICE free-slip and lateral-drag development branch

This branch develops and evaluates free-slip coastal boundary conditions and lateral-drag parameterisations for representing Antarctic landfast sea ice in CICE.

The central physical problem is that Antarctic landfast sea ice is not simply slow pack ice. It forms where coastal geometry, ice shelves, islands, grounded icebergs, internal ice stress, and local deformation interact to immobilise otherwise mobile sea ice. This branch tests whether high-resolution coastline and grounded-iceberg form factors can provide that geometric anchoring in CICE without imposing a no-slip boundary condition that unrealistically suppresses coastal pack-ice motion.

## 1. Boundary-condition motivation

In an idealised no-slip configuration, ice velocity is forced to vanish at solid boundaries. This strongly constrains coastal motion and can promote artificial immobilisation near land.

A free-slip boundary condition removes the tangential no-slip constraint. This allows ice to slide along coastlines and grounded-iceberg boundaries, which is more appropriate when the goal is to let landfast ice emerge from the resolved momentum balance rather than be imposed directly by the boundary condition.

However, free-slip alone removes an important source of unresolved lateral resistance. This motivates an explicit lateral-drag term.

## 2. C-grid free-slip implementation

The C-grid implementation uses staggered E- and N-face velocity components. The free-slip boundary treatment is applied in the strain-rate and stress-divergence calculations used by the EVP dynamics.

The relevant CICE dynamics files are:

- `cicecore/cicedyn/dynamics/ice_dyn_evp.F90`
- `cicecore/cicedyn/dynamics/ice_dyn_shared.F90`

The EVP driver in `ice_dyn_evp.F90` computes strain rates, internal stress, stress divergence, and then updates the E- and N-face velocity components through `stepu_C` and `stepv_C` in `ice_dyn_shared.F90`.

## 3. Lateral-drag form factor

Lateral drag is applied where high-resolution coastline and grounded-iceberg geometry indicate unresolved lateral contact or anchoring potential.

The static geometric factor is represented through E- and N-face form factors:

- `F2E`
- `F2N`

The effective lateral-drag mass/form factor is:

```text
Ku = M_x F2
```

where M_x is the velocity-grid ice/snow mass and F2 is the local form factor. This makes the lateral stress scale with both the available ice mass and the local geometric exposure to coastline or grounded-iceberg anchoring.

The lateral-drag stress is applied as a dissipative momentum sink:

tau_LD = - Ku phi u

where phi is the selected lateral-drag form function.

4. Lateral-drag form functions

The implementation supports several form functions for evaluating the damping rate phi.

static
phi_static = Cs / (|u| + u0)

This is the Liu-style reference branch. It produces a strong low-speed locking tendency while remaining regularised by u0.

quadratic
phi_quad = Cq |u|

This is a mobile, velocity-scaled branch. It becomes weak at low speed and stronger for faster moving ice.

linear
phi_linear = C_L

This is a Rayleigh-style damping branch.

blend_strain

The corrected blend_strain formulation combines a static locking branch and a quadratic mobile branch using both speed and strain-rate gates.

w_eps  = 1 / (1 + (eps_eff / eps_blend)^p)
w_u    = 1 / (1 + (|u| / u_blend)^p)
w_lock = w_eps w_u

phi_blend = w_lock phi_static + (1 - w_lock) phi_quad

Thus:

low speed + low strain-rate  -> static/locking branch
high speed or high strain    -> quadratic/mobile branch

This formulation is designed to permit landfast ice to emerge where ice is slow, coherent, and geometrically anchored, while allowing mobile or deforming pack ice to remain dynamically mobile.

5. Strain-rate normalisation

The C-grid strain-rate diagnostic deltaU is a deformation invariant multiplied by U-cell area. Therefore the effective strain rate used by blend_strain is area-normalised before comparison with eps_blend.

For E-face velocity points:

eps_eff = (deltaU(i,j) + deltaU(i,j-1)) / (uarea(i,j) + uarea(i,j-1))

For N-face velocity points:

eps_eff = (deltaU(i,j) + deltaU(i-1,j)) / (uarea(i,j) + uarea(i-1,j))

This produces a strain-rate scale in s^-1.

6. Diagnostics

The branch includes lateral-drag diagnostics intended to support process-based evaluation rather than tuning against fast-ice area alone.

Key diagnostics include:

KuxE, KuyE, KuxN, KuyN: lateral-drag stress components
F2E, F2N: geometric form factors
ldphiE, ldphiN: realised damping rate
ldwgtE, ldwgtN: static/locking branch weight
ldepsE, ldepsN: effective strain rate used by blend_strain
ldspdE, ldspdN: speed used by the form function
ldpstatE, ldpstatN: static branch damping-rate diagnostic
ldpquadE, ldpquadN: quadratic branch damping-rate diagnostic
ldplinE, ldplinN: linear branch damping-rate diagnostic

These diagnostics are written only on the final EVP subcycle to avoid excessive memory traffic within the dynamics loop.

7. Recommended interpretation

The intended public-facing CICE options are:

static: simple Liu-style lateral-drag reference
blend_strain: physically gated extension for Antarctic landfast ice

The quadratic and linear branches are retained as diagnostic comparators and sensitivity tools.

8. Current experiment focus

The current Sep-Dec 1994 experiment set tests whether corrected blend_strain can improve Antarctic fast-ice formation, persistence, and retreat relative to heavier-handed static, quadratic, and linear form functions.

The primary corrected blend_strain hypothesis is:

u_blend   = 5.0e-4 m s^-1
eps_blend = 5.0e-8 s^-1
blend_exp = 3.0

For representative Antarctic fast-ice cell scales of approximately 12 km, this aligns the speed and strain-rate gates across the fast-ice velocity range.


---

Present status: implementing the `blend_exp_int` change, smoke-test `blend_exp = 3.0` for a few days, then launch `LD-blend-base`, `smooth`, and `sharp` before the stricter/permissive cases.

---

<!--- [![Travis-CI](https://travis-ci.org/CICE-Consortium/CICE.svg?branch=main)](https://travis-ci.org/CICE-Consortium/CICE) --->
[![GHActions](https://github.com/CICE-Consortium/CICE/workflows/GHActions/badge.svg)](https://github.com/CICE-Consortium/CICE/actions)
[![Documentation Status](https://readthedocs.org/projects/cice-consortium-cice/badge/?version=main)](http://cice-consortium-cice.readthedocs.io/en/main/?badge=main)
[![lcov](https://img.shields.io/endpoint?url=https://apcraig.github.io/coverage.json)](https://apcraig.github.io)

<!--- [![codecov](https://codecov.io/gh/apcraig/Test_CICE_Icepack/branch/master/graph/badge.svg)](https://codecov.io/gh/apcraig/Test_CICE_Icepack) --->

## The CICE Consortium sea-ice model
CICE is a computationally efficient model for simulating the growth, melting, and movement of polar sea ice. Designed as one component of coupled atmosphere-ocean-land-ice global coupled models, today’s CICE model is the outcome of more than two decades of community collaboration in building a sea ice model suitable for multiple uses including process studies, operational forecasting, and Earth system simulation.


This repository contains the files and code needed to run the CICE sea ice numerical model starting with version 6. CICE is maintained by the CICE Consortium. 
Versions prior to v6 are found in the [CICE-svn-trunk repository](https://github.com/CICE-Consortium/CICE-svn-trunk).

CICE consists of a top level driver and dynamical core plus the [Icepack][icepack] column physics code], which is included in CICE as a Git submodule.  Because Icepack is a submodule of CICE, Icepack and CICE development are handled independently with respect to the GitHub repositories even though development and testing may be done together.  

[icepack]: https://github.com/CICE-Consortium/Icepack

The first point of contact with the CICE Consortium is the Consortium Community [Forum][forum]. 
This forum is monitored by Consortium members and also opened to the whole community.
Please do not use our issue tracker for general support questions.

[forum]: https://xenforo.cgd.ucar.edu/cesm/forums/cice-consortium.146/

If you expect to make any changes to the code, we recommend that you first fork both the CICE and Icepack repositories. 
In order to incorporate your developments into the Consortium code it is imperative you follow the guidance for Pull Requests and requisite testing.
Head over to our [Contributing][contributing] guide to learn more about how you can help improve CICE.

[contributing]: https://github.com/CICE-Consortium/About-Us/wiki/Contributing

## Useful links
* **CICE wiki**: https://github.com/CICE-Consortium/CICE/wiki

   Information about the CICE model

* **CICE Release Table**: https://github.com/CICE-Consortium/CICE/wiki/CICE-Release-Table

   Numbered CICE releases since version 6 with associated documentation and DOIs. 
   
* **Consortium Community Forum**: https://xenforo.cgd.ucar.edu/cesm/forums/cice-consortium.146/

   First point of contact for discussing model development including bugs, diagnostics, and future directions.   

* **Resource Index**: https://github.com/CICE-Consortium/About-Us/wiki/Resource-Index

   List of resources for information about the Consortium and its repositories as well as model documentation, testing, and development.

## License
See our [License](LICENSE.pdf) and [Distribution Policy](DistributionPolicy.pdf).
