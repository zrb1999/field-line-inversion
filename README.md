# Code for "Remote sensing of the inner-magnetospheric magnetic field from field-line geometry: a 2D Hamiltonian inversion method"

MATLAB implementation of the two-dimensional axisymmetric Hamiltonian level-set
inversion described in the manuscript submitted to *Geophysical Research Letters*.
The archive generates the synthetic target fields, performs the anchor-differenced
inversion, evaluates the field metrics, runs the sensitivity experiments, and
reproduces Figures 1--4, Figures S1--S2, and Tables S1--S7 of the paper and its
Supporting Information.

## Contents

| Path | Content |
|---|---|
| `paper1_2d_hamiltonian/` | The scripts and saved MAT results needed to reproduce the paper |
| `paper1_2d_hamiltonian/manuscript/figures/` | Figures produced by the `make_*.m` scripts (PDF and PNG) |
| `3d/igrf13_gh.m` | IGRF-13 coefficient reader required by `igrf_m0_coeffs.m` |

The `3d/` folder is a sibling of `paper1_2d_hamiltonian/` because
`igrf_m0_coeffs.m` adds it to the MATLAB path with
`addpath(fullfile(thisdir,'..','3d'))`. Keep this layout when running the code.

## Requirements

- MATLAB R2020a or later (the figure scripts use `exportgraphics`).
- Aerospace Toolbox: `igrf13_gh.m` loads the IGRF-13 coefficient table
  `aeroigrf13dataTable.mat` from the toolbox.
- No other toolboxes or third-party packages.

## Configuration (defaults in the scripts)

- Internal baseline: IGRF-13, zonal (`m=0`), 2020.0, degree 10.
- Target fields: 20 realizations, each the baseline plus a five-year difference
  of the historical IGRF-13 coefficient series (start years 1920--2015 in
  five-year steps, degrees 1--10).
- External field: Gaussian-weighted superposition of nine circular loops
  (`a_0 = 6 R_E`, `sigma_a = 0.8 R_E`, `B_c = -50 nT`), represented in the
  inversion by four external harmonics.
- Constellation: five L-shells at 1.2, 1.5, 1.9, 2.4, and 2.7 `R_E`,
  100 points per line, 0.1 km radial Gaussian path noise.
- Inversion: 10 internal perturbations plus 4 external moments (14 unknowns),
  solved by column-scaled damped least squares (`lambda = 1e-4`); endpoint
  `|B|` enters through Gauss--Newton iterations with `sigma_B = 1 nT`.
- Evaluation grid: `30 x 30` over `r = 1.05--3.10 R_E` and
  `theta = 10--90` degrees.

## Reproducing the paper

Run the scripts from inside `paper1_2d_hamiltonian/` in this order:

| Script | Paper element |
|---|---|
| `run_g1free_nline_grid.m` | `N=1`--5 field metrics (Table S2, Figures 2 and S2) |
| `run_g1free_nline_grid_Ne5.m` | `N_e=5` variant of the `N` sweep (Table S5) |
| `run_grid_field_cases.m` | Ring / no-ring control on the fixed grid (Table S3) |
| `run_sensitivity_structure.m` | External truncation and leave-one-out coverage (Table S5) |
| `run_sensitivity_obs.m` | Path noise, endpoint error, ring amplitude (Table S6, Figure 4) |
| `run_sensitivity_extended.m` | 2--5 km path noise and the quiet-to-storm sweep (Table S6) |
| `run_tracing_accuracy_mitigation.m` | `N=5/10/15` and points-per-line sweep at 5 km (Table S7) |
| `run_noise_geometry_test.m` | Radial, tangential, isotropic, and correlated 5 km errors (Table S7) |
| `run_model_complexity_test.m` | Internal-order sweep at 5 km (Table S7) |
| `make_all_figures.m` | All main-text and supporting figures |

Individual figures can be regenerated with `make_fig1_concept.m`,
`make_fig2_nline_field.m`, `make_fig3_gridmaps.m`, `make_fig4_sensitivity.m`,
`make_figS1_svd_g1free.m`, and `make_figS2_local_global.m`.

Core numerical routines:

| File | Purpose |
|---|---|
| `flux_basis.m`, `flux_basis_ext.m` | Internal and external zonal-harmonic basis functions and derivatives |
| `invert_flux_perturbation_full.m` | 14-parameter linear (path-only) inversion |
| `invert_flux_perturbation_B_iter_full.m` | Iterative Gauss--Newton inversion with endpoint `|B|` |
| `igrf_m0_coeffs.m` | IGRF-13 zonal Gauss coefficients `g_10`--`g_100` |
| `paper_style.m` | Shared figure style |

The archive contains only the scripts and MAT files needed to reproduce the
published results; development and diagnostic variants used during the study are
not included.

## Reproducibility

- The ensemble scripts fix the random seed (`rng(8000+iset)` and related
  seeds), so the results are deterministic.
- The saved `*.mat` files let the figures be regenerated without repeating the
  inversions. Running the `run_*.m` scripts overwrites them.

## License

Released under the MIT License; see `LICENSE`.

## Citation

If you use this code, please cite the software and the paper:

> Zhu, R., & Zhang, X. (2026). *Code for "Remote sensing of the
> inner-magnetospheric magnetic field from field-line geometry: a 2D
> Hamiltonian inversion method"* (Version 1.0.0) [Computer software]. Zenodo.
> https://doi.org/10.5281/zenodo.XXXXXXX
