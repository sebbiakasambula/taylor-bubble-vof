# Solver configuration in VOF simulations of Taylor bubbles on film-resolving stretched meshes

Sebbi Adam Kasambula

A Taylor bubble is a long gas bubble rising through liquid in a vertical pipe, with a thin film of liquid falling between the gas and the wall. To resolve that film in a volume of fluid simulation the mesh has to be stretched towards the wall, so the cells at the interface are long and thin. I ran three rising bubbles (Eo = 300, 100 and 20 at Nf = 1000, 50 and 500) and a static bubble in OpenFOAM v2412 and in Ansys Fluent 2024 R1, on the same axisymmetric meshes, with the same fluid properties and time step, and asked how much of the difference between the two codes comes from their solver settings.

On identical inputs Fr in Fluent differed from OpenFOAM by +2.5%, +20.8% and -10.7%, all inside the scatter of the correlation of Viana et al. At Eo = 100 I changed the differing settings one at a time and compared with the sharp-interface simulations of Abubakar and Matar and the nearest measurement, which the OpenFOAM baseline underestimated by 13.7% and 18.2%. Solving the momentum equation before the pressure correction raised Fr by 5.2%, close to the 4.1% from halving the time step. Geometric advection with isoAdvector then brought Fr and the film to within 1.0% and 0.4% of the sharp-interface values. Changing the surface tension term of OpenFOAM to the density-weighted form that Fluent uses lowered Fr by 29.5%. In every run at Eo = 100, Fr lay within 3.3% of the value that laminar film theory gives for its film. The same settings lowered Fr at Eo = 300 by 3.8% and left it within the fit error at Eo = 20. The paper has the method, the grid studies and the uncertainty budget.

## Contents

- `paper/`, the paper (`paper.tex`, `paper.pdf`, `references.bib`, `highlights.txt`) and its figures as EPS files in `paper/figures/`.
- `solver/interFoamRhoST/`, interFoam v2412 with one switch added, `rhoWeightedSurfaceTension` in the PIMPLE dictionary. When it is on, the surface tension force on each face is multiplied by the face density over the mean density, in both the momentum and the pressure equation. With the switch off it matched the standard interFoam over 64 time steps.
- `scripts/`, the scripts that write and run the cases, the first analysis chain (arrays, per-snapshot tables, comparison tables) and the field figures. `scripts/paraview/` holds the ParaView scripts.
- `audit/`, the scripts behind the numbers of the paper as it stands. Each `set*.py` recomputes one group of numbers and writes its result next to itself as `set*_out.md`. `of_metrics.py` and `fl_metrics.py` make the per-snapshot tables in `m/` from the OpenFOAM run folders and the Fluent exports, and `intact_check.py` the tables in `intact/` (whether the bubble stayed in one piece). `run_new.py` makes the tables of the Fluent tests in `flcopy/`, which `set8_fluent_out.md` sums up. The other folders hold stored inputs: `a5/` arrays of the geometric advection runs and the Eo = 20 baseline, `cmp/` the comparison tables, `axisgas/` the gas on the axis, `hysing/` the rising bubble benchmark of Hysing et al., `xb9/` the check of the density weighting at Eo = 20, and `profiles.npz` the film velocity profiles. The `fig_*.py` scripts here make six of the paper figures.
- `abu/`, the digitising of two figures of Abubakar and Matar (Fig. 6a, Fr against Eo, and Fig. 7b, the film) and the measured Fr from the appendix of Viana et al., typed from their tables. The scripts read crops of the published figures, which are not included here. The digitised values are in the `*_out.txt` files and `RESULTS_abubakar_T2.md`.
- `data/`, the stored arrays (`.npy`) and per-snapshot and comparison tables (`.csv`) of the first chain, so the tables and several figures can be remade without the run folders.

There is no simulation data in this repository. The OpenFOAM run folders are large and the scripts rebuild them. The Fluent exports are not included either.

## Requirements

- OpenFOAM v2412 through the official ESI Docker image `opencfd/openfoam-run:2412`. The scripts that write and run cases call Docker.
- Python 3 with numpy, scipy and matplotlib, h5py for the Fluent `.dat.h5` files read by `audit/run_new.py`, and Pillow for the digitising in `abu/`.
- ParaView 6.2 (`pvbatch`) for the field figures.
- Ansys Fluent 2024 R1 for the Fluent side of the comparison.
- Time. The medium mesh at Eo = 300 took 5.7 h to 39.4 s on 2 processes of a laptop with an Intel i3-8145U, and the run at Eo = 20 took 25.1 h.

## Paths

All paths are read in `scripts/config.py` from three environment variables.

| variable | content | default |
|---|---|---|
| `P2_RUNS` | the OpenFOAM run folders, one per case | `~/cfdml/runs_p2` |
| `P2_FLUENT` | the Fluent exports | `../Fluent_Data` next to the repository |
| `P2_WORK` | arrays, per-snapshot tables, VTK files, render strips and figure outputs | `./work` |

Run every script from the repository root as `python scripts/<name>.py` or `python audit/<name>.py`, and the ParaView scripts as `pvbatch scripts/paraview/<name>.py`. No script writes into `paper/figures`. The figure scripts write into `P2_WORK`, and the finished figures are copied into `paper/figures` by hand. Figures written as PDF were turned into EPS with `pdftops -eps`.

## How to reproduce

Run the steps in this order, since later steps read the output of earlier ones.

1. Cases and runs
   `python scripts/vv_setup.py` writes and meshes the static bubble and the Eo = 300 cases (coarse, medium and fine meshes, the halved time step and the long domain). `python scripts/mk_t.py` writes the Eo = 100 and Eo = 20 cases. `python scripts/diag3.py` writes and runs the uniform 96 by 12 static mesh. Each case is then run with `python scripts/run_case.py <case folder>`, which uses 1 process up to 20 000 cells and 2 above.
2. Runs that are copies of a base case
   These have no script of their own, since each one changes a few lines of a case from step 1. The tests at Eo = 100 and Eo = 20 run to `endTime 8.4` with `writeInterval 0.2`.
   - Square mesh at Eo = 300: `write_case` from `tb_case.py` with `Eo=300, Nf=1000, rho_r=1000, nx=1920, nr=48, grading=1.0`.
   - Uniform 384 by 48 static mesh: the case of `diag3.py` with 384 axial and 48 radial cells.
   - Momentum predictor: set `momentumPredictor yes;` in PIMPLE, with a `"U.*"` entry in the solvers dictionary. Eo = 100 runs to 8.4 s, the medium mesh at Eo = 300 to 10 s and the square mesh to 15 s.
   - Three outer loops at Eo = 100: the momentum predictor run with `nOuterCorrectors 3;`.
   - Geometric advection: the momentum predictor run with `application interIsoFoam;` and, in the `alpha.water` dictionary of `fvSolution`, `reconstructionScheme isoAlpha;` (or `plicRDF;`), `nAlphaBounds 3;`, `isoFaceTol 1e-8;`, `surfCellTol 1e-8;`, `snapTol 1e-12;` and `clip true;`. At Eo = 300 it runs to 10 s with `writeInterval 0.25`.
   - Curvature smoothing at Eo = 100: `nAlphaSmoothCurvature 2;` in the `alpha.water` dictionary.
   - Half the time step at Eo = 100: `maxDeltaT 0.0004262375` in `controlDict`.
   - Wider wedge at Eo = 100: the wall vertices of `blockMeshDict` at a wedge of 10 degrees instead of 5.
   - Coarse and fine meshes at Eo = 100: the Eo = 100 case on the 343 by 69 and the 672 by 134 blocks of the Eo = 300 grid study, with `maxDeltaT 0.00140711` and `0.000515677`. The geometric advection and momentum predictor tests on these meshes change them as above.
   - Density weighting at Eo = 100 and Eo = 20: set `application interFoamRhoST` in `controlDict` and add `rhoWeightedSurfaceTension yes;` to the PIMPLE dictionary in `fvSolution`. The static bubble with the weight is the static case with the same two dictionary changes and its own end time of 3.95 s.
   - Cell centres and volumes: the analysis reads `0/C` and `0/V`, written with `postProcess -func writeCellCentres` and `postProcess -func writeCellVolumes` inside the container.
   - The solver: build `solver/interFoamRhoST` with `wmake` inside the OpenFOAM v2412 container. On Windows, copy the source into the container's `/tmp` first, since the mounted folder is case-insensitive and `CorrectPhi.H` clashes with `correctPhi.H`.
3. Fluent
   The Fluent runs use a mesh that matches the OpenFOAM medium mesh to within 1%, and the settings in the paper's settings table, with 20 iterations per time step. Each run exports, every few hundred steps, an ASCII file of cell centre coordinates, pressure, axial and radial velocity, cell volume and gas volume fraction, plus report files of the gas centroid and gas volume. `f_t1.py`, `f_t1sq2.py`, `f_t2.py` and `f_t3.py` read these exports.
4. Arrays and tables
   `python scripts/t_t23.py <case> ...` for the OpenFOAM arrays, the four `f_t*.py` scripts for the Fluent arrays, then `run_of.py` and `run_fl.py` for the per-snapshot tables and `cmp_table.py` for the comparison tables. The weighting tests are read with `t2rho.py`, `t2rho_derive.py`, `t3rho_derive.py`, `b11_srho.py` and `b11_report.py`. The grid study numbers come from `v_film.py` and `v_grid.py`, the mesh table from `meshfacts.py` and the gas drift from `drift_all.py`. To start from the stored arrays and tables instead, copy the files in `data/` into the work folder.
5. The numbers of the paper
   `python audit/of_metrics.py`, `audit/fl_metrics.py` and `audit/intact_check.py` make the tables in `audit/m/` and `audit/intact/`, then each `audit/set*.py` writes its `set*_out.md`. `set1_static.py`, `set9_t3tests.py`, `set16_m2f.py`, `set17_film_theory.py`, `set18_budget.py` and `hysing/hysing_post.py` run from the stored files alone. The other `set*.py` also read the run folders or the Fluent exports. `set17_profiles.py` writes `profiles.npz`, the film profiles of the film figure, and `set18_budget.py` gives the uncertainty budget.
6. Figures

| figure | scripts, in order |
|---|---|
| `fig_t2_split` (the bubble split on the fine mesh) | `prep_split.py`, `paraview/render_split.py`, `assemble_split.py` |
| `fig_fr_sq` (Fr over 5 s windows) | `fig_fr_split.py` |
| `fig_film_mesh` (cells across the film) | `prep_film_mesh.py`, `paraview/render_film_mesh.py`, `fig_film_mesh.py` |
| `fig_s_fields`, `fig_t3_fields` (static bubble and Eo = 20 fields) | `tovtk.py`, `paraview/render_flat.py`, `assemble_s_flat.py`, `assemble_t3_flat.py` |
| `fig_t1_fields` (Eo = 300 fields on both meshes) | `tovtk_t1.py`, `paraview/render_t1_flat.py`, `fig_t1_fields_flat.py` |
| `fig_t1mesh` (mesh shape at Eo = 300) | `audit/fig_t1mesh_bars.py` |
| `fig_t1_outline` (Eo = 300 interface outlines) | `fig_t1_outline.py` |
| `fig_cases`, `fig_cases2` (the code comparison) | `audit/fig_cmp2.py` |
| `fig_shapes` (bubble shapes) | `fig_shapes_compact2.py` |
| `fig_t2_tests` (the settings at Eo = 100) | `audit/fig_t2_ladder2.py` |
| `fig_t2_film` (the film against laminar film theory) | `audit/fig_t2_film.py` |
| `fig_transfer` (the same settings at all three cases) | `audit/fig_transfer.py` |

The render scripts draw the field without lighting, so the colours in the strips match the colour bar. `fig_t3rho_flat.py` and `fig_cases_compact.py` made figures of an earlier version of the paper.

## Licence

MIT, see `LICENSE`.
