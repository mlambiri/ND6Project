# Q-site ROS workflow (why MD engines aren’t enough + what to do instead)

That “important: none…” means:

Classical MD engines (OpenMM/GROMACS/NAMD) use **molecular mechanics force fields**: fixed bonds (or harmonic bonds), fixed atom types, and usually fixed partial charges. They do **not** have electrons, so they cannot natively model:

- **changing oxidation states** (Q ↔ semiquinone ↔ QH2),
- **electron transfer** from the Q-site redox chain to **O₂**,
- formation of **radicals** (e.g., superoxide O₂•−),
- the coupled electronic/protonic rearrangement that governs whether O₂ gets reduced.

They *can* still tell you the structural/diffusional prerequisites for Q-site ROS (Q pose, pocket hydration, O₂ access, gating motions), but the actual redox chemistry step needs an electronic-structure method (typically **QM/MM**).

---

## a) Steps that take you to the “important” step (QM/MM)

1) **Define the exact Q-site ROS reaction you mean**

- Example target: “electron leak from a semiquinone-like state to O₂ → O₂•−”.
- Decide whether you’re studying:
  - O₂ access only, or
  - semiquinone stabilization + O₂ access, or
  - the full ET/chemistry rate estimate.
- See `q_site_ros_step_a1_options.md` for a menu of reaction/stoichiometry choices (and pros/cons).

2) **Build two matched systems (WT and ND6 M64V)**

- Same membrane composition/size, water, ions, temperature, constraints, etc.
- Keep the Q-site ligand present (your PDB has `8Q1`; ensure it’s actually in the pocket and stays there).

3) **Decide how you will represent redox states in classical MD**

Because MD can’t change redox state by itself, you typically run *separate* classical systems with different parameters for:

- oxidized Q,
- reduced QH2,
- semiquinone-like Q•−/QH• (harder; requires a defined charge/protonation model).

4) **Equilibrate**

- Minimization → restrained equilibration → unrestrained equilibration.
- Verify Q-site stability (no ligand drifting out, no broken coordination, etc.).

5) **Production MD (multiple replicas)**

- Collect long trajectories so you can compare WT vs mutant statistically (not just one run).

6) **O₂ access sampling**

- Add many O₂ molecules (“O₂ cosolvent”) and run MD to map where O₂ accumulates/resides.
- Output: O₂ density/residence times near the Q pocket, and how often O₂ is in a “reactive pose” near the quinone headgroup.

7) **Select “reactive snapshots”**

- Frames where:
  - O₂ is within a cutoff of the quinone reactive region, and
  - Q-site geometry/hydration matches the state you want to test.
- Cluster them so you don’t QM/MM 1000 near-duplicates.

8) **Now you’re at the “important” step**

You have representative WT and mutant snapshots with O₂ in/near the Q pocket and a chosen redox/protonation model to evaluate.

---

## b) The “important” step itself: what you actually do in QM/MM (detailed)

### 1) Choose the *specific* QM target calculation

You usually break “ROS generation” into components you can compute:

- **O₂ binding/placement** in the pocket (thermodynamics)
- **Electron transfer (ET) feasibility**: is electron-on-O₂ energetically favorable in that environment?
- **ET rate** (kinetics): depends on coupling + reorganization + driving force
- **Proton-coupled effects** (if you include QH• vs Q•−, nearby proton donors, waters)

A very common practical target is: estimate how “easy” it is for O₂ to become **O₂•−** given the pocket electrostatics and quinone state.

### 2) Define the QM region (keep it small, but chemically complete)

For Q-site ROS, a typical QM region includes:

- the **quinone** headgroup (your `8Q1`, or whatever Q species you model),
- the **O₂** molecule (the one in the reactive pose),
- **key nearby residues** that H-bond to / polarize the quinone (sidechains only if possible),
- **a few nearby waters** if they directly contact the quinone/O₂.

You may *optionally* include the nearest **Fe–S cluster** only if your question requires explicitly modeling the donor redox center. For “electron leak from semiquinone to O₂”, you can often treat the semiquinone as the donor without including a Fe–S cluster in QM.

### 3) Define the MM environment and boundary handling

Two common approaches:

**(A) True QM/MM**

- Whole protein/membrane/water stays as MM.
- QM region sits inside and “feels” MM point charges (electrostatic embedding).
- You cut covalent bonds at the QM/MM boundary using link atoms/caps.

**(B) Embedded cluster (cheaper, common on workstations)**

- Carve a sphere around the Q-site (e.g., 15–25 Å).
- Freeze outer atoms; keep inner atoms flexible.
- Include point charges from the outer region (or a dielectric continuum) to mimic the environment.

(A) is more rigorous; (B) is often much more tractable on a single workstation.

### 4) Assign QM charge and spin (critical for radicals)

You must set:

- total **charge** of the QM region (depends on Q state + residue protonation),
- **spin multiplicity** (open-shell) for:
  - semiquinone state,
  - O₂ (triplet) and superoxide (doublet),
  - any radical intermediates.

This is one of the biggest “gotchas”: wrong charge/spin → meaningless ROS energetics.

### 5) Pick a QM method appropriate for open-shell ET

Typically:

- **DFT** (density functional theory) with dispersion corrections for pocket interactions.
- Open-shell settings (unrestricted / broken-symmetry where needed).
- A basis set that’s not too small (otherwise ET energetics are unreliable).

You’ll likely do:

- local geometry optimization of the QM region with the environment fixed/partly flexible,
- single-point energy evaluations across states/geometries.

### 6) Compute ET/ROS metrics from an ensemble (not one snapshot)

For each representative snapshot (WT and mutant), compute something like:

**(i) Driving force (ΔG-like) for electron on O₂**

- Compare energies of two “states”:
  - electron on quinone (donor state),
  - electron on O₂ (acceptor state → superoxide-like).
- Methods vary (constrained DFT, diabatic state approximations, etc.), but goal is consistent WT vs mutant comparison.

**(ii) Reorganization energy (λ)**

- How much the environment must reorganize for ET to occur.
- Estimated by vertical energy-gap calculations across geometries.

**(iii) Electronic coupling (V)**

- How strongly donor and acceptor are electronically connected (geometry dependent).

Then you can estimate an **ET rate** using a Marcus-type model (rate depends on ΔG, λ, V), and compare WT vs M64V distributions.

### 7) Combine chemistry + access to get a ROS “propensity”

A practical comparison metric is:

ROS propensity ≈ (how often O₂ is in the right place from MD) × (ET rate / ET favorability from QM/MM)

Do this for many snapshots/replicas and compare WT vs mutant with uncertainty bars.

### 8) Sensitivity checks (you must do these)

ROS conclusions can flip if you change:

- quinone redox/protonation model,
- key residue protonation states near the pocket,
- which waters are included,
- which snapshots you sampled.

So you run a small matrix of alternatives to see if the WT vs mutant *direction* is robust.

---

If you want, tell me:

- which MD engine you’re leaning toward (OpenMM vs GROMACS/WSL2),
- whether you want to keep `8Q1` exactly as-is,
- and whether you want the “important step” to be **true QM/MM** or an **embedded cluster** on your PC,

and I’ll outline a concrete, minimal workflow for Q-site ROS that’s actually feasible on a workstation.
