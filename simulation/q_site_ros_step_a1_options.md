# Step a1: Define the Q-site ROS reaction (options + stoichiometry)

This file expands step **a1** from `simulation_steps.md`: you must choose **what exact chemical event you mean by “Q-site ROS generation”** before you can design (or interpret) simulations.

Why this matters: different “ROS definitions” imply different **redox states**, **protonation states**, **spin states**, and therefore different **modeling methods** (MD-only vs QM/MM).

---

## Notation (used below)

### Quinone (Q) redox/protonation states

- `Q` = oxidized quinone (ubiquinone-like)
- `Q^.-` = anionic semiquinone radical (one-electron reduced)
- `QH^.` = neutral semiquinone radical (one-electron reduced + one proton)
- `QH2` = fully reduced quinol (two-electron reduced + two protons)

Related half-reactions (schematic):

```text
Q   + e-        -> Q^.-
Q^.- + H+       -> QH^.
QH^. + e- + H+  -> QH2
```

### Oxygen / ROS species

- `O2` = molecular oxygen (ground state is triplet; spin matters in QM)
- `O2^.-` = superoxide (radical anion)
- `HO2^.` = hydroperoxyl radical (neutral)
- `H2O2` = hydrogen peroxide

Acid-base equilibrium (aqueous; pKa ~ 4.8, environment-dependent):

```text
HO2^.  <->  H+ + O2^.-
```

At physiological pH, the dominant species in bulk water is usually `O2^.-` (but pocket pKa can shift).

---

## You choose 3 things (pick one from each category)

### 1) What do you call “ROS” at the Q-site?

1. **Superoxide formation**: `O2^.-` is produced near the Q-site.
2. **Hydroperoxyl formation**: `HO2^.` is produced (proton-coupled cases).
3. **Net peroxide output**: `H2O2` is formed (often via multiple steps, not a single elementary event).

### 2) Who is the electron (and/or H) donor in your model?

Common “Q-site-local” donors:

- `Q^.-` (anionic semiquinone)
- `QH^.` (neutral semiquinone)
- `QH2` (quinol)

More mechanistic donors (closer to the actual enzyme pathway, but harder):

- terminal **Fe-S cluster** near the Q-site (often called “N2” in Complex I literature) donating an electron to `O2` directly.

### 3) What mechanism do you assume?

- **Outer-sphere ET** (electron transfer; no bonds formed/broken)
- **PCET** (proton-coupled ET: `e-` and `H+` move in a coupled fashion)
- **HAT** (hydrogen-atom transfer: transfers `H^. = e- + H+` together)
- **Stepwise 2e- chemistry** (builds to `H2O2` by two sequential 1e- steps)

---

## Reaction options (equations + pros/cons)

Below is a menu of definitions you can adopt. Pick the one that matches your scientific question and compute budget.

### Option A — “O2 access only” (no chemistry; a prerequisite metric)

Definition (not a reaction):

```text
Measure P(O2 in reactive pose near Q-site) and residence times.
```

Pros:
- **MD-only** (GPU-friendly), easiest to do robustly with replicas.
- Strongly informative: if O2 cannot access the pocket, ROS chemistry is unlikely.

Cons:
- Not chemistry: cannot tell you if ET to O2 is favorable or fast.
- A mutation can change ET energetics without changing O2 access much (and vice versa).

When to pick:
- As the first screen for WT vs mutant differences.

---

### Option B — 1e- leak from anionic semiquinone to O2 (superoxide)

Reaction:

```text
Q^.-  +  O2   ->  Q  +  O2^.-
```

Pros:
- Clean, charge-balanced elementary ET step.
- Directly targets `O2^.-` formation (a primary ROS species).
- Often the most interpretable “Q-site ROS” definition.

Cons:
- Requires **open-shell QM/MM** (or another electronic method) to treat radicals/ET meaningfully.
- You must justify that `Q^.-` exists with sufficient population/lifetime (which the mutation may affect).

Good outputs to compare WT vs mutant:
- distribution of ET driving force (energy gap) across snapshots
- dependence on O2 position/orientation

---

### Option C — PCET: semiquinone donates e- while a proton is available (hydroperoxyl)

Reaction:

```text
Q^.-  +  O2  +  H+   ->   Q  +  HO2^.
```

Pros:
- Captures the fact that proton availability can change ROS speciation (`HO2^.` vs `O2^.-`).
- More chemically complete than Option B if the pocket provides proton donors (waters, sidechains).

Cons:
- Harder to model: you must specify **where the proton comes from** and include it in the reactive region.
- Protonation states in/near the pocket become a major uncertainty.

When to pick:
- If you believe the mutation changes pocket hydration/proton wires and you want that explicitly.

---

### Option D — HAT from quinol to O2 (hydroperoxyl + semiquinone)

Reaction (one plausible net form):

```text
QH2  +  O2   ->   QH^.  +  HO2^.
```

Pros:
- Represents a concerted transfer of `H^.` (electron + proton together), which is common in radical chemistry.
- Connects ROS directly to the presence of a reduced quinol state.

Cons:
- Requires modeling **both** the quinone redox/protonation state and radical chemistry (QM/MM).
- Whether `QH2` is present in the reactive configuration at the Q-site depends on catalytic turnover state assumptions.

When to pick:
- If your mechanistic hypothesis is “quinol-like state leaks H^./e- to O2”.

---

### Option E — “Net H2O2 formation” (two-electron product; usually multi-step)

Convenient net stoichiometry:

```text
O2  +  2 e-  +  2 H+   ->   H2O2
```

Or expressed via quinol (net, schematic):

```text
QH2  +  O2   ->   Q  +  H2O2
```

Pros:
- Closer to a common experimental readout (H2O2 production).

Cons:
- Usually **not** a single elementary step in a protein pocket; often proceeds via `O2^.-` / `HO2^.` intermediates
  and/or dismutation elsewhere.
- Hard to attribute to “Q-site” without a detailed mechanistic model and a pathway definition.

When to pick:
- If your end goal is to connect to H2O2 signals, but be prepared for a multi-step modeling campaign.

---

### Option F — Include downstream dismutation explicitly (often outside the pocket)

Reaction:

```text
2 O2^.-  +  2 H+   ->   H2O2  +  O2
```

Pros:
- Converts superoxide to peroxide explicitly (useful for mapping to assays).

Cons:
- Typically occurs in bulk solvent and/or is enzyme-catalyzed (SOD), not necessarily at the Q-site.
- Adds another layer of modeling that may not be mutation-specific to ND6.

When to pick:
- Usually as a post-processing/interpretation step, not as the primary “Q-site ROS” definition.

---

## Practical recommendation (starting point for WT vs ND6 M64V)

For a first, defensible computational study that’s feasible on a workstation:

1. Use **Option A** (O2 access/residence) with classical MD (replicas).
2. Combine with **Option B** (ET: `Q^.- + O2 -> Q + O2^.-`) on a **small set of clustered snapshots** via QM/MM or an embedded-cluster approach.

Why:
- It separates “can O2 get there?” from “is ET favorable?”, which are often the two dominant determinants of ROS propensity.
- It avoids having to fully model multi-step peroxide formation while still targeting the primary ROS event.

---

## What to record when you pick an option (so you can reproduce it later)

Write down:

- which option(s) you chose (A/B/C/D/E/F)
- which quinone state you assume (`Q`, `Q^.-`, `QH^.`, `QH2`)
- which oxygen product you treat as ROS (`O2^.-`, `HO2^.`, `H2O2`)
- whether proton donors are included and where the proton comes from (if applicable)
- whether you treat the donor as “the semiquinone” vs explicitly include the nearby Fe-S cluster as the donor

