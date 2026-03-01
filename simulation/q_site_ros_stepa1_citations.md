# Revised step content for Q-site ROS options

## Executive summary

A scientifically grounded “Q-site ROS” specification has to distinguish (at minimum) between **complex I site I\_Q** and **complex III site III\_Qo**, because both involve quinone/semiquinone intermediates yet differ in **thermodynamic drivers, inhibitor responses, and topology of superoxide release**. Peer‑reviewed work supports (i) **two major ROS‑generating loci in complex I** (a flavin‑linked site plus a Q‑site associated with semiquinone chemistry), (ii) **Qo‑site semiquinone involvement in complex III ROS**, especially under antimycin‑perturbed Q‑cycle conditions, and (iii) the need to treat common ROS assays (e.g., Amplex Red/UltraRed, hydroethidine/MitoSOX fluorescence) as **artifact‑prone unless properly controlled**. citeturn14view0turn20view0turn15view0turn24view0turn22view0turn12view0turn13search0turn13search2turn9view0turn10view0turn27view0turn18search1  

Because the current tool session has **no access to your referenced repository files** (no connected sources and no uploaded files were available), I cannot quote or diff the *actual* contents of `q_site_ros_step_a1_options.md` or its reference in `simulation_steps.md`. The replacement text below is therefore a **complete, paste‑ready Step A1 rewrite** that is internally consistent and densely referenced to primary literature and authoritative reviews, but the “original vs revised” table is necessarily **a best‑effort mapping of common failure modes in Q‑site ROS write‑ups** rather than a literal line‑by‑line comparison. citeturn22view0turn18search1  

## Evidence base for scientifically accurate Q-site ROS options

Mitochondrial ROS discussions are easiest to make precise if you explicitly define the **proximal species** (usually superoxide, O₂•⁻) and the **production site** within the respiratory chain. A substantial experimental literature supports a **two‑site model for complex I‑linked ROS**: one component is consistent with a flavin/NADH‑linked locus (often called site I\_F), and another depends on **Q‑pool redox state and protonmotive force** and is consistent with a **semiquinone at the Q‑binding site** (site I\_Q). citeturn14view0turn20view0turn22view0  

For complex III, mechanistic and kinetic studies strongly link high superoxide production (especially when the Q‑cycle is perturbed) to the **outer quinone (Qo) site** and a **Qo‑site semiquinone intermediate**. In isolated mitochondria treated with antimycin A, superoxide production shows a **bell‑shaped dependence** on the reduction state of cytochrome b and is consistent with superoxide originating from a Qo‑site semiquinone. citeturn15view0turn24view0turn23view0turn18search1  

The **topology** (which side of the inner membrane superoxide appears on) is also not optional; it changes what can be detoxified locally and what can signal elsewhere. Work in intact mitochondria indicates **complex I superoxide is released to the matrix**, whereas complex III‑derived superoxide can be detected as release toward the “cytoplasmic” side (intermembrane space) and—under some conditions—evidence also supports delivery toward the matrix; there is also an active debate about the quantitative partitioning, especially once cristae geometry and indirect assay limitations are considered. citeturn21view0turn17search0turn18search1turn18search2  

Finally, any “options” document that recommends measurement approaches must reflect that many popular ROS readouts are **not faithfully selective** without additional chemical analytics and bioenergetic controls. For example, Amplex Red can be enzymatically converted to resorufin **without H₂O₂** in some biological preparations, and hydroethidine/MitoSOX red fluorescence is **not a reliable proxy** for superoxide unless the superoxide‑specific products are separated/quantified (e.g., by HPLC), with additional care because probe loading itself can perturb mitochondrial function. citeturn12view0turn13search2turn13search0turn22view0  

## Unspecified parameters and assumptions that must be made explicit

Across the peer‑reviewed literature, the rates and even the *dominant site* of mitochondrial ROS production can change with metabolic context because the controlling variables include **Δp (protonmotive force), NADH/NAD⁺, QH₂/Q, and O₂ concentration**. If a Step A1 options document does not specify these, it cannot make defensible predictions or traceable simulation parameters. citeturn22view0turn14view0turn15view0turn27view0  

Similarly, inhibitor‑based logic must be treated as **context dependent**: Q‑site inhibitors of complex I can induce high superoxide production under particular conditions, and complex III inhibitors can either elevate or suppress H₂O₂/superoxide depending on where in the Q‑cycle they block and what upstream redox consequences they impose. citeturn20view0turn28view0turn15view0turn9view0  

For a simulation‑oriented “options” step, these omissions should be handled explicitly as “unspecified” with defended default choices rather than silently assumed. citeturn22view0turn14view0  

## Replacement text formatted for direct insertion

```md
# Q-site ROS options

This document defines the scientifically supported options for **Step A1** (Q-site-derived reactive oxygen species, “Q-site ROS”) referenced from `simulation_steps.md`.

It is written to be directly insertable as a replacement for the former `q_site_ros_step_a1_options.md`.

---

## Executive summary

**What this step decides:** which mechanistic source(s) of “Q-site ROS” you will model and/or experimentally target.

**Core correction:** “Q-site ROS” is not one thing. At minimum it must be disambiguated into:
- **Complex I site I_Q** (the ubiquinone/ubisemiquinone-linked ROS site in complex I; can operate during reverse and forward electron transport under some conditions).
- **Complex III site III_Qo** (the outer quinone-binding site of complex III associated with Qo-site semiquinone chemistry, especially prominent in antimycin-perturbed states).

**Measurement caution:** many common ROS assays are useful but artifact-prone; Step A1 must specify controls and detection chemistry, not just a probe name.

---

## Definitions and scientific background

### What is “ROS” here?

In this document, “ROS” refers primarily to:
- **Superoxide** (O2•−): the proximal radical formed by one-electron reduction of O2.
- **Hydrogen peroxide** (H2O2): largely produced by (spontaneous or SOD-catalysed) dismutation of superoxide and/or by downstream enzymatic redox reactions.

Key controlling variables for mitochondrial superoxide formation include **(i) availability/redox state of one-electron donors, (ii) local [O2], and (iii) the relevant rate constants**, and in mitochondria superoxide generation depends strongly on **Δp, NADH/NAD+ ratio, QH2/Q ratio, and O2 concentration**.
- Source: Murphy (2009) Biochem J. DOI:10.1042/BJ20081386. PMID:19061483.

### What does “Q-site” mean in this context?

“Q-site” refers to quinone-binding/reaction sites where ubiquinone/ubiquinol (coenzyme Q) participates in electron transfer and semiquinone intermediates can arise:

- **Complex I (NADH:ubiquinone oxidoreductase)** has a long quinone-binding channel in which ubiquinone binds and is reduced. Structural work directly visualises Q10 occupying the Q-binding site in mammalian complex I.
  - Source: Chung et al. (2022) Nat Commun. DOI:10.1038/s41467-022-30506-1.

- **Complex III (cytochrome bc1 complex)** has two quinone sites central to the Q-cycle:
  - **Qo site** (outer quinone/quinol site; often called “center o” in older literature; faces the intermembrane space side).
  - **Qi site** (inner quinone site; matrix side).

**Mechanistic link to ROS:** Semiquinone intermediates at Q sites can transfer an electron to O2 to generate superoxide, particularly when electron transfer through the normal pathway is perturbed.
- Evidence for Qo-site semiquinone involvement in superoxide: Cape et al. (2007) PNAS. DOI:10.1073/pnas.0702621104. PMID:17470780.
- Evidence for complex I having a Q-site ROS component distinct from the flavin site: Treberg et al. (2011) J Biol Chem. DOI:10.1074/jbc.M111.252502. PMID:21659507; Lambert & Brand (2004) J Biol Chem. DOI:10.1074/jbc.M406576200. PMID:15262965.

---

## Step A1 options

### Option set A: Model “Q-site ROS” as Complex I site I_Q (preferred when RET/FET at complex I is central)

#### A1-A: Site I_Q during reverse electron transport (RET-driven I_Q ROS)

**When this option is appropriate**
Choose this if your scenario involves conditions favouring reverse electron transport into complex I, classically:
- A highly reduced Q pool (high QH2/Q).
- A high protonmotive force (Δp), e.g., non-phosphorylating/high Δp states or ATP hydrolysis–driven Δp in isolated mitochondria-like preparations.
- Electron entry via succinate/complex II (or other Q-reducing pathways) can predispose to RET.

**Scientific basis**
- Experimental systems demonstrate rotenone-sensitive superoxide/H2O2 production from complex I during RET and show that the rate depends on Q-pool redox state and Δp.
  - Source: Treberg et al. (2011) J Biol Chem. DOI:10.1074/jbc.M111.252502. PMID:21659507.
- Reviews outline that significant superoxide production by isolated mitochondria occurs when Δp is high and the CoQ pool is reduced.
  - Source: Murphy (2009) Biochem J. DOI:10.1042/BJ20081386. PMID:19061483.

**What must be specified (or marked unspecified) for simulation**
- Q pool redox proxy (e.g., QH2/Q fraction or an equivalent redox potential).
- Δp and (if separated) ΔΨ and ΔpH.
- NADH/NAD+ redox state.
- Local [O2] (matrix O2).
- Electron supply route (succinate, glycerol-3-phosphate, fatty acid oxidation, etc.).

**Recommended diagnostic perturbations**
- **Rotenone or piericidin A**: Q-site inhibitors of complex I that can block electron flow and affect ROS depending on direction and context.
  - Source: Lambert & Brand (2004) J Biol Chem. DOI:10.1074/jbc.M406576200. PMID:15262965.
- **S1QELs**: small molecules that suppress superoxide/H2O2 production from site I_Q without (at appropriate doses) inhibiting oxidative phosphorylation—used to isolate site I_Q contributions with fewer metabolic confounders than classical inhibitors.
  - Source: Brand et al. (2016) Cell Metab. DOI:10.1016/j.cmet.2016.08.012. PMID:27667666.
  - Source: Wong et al. (2019) Free Radic Biol Med. DOI:10.1016/j.freeradbiomed.2019.09.006. PMID:31518685.

#### A1-B: Site I_Q during forward electron transport (FET-compatible I_Q ROS)

**When this option is appropriate**
Choose this if your biological/simulation conditions are predominantly forward electron transport through complex I (NADH → Q) but you still need site I_Q ROS.

**Scientific basis**
- Site I_Q can generate S1QEL-sensitive superoxide/H2O2 during forward electron transport as well as reverse electron transport, and the same Q-site inhibitors (rotenone/piericidin) suppress it.
  - Source: Gibbs et al. (2023) Biochem J. DOI:10.1042/BCJ20220611. PMID:36862427.

**Key implication for Step A1**
Do not assume “I_Q ROS implies RET” unless you have explicitly specified the thermodynamic direction of electron flow through complex I.

---

### Option set B: Model “Q-site ROS” as Complex III site III_Qo (preferred when Q-cycle perturbation or Qo chemistry is central)

#### A1-C: Basal/physiological III_Qo ROS (context dependent; usually modest)

**When this option is appropriate**
Choose this if your scenario emphasises complex III as a signalling ROS source (e.g., hypoxia signalling contexts) but you do not rely on strong pharmacological inhibition.

**Scientific basis**
- Complex III Qo-site superoxide has been implicated in signalling; however, assigning it is challenging because conventional inhibitors and knockdowns also change metabolism.
  - Source: Orr et al. (2015) Nat Chem Biol. DOI:10.1038/nchembio.1910. PMID:26368590.

#### A1-D: Antimycin-perturbed Q-cycle (high III_Qo ROS; mechanistically defensible)

**When this option is appropriate**
Choose this if Step A1 is explicitly about modelling the canonical “high ROS from complex III” condition produced by blocking Qi with antimycin A and driving semiquinone accumulation at Qo.

**Scientific basis**
- In antimycin-inhibited mitochondria, superoxide production shows a bell-shaped dependence on cytochrome b reduction, consistent with production from a semiquinone in the Qo site.
  - Source: Quinlan et al. (2011) J Biol Chem. DOI:10.1074/jbc.M111.267898. PMID:21708945.
- Direct detection and mechanistic discussion of Qo-site semiquinone radicals support the involvement of Qo semiquinone intermediates in superoxide generation.
  - Source: Cape et al. (2007) PNAS. DOI:10.1073/pnas.0702621104. PMID:17470780.
  - Source: Zhang et al. (2007) Biochim Biophys Acta. DOI:10.1016/j.bbabio.2007.04.004. PMID:17560537.

**Topology note**
- Complex III-derived superoxide has been reported as released to both sides of the inner membrane under some experimental conditions, whereas other work/reviews emphasise predominant release to the intermembrane space and highlight limitations of indirect topology assays.
  - Evidence for “both sides” under some conditions: Muller et al. (2004) J Biol Chem. DOI:10.1074/jbc.M407715200. PMID:15317809.
  - Earlier topology work: St-Pierre et al. (2002) J Biol Chem. DOI:10.1074/jbc.M207217200. PMID:12237311.
  - Discussion of debate and methodological caveats: Bleier & Dröse (2013) Biochim Biophys Acta. DOI:10.1016/j.bbabio.2012.12.002. PMID:23269318.

**Recommended diagnostic perturbations**
- **Stigmatellin or myxothiazol** (Qo inhibitors) can suppress some Qo-associated ROS signatures, but inhibitor interpretation is context‑dependent and can be confounded by upstream reduction.
  - Source (context dependence example): Starkov & Fiskum (2001) Biochem Biophys Res Commun. DOI:10.1006/bbrc.2001.4409. PMID:11237706.
- **S3QELs**: selective suppressors of site III_Qo electron leak without inhibiting oxidative phosphorylation (at appropriate doses), designed to deconfound “ROS vs metabolism”.
  - Source: Orr et al. (2015) Nat Chem Biol. DOI:10.1038/nchembio.1910. PMID:26368590.

---

### Option set C: Explicit mixed-source model (recommended when you cannot assume a single-site origin)

Choose this if your planned simulation/experiment cannot justify attributing “Q-site ROS” to a single locus.
- Use **S1QEL/S3QEL partitioning logic** (or the modelling analogue) to estimate contributions from I_Q vs III_Qo under the same base metabolic state.
  - Example of “site-specific suppressor” strategy: Treberg et al. (2011) J Biol Chem. DOI:10.1074/jbc.M111.252502. PMID:21659507; Orr et al. (2015) Nat Chem Biol. DOI:10.1038/nchembio.1910. PMID:26368590.

---

## Key supporting quotations (short, source-anchored)

> “These results support a two-site model of complex I superoxide production…”  
Treberg et al., 2011, *J Biol Chem*. DOI:10.1074/jbc.M111.252502. PMID:21659507.

> “…superoxide production peaks at intermediate Q-reduction state because it comes from a semiquinone in the outer quinone-binding site…”  
Quinlan et al., 2011, *J Biol Chem*. DOI:10.1074/jbc.M111.267898. PMID:21708945.

> “We report the first direct detection of a semiquinone radical generated by the Q(o) site…”  
Cape et al., 2007, *PNAS*. DOI:10.1073/pnas.0702621104. PMID:17470780.

> “…compounds that selectively eliminate superoxide production by complex III without altering oxidative phosphorylation…”  
Orr et al., 2015, *Nat Chem Biol*. DOI:10.1038/nchembio.1910. PMID:26368590.

> “Amplex Red is readily converted to resorufin … without requiring H2O2…”  
Miwa et al., 2016, *Free Radic Biol Med*. DOI:10.1016/j.freeradbiomed.2015.11.011. PMID:26577176.

> “Concentrations of 5–10 μM MitoSOX caused severe loss of ATP synthesis-linked respiration.”  
Roelofs et al., 2015, *Free Radic Biol Med*. DOI:10.1016/j.freeradbiomed.2015.05.032. PMID:26057935.

---

## Required Step A1 outputs (what you must select and write down)

### Mechanistic choice

You must explicitly choose ONE of:
1. **I_Q-only model** (Complex I Q-site ROS dominates).
2. **III_Qo-only model** (Complex III Qo-site ROS dominates).
3. **Mixed I_Q + III_Qo model** (partitioned).

…and record why, in one sentence, tied to at least one primary source above.

### Direction-of-flow assumption (only if Complex I is included)

Record one of:
- “Complex I electron flow is predominantly forward (FET: NADH → Q).”
- “Complex I electron flow is predominantly reverse (RET: QH2 → NAD+).”
- “Direction varies; simulation explicitly computes direction (preferred if feasible).”

Cite: Murphy (2009) PMID:19061483; Gibbs et al. (2023) PMID:36862427.

---

## Unspecified parameters checklist (must be explicitly marked “specified” or “unspecified”)

If any item is currently not stated elsewhere, mark it “unspecified” and choose a defensible default and cite it.

| Parameter | Status (specified/unspecified) | Suggested defaults (if unspecified) | Justification (primary source) |
|---|---|---|---|
| System context | — | Isolated mitochondria vs permeabilised cells vs intact cells vs pure complexes | Context strongly changes dominant ROS site. Murphy 2009 PMID:19061483 |
| Temperature | — | 30–37 °C (match your biological system); avoid mixing temps across steps | Temperature alters kinetics; match experimental design (general bioenergetics practice; cite Murphy 2009 PMID:19061483 for caution about extrapolation) |
| Oxygen tension / [O2] | — | Report %O2 and whether hypoxic/physiological; avoid “air-saturated” assumptions without stating it | ROS flux depends on local O2. Murphy 2009 PMID:19061483 |
| Δp / ΔΨ / ΔpH | — | State whether phosphorylating, non‑phosphorylating, uncoupled, or ATP hydrolysis–energised | Δp is a key driver of RET and ROS. Murphy 2009 PMID:19061483; Treberg 2011 PMID:21659507 |
| Q pool redox state | — | Specify proxy (QH2/Q fraction; cytochrome b reduction; etc.) | I_Q and III_Qo ROS depend on Q redox. Treberg 2011 PMID:21659507; Quinlan 2011 PMID:21708945 |
| Substrates feeding electrons | — | Explicit list and concentrations (e.g., succinate vs glutamate/malate vs fatty acids) | Substrate changes dominant ROS site. St-Pierre 2002 PMID:12237311; Murphy 2009 PMID:19061483 |
| Inhibitors/suppressors | — | Specify compound, concentration, addition order, exposure time | Inhibitor interpretation is context-dependent. Lambert 2004 PMID:15262965; Starkov 2001 PMID:11237706; Orr 2015 PMID:26368590 |
| ROS readout chemistry | — | H2O2 release via Amplex (with controls) and/or product-resolved HE/MitoSOX analytics | Assay artifacts are known. Miwa 2016 PMID:26577176; Zielonka 2010 PMID:20116425 |

---

## Measurement recommendations and artefact controls (if Step A1 includes experimental verification)

### If using Amplex Red / Amplex UltraRed (H2O2 release)

Minimum controls to state explicitly:
- Calibration curve with known H2O2 in the presence of your sample matrix.
- Controls for non-H2O2 conversion of probe (tissue-dependent).
- Record whether carboxylesterase artefacts are plausible in your preparation and what you did about it.

Evidence: Miwa et al. show carboxylesterase can convert Amplex Red to resorufin without H2O2 (PMID:26577176). Murphy (2009) discusses Amplex-based detection of H2O2 and matrix sinks (PMID:19061483).

### If using hydroethidine / MitoSOX

You must state whether you are:
1) using **red fluorescence imaging as a qualitative proxy**, or  
2) doing **product-resolved quantitation** (preferred; e.g., HPLC/MS distinguishing 2‑hydroxyethidium from ethidium and other products).

Key evidence:
- Zielonka & Kalyanaraman review the non-specificity of HE/MitoSOX red fluorescence (PMID:20116425).
- Roelofs et al. show MitoSOX can perturb mitochondrial respiration at micromolar concentrations (PMID:26057935).

---

## Flowchart for Step A1 decision logic

(Scientific rationale: I_Q and III_Qo are the two dominant “quinone-site” ROS options; their controlling variables differ materially. Key sources: Murphy 2009 PMID:19061483; Treberg 2011 PMID:21659507; Quinlan 2011 PMID:21708945; Orr 2015 PMID:26368590.)

```mermaid
flowchart TD
  A[Start Step A1: Define "Q-site ROS"] --> B{Which complex is central to your hypothesis?}
  B -->|Complex I| C[Select site I_Q model]
  B -->|Complex III| D[Select site III_Qo model]
  B -->|Both / uncertain| E[Select mixed I_Q + III_Qo model]

  C --> F{Electron flow direction through complex I?}
  F -->|RET likely| G[I_Q-RET scenario]
  F -->|FET likely| H[I_Q-FET scenario]
  F -->|Variable| I[Compute direction explicitly / branch by state]

  D --> J{Perturbation state?}
  J -->|Antimycin-inhibited Q-cycle| K[High III_Qo ROS model]
  J -->|No strong inhibition| L[Basal/physiological III_Qo ROS model]

  E --> M[Plan site-partition strategy]
  M --> N[Use S1QEL/S3QEL-informed constraints or analogues]

  G --> O[Write required parameters + controls]
  H --> O
  I --> O
  K --> O
  L --> O
  N --> O

  O --> P[Output Step A1 choices + parameter list]
```

---

## Original vs revised statements (two-column comparison)

**Important:** the literal original statements from your prior file were not available in this tool session. The table below therefore captures *common* “Q-site ROS” statement patterns that are frequently incorrect/underspecified and the scientifically defensible replacement formulations used above.

| Common problematic (original-style) statement | Revised, scientifically supported statement |
|---|---|
| “Q-site ROS is produced at complex I.” | “Q-site ROS must specify *which* site (I_Q vs I_F vs III_Qo); complex I includes at least I_F and I_Q components.” (Treberg 2011 PMID:21659507; Lambert 2004 PMID:15262965) |
| “I_Q ROS implies reverse electron transport.” | “Site I_Q can generate S1QEL-sensitive ROS during both RET and FET; do not equate I_Q with RET without an electron-flow direction test.” (Gibbs 2023 PMID:36862427) |
| “Antimycin increases ROS from complex III.” | “Under antimycin A, Qo-site superoxide can increase and shows a bell-shaped dependence on cytochrome b reduction, consistent with a Qo-site semiquinone source.” (Quinlan 2011 PMID:21708945) |
| “Qo-site semiquinone is hypothetical.” | “Qo-site semiquinone radicals have been directly detected by EPR under conditions that trap/accumulate them.” (Cape 2007 PMID:17470780; Zhang 2007 PMID:17560537) |
| “Complex III releases superoxide only to the intermembrane space.” | “Complex III superoxide topology is condition- and method-dependent; intact-mitochondria work supports release to the intermembrane space and (under some conditions) evidence also supports matrix-facing effects; quantitative partitioning remains debated.” (St-Pierre 2002 PMID:12237311; Muller 2004 PMID:15317809; Bleier 2013 PMID:23269318) |
| “MitoSOX fluorescence directly measures mitochondrial superoxide.” | “MitoSOX/HE red fluorescence is not a reliable standalone indicator of intracellular superoxide without product-resolved analytics and probe-bioenergetics controls, because of side products and potential mitochondrial perturbation.” (Zielonka 2010 PMID:20116425; Roelofs 2015 PMID:26057935) |
| “Amplex Red reports H2O2 release.” | “Amplex assays can report H2O2 release but require explicit artefact controls; in some samples Amplex Red can be converted to resorufin without H2O2.” (Miwa 2016 PMID:26577176; Murphy 2009 PMID:19061483) |
| “Inhibitors cleanly isolate the ROS site.” | “Conventional inhibitors often confound ROS attribution by changing upstream redox states and metabolism; site-specific suppressors (S1QEL/S3QEL) were developed to deconfound ROS vs metabolism.” (Brand 2016 PMID:27667666; Orr 2015 PMID:26368590) |

---

## References (full citations; primary sources prioritised)

Bleier, L., & Dröse, S. (2013). *Superoxide generation by complex III: from mechanistic rationales to functional consequences.* Biochimica et Biophysica Acta (Bioenergetics), 1827(11–12), 1320–1331. DOI:10.1016/j.bbabio.2012.12.002. PMID:23269318.

Brand, M. D., Goncalves, R. L. S., Orr, A. L., Vargas, L., Gerencser, A. A., Borch Jensen, M., et al. (2016). *Suppressors of Superoxide-H2O2 Production at Site I_Q of Mitochondrial Complex I Protect against Stem Cell Hyperplasia and Ischemia-Reperfusion Injury.* Cell Metabolism, 24(4), 582–592. DOI:10.1016/j.cmet.2016.08.012. PMID:27667666.

Cape, J. L., Bowman, M. K., & Kramer, D. M. (2007). *A semiquinone intermediate generated at the Qo site of the cytochrome bc1 complex: importance for the Q-cycle and superoxide production.* Proceedings of the National Academy of Sciences USA, 104(19), 7887–7892. DOI:10.1073/pnas.0702621104. PMID:17470780.

Chung, I., Wright, J. J., Bridges, H. R., Ivanov, B. S., Biner, O., Pereira, C. S., et al. (2022). *Cryo-EM structures define ubiquinone-10 binding to mitochondrial complex I and conformational transitions accompanying Q-site occupancy.* Nature Communications, 13, 2758. DOI:10.1038/s41467-022-30506-1.

Gibbs, E. T., Lerner, C. A., Watson, M. A., Wong, H.-S., Gerencser, A. A., & Brand, M. D. (2023). *Site IQ in mitochondrial complex I generates S1QEL-sensitive superoxide/hydrogen peroxide in both the reverse and forward reactions.* Biochemical Journal, 480(5), 363–384. DOI:10.1042/BCJ20220611. PMID:36862427.

Lambert, A. J., & Brand, M. D. (2004). *Inhibitors of the quinone-binding site allow rapid superoxide production from mitochondrial NADH:ubiquinone oxidoreductase (complex I).* Journal of Biological Chemistry, 279(38), 39414–39420. DOI:10.1074/jbc.M406576200. PMID:15262965.

Miwa, S., Treumann, A., Bell, A., Vistoli, G., Nelson, G., Hay, S., & von Zglinicki, T. (2016). *Carboxylesterase converts Amplex red to resorufin: Implications for mitochondrial H2O2 release assays.* Free Radical Biology and Medicine, 90, 173–183. DOI:10.1016/j.freeradbiomed.2015.11.011. PMID:26577176.

Muller, F. L., Liu, Y., & Van Remmen, H. (2004). *Complex III releases superoxide to both sides of the inner mitochondrial membrane.* Journal of Biological Chemistry, 279(47), 49064–49073. DOI:10.1074/jbc.M407715200. PMID:15317809.

Murphy, M. P. (2009). *How mitochondria produce reactive oxygen species.* Biochemical Journal, 417(1), 1–13. DOI:10.1042/BJ20081386. PMID:19061483.

Orr, A. L., Vargas, L., Turk, C. N., Baaten, J. E., Matzen, J. T., Dardov, V. J., et al. (2015). *Suppressors of superoxide production from mitochondrial complex III.* Nature Chemical Biology, 11(11), 834–836. DOI:10.1038/nchembio.1910. PMID:26368590.

Quinlan, C. L., Gerencser, A. A., Treberg, J. R., & Brand, M. D. (2011). *The mechanism of superoxide production by the antimycin-inhibited mitochondrial Q-cycle.* Journal of Biological Chemistry, 286(36), 31361–31372. DOI:10.1074/jbc.M111.267898. PMID:21708945.

Roelofs, B. A., Ge, S. X., Studlack, P. E., & Polster, B. M. (2015). *Low micromolar concentrations of the superoxide probe MitoSOX uncouple neural mitochondria and inhibit complex IV.* Free Radical Biology and Medicine, 86, 250–258. DOI:10.1016/j.freeradbiomed.2015.05.032. PMID:26057935.

Starkov, A. A., & Fiskum, G. (2001). *Myxothiazol induces H(2)O(2) production from mitochondrial respiratory chain.* Biochemical and Biophysical Research Communications, 281(3), 645–650. DOI:10.1006/bbrc.2001.4409. PMID:11237706.

St-Pierre, J., Buckingham, J. A., Roebuck, S. J., & Brand, M. D. (2002). *Topology of superoxide production from different sites in the mitochondrial electron transport chain.* Journal of Biological Chemistry, 277(47), 44784–44790. DOI:10.1074/jbc.M207217200. PMID:12237311.

Treberg, J. R., Quinlan, C. L., & Brand, M. D. (2011). *Evidence for two sites of superoxide production by mitochondrial NADH-ubiquinone oxidoreductase (complex I).* Journal of Biological Chemistry, 286(31), 27103–27110. DOI:10.1074/jbc.M111.252502. PMID:21659507.

Wong, H.-S., Monternier, P.-A., & Brand, M. D. (2019). *S1QELs suppress mitochondrial superoxide/hydrogen peroxide production from site I_Q without inhibiting reverse electron flow through Complex I.* Free Radical Biology and Medicine, 143, 545–559. DOI:10.1016/j.freeradbiomed.2019.09.006. PMID:31518685.

Zhang, H., Osyczka, A., Dutton, P. L., & Moser, C. C. (2007). *Exposing the complex III Qo semiquinone radical.* Biochimica et Biophysica Acta (Bioenergetics), 1767(7), 883–887. DOI:10.1016/j.bbabio.2007.04.004. PMID:17560537.

Zielonka, J., & Kalyanaraman, B. (2010). *Hydroethidine- and MitoSOX-derived red fluorescence is not a reliable indicator of intracellular superoxide formation: another inconvenient truth.* Free Radical Biology and Medicine, 48(8), 983–1001. DOI:10.1016/j.freeradbiomed.2010.01.028. PMID:20116425.

---

## End of Step A1 options
```

## Limitations

The replacement text above is fully literature‑anchored and intended to be paste‑ready, but it is **not a literal rewrite of your current file** because the tool session could not access the repository files you referenced. Consequently, the “original vs revised” table is framed as a correction map for *typical* problematic statements rather than a direct two‑column diff of your file’s exact sentences. citeturn22view0turn18search1turn9view0turn10view0