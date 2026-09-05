---
type: summary
status: active
tags: [collagen-hydrolysis, activation-energy, thermal-age, radiocarbon, audit, palaeoproteomics, diagenesis]
relatedTo: [ea-collagen-thermal-age-space-is-circular, palaeoprot-publications, thermal-age-specification]
---

# Handover — Manuscript Revision and Statistical Audit
**Date & Time:** 2026-09-05 12:05 (+02:00)

Supersedes the analytical claims in `PROJECT_HANDOVER_COLLAGEN_HYDROLYSIS_EA.md` (2026-09-05 10:25).
That document remains valid for data locations, file paths and literature inventory.

---

## 1. Read this first

The manuscript was revised and then audited. **The audit did not confirm the
headline result.** Anyone picking this up should understand four things before
touching the text:

1. **$E_a = 133.4 \pm 3.6\ \mathrm{kJ\,mol^{-1}}$ is not supportable.** The stated
   precision is roughly an order of magnitude too tight, and the site-clustered
   confidence interval $[85, 192]$ contains the laboratory value of $173$.
2. **The estimate is not stable under reclassification.** Admitting undifferentiated
   bone as probable collagen moves it to $116.8$ (bin maxima) or $94.6$ ($Q_{95}$).
   Every correction so far has moved it downward: $133.4 \to 130.8 \to 116.8/94.6$.
3. **The matched-control design described in earlier drafts does not exist.** The
   cohort is size-matched, not site-matched. Only 388 of 4,695 collagen sites carry
   any control.
4. **The strongest quantitative support is no longer the field inversion.** It is
   Ortner et al. (1972), $E_a = 132.1\ \mathrm{kJ\,mol^{-1}}$ on *intact* bone,
   which is independent of the radiocarbon data entirely.

What survives robustly: bone collagen shows a sharp, **collagen-specific** ceiling
that plant controls do not share, and its form is inconsistent with $173$ anchored
to temperate radiocarbon limits. The direction is secure; the number is not.

---

## 2. Current state of files

### Manuscript
| File | State |
| :--- | :--- |
| `manuscript_thermal_age_collagen_hydrolysis.md` | **Current.** 608 lines, revised + audited |
| `manuscript_thermal_age_collagen_hydrolysis.pdf` | Current (9.2 MB, rebuilt 2026-09-05 12:00) |
| `manuscript_thermal_age_collagen_hydrolysis.html` | Current |
| `manuscript_thermal_age_collagen_hydrolysis.md.bak` | Pre-revision original (469 lines) |
| Mirror | `GitHub/manuscripts/07_Thermal_Age_and_Collagen_Hydrolysis_Ea/` |

### Analysis added this session
| File | Purpose |
| :--- | :--- |
| `sensitivity_charcoal_only.py` / `_report.md` | Charcoal-only restriction + full robustness audit |
| `probable_collagen_reclassification.py` / `_report.md` | Undifferentiated bone admitted as probable collagen |
| `censored_regression_ea.py` | Matched-control logistic estimator — **documents why it is not identified** |
| `plot_comprehensive_deconvolution_4panel.py` | Figure 2 regenerated; `.bak` holds the original |
| `build_manuscript_html_pdf.py` | PDF step fixed (see §6) |

---

## 3. What was changed in the manuscript, and why

### 3.1 Arithmetic that did not reproduce
Table 3's thermal ages could not be recovered from the stated $E_a$ and temperatures.
Recomputed from Equation (1) and now reproducible from the adjacent columns:

| Specimen | published | recomputed |
| :--- | ---: | ---: |
| Ellesmere bear | 31.8 ka | **46.8 ka** |
| Arctic camel | 27.7 ka | **40.8 ka** |
| Yukon horse | 8.2 ka | **12.5 ka** |
| Sima / Boxgrove | 298 ka | **234 ka** |
| Dmanisi dentin | 2.05 Ma | **2.16 Ma** |
| Broken River karst | 1.82 Ma | **5.14 Ma** |
| Swartkrans | 4.8 Ma | **6.78 Ma** |

### 3.2 The "unification" claim was false
Dmanisi dentin at $1.77$ Ma / $+11\,^\circ$C is $2.16$ Ma @ 10 °C — **38× beyond
$Q_{99}$**. Broken River is 90× beyond. No admissible $E_a$ absorbs these: lowering
it destroys the tropical collapse, raising it reinstates the Arctic paradox.
Section 5 was rewritten so benchmarks are scored on one reproducible axis, with
explicit statements of where the model holds and where it does not.

Crucially, **all discrepancies are unexpected survival, never unexpected loss** —
the signature of a correct rate law applied to the wrong compartment.

### 3.3 The occluded-niche mechanism (M. Collins)
Salamon, Tuross, Arensburg & Weiner (2005, *PNAS* 102:13783) showed fossil bone
contains NaOCl-resistant intergrown crystal aggregates holding better-preserved
DNA — a privileged intra-crystalline niche. §5.3 now uses this as the unifying
explanation for Dmanisi, Broken River and the enamel record: one occlusion axis,
not three odd sites. It also explains why the effect is rare *in practice* —
standard acid-demineralisation-and-gelatinisation destroys that fraction before
assay. §5.3.2 states the test: NaOCl-resistant fraction vs whole-bone extract.

### 3.4 Circularity of thermal-age-space fitting (§4.4, new)
`empirical_ea_fitting_report.md` "validated" 173 by sweeping $E_a$, recomputing
thermal age under each candidate, and picking the sharpest boundary. **That cannot
identify $E_a$** — thermal age is defined by an integral containing the assumed
value, so the regression abscissa is a function of the parameter being estimated.
The sweep has no interior optimum (23.6 ka at $E_a{=}100$ → 29.4 ka at 200, monotonic).
This is likely how 173 survived twenty years of apparent testing.

### 3.5 Other corrections
- Ortner (1972, intact bone, 132.1) vs Von Endt & Ortner (1984, powdered, 183.3) —
  same temperatures, 51 kJ/mol apart. **The discriminant is maceration, not temperature.**
- Ref 44 was wrong in journal, title, volume and DOI → Peters et al. 2023,
  *Commun. Earth Environ.* **4**, 438, doi 10.1038/s43247-023-01114-8.
- Added Salamon [45] and Ni [46]. $R^2 \to$ Koenker–Machado $R^1$.
- $A = 1.98\times10^{16} \to 3.11\times10^{16}$ ($e^{37.976}$).
- New §5.4: bulk-yield (1%) vs LC-MS/MS ($10^{-6}$) ceilings differ by exactly 3×.
- New §6.3 (limitations — there were none) and §6.4 (six falsifiable predictions).
- §3.1.1: control materials are **not** interchangeable chronometers.

---

## 4. The audit — full numbers

### 4.1 The binned fit rests on eight bins
No admissible bins at 21–24 °C. Effective sample size is 8, not 18,101.

| Test | Result |
| :--- | :--- |
| Drop the 25 °C bin (40 dates, **5 sites**) | 130.8 → **86.0** |
| Bin max / $Q_{95}$ / median | 130.8 / 124.5 / **77.5** |
| Equal-*n* resampling | **119.9** [96.4, 140.4] |
| $T_c$ from 12–16 °C | **113.7–133.4** |
| Site-clustered bootstrap | **[85.1, 191.6]** — contains 173 |
| Charcoal-witnessed sites only | **88.6** (SE 52; 18.7% retained) |

Sample size per bin is anti-correlated with temperature ($r = -0.76$) while bin
maximum rises mechanically with bin size ($r = +0.56$): **part of the collapse is
sampling, not chemistry.**

### 4.2 The matched-control estimator is not identified
Proposed model: with controls as an availability denominator,
$\mathrm{logit}\,P(\text{collagen} \mid a,T) = c + b(1000/T) - g\ln a$, giving
$E_a = 1000R\,b/g$.

It fails. The collagen **share rises** with age:

| median age (yr) | 840 | 2,170 | 4,260 | 11,870 | 31,540 |
| :--- | ---: | ---: | ---: | ---: | ---: |
| collagen share | 23% | 30% | 53% | 76% | **90%** |

The age coefficient takes the wrong sign; the model is rejected pooled and with
site fixed effects. Cause: Holocene settlement archaeology is dated on charcoal
and seeds, Palaeolithic and megafaunal contexts on bone. **Material choice tracks
research tradition, not preservation, and that confound exceeds the signal.**
A Tobit on $\ln$(age) right-censored at 42 ka returns $E_a = 8.2$ — confirming the
age distribution is archaeology, not chemistry.

### 4.3 Reclassification (probable collagen)
Excluding cremated/calcined/burnt (1,067) leaves 19,444 probable-collagen records.

| | baseline | + probable collagen |
| :--- | ---: | ---: |
| collagen records | 18,084 | **37,528** |
| admissible warm bins | 8 | **12** |
| $E_a$ (bin max) | 130.8 | **116.8** |
| $E_a$ ($Q_{95}$) | 124.5 | **94.6** |

**The 21–24 °C gap fills completely**; single-bin leverage disappears. Warm-bin
maxima rise sharply — 16 °C from 13,070 → **33,550** yr; new 21 °C bin at **33,490** yr.
The sharp collapse above 14 °C in Figure 2A is therefore **substantially an artefact
of which bone dates were counted as collagen**. It remains real and collagen-specific,
but it is gentler, and the implied barrier lower.

### 4.4 Material inventory (185,189 master records)
| Class | N | Note |
| :--- | ---: | :--- |
| `NON_COLLAGEN_ORGANIC_CONTROL` | 100,809 | charcoal 13.4k, wood 1.5k, charred grain, seeds |
| `OTHER` | 29,786 | **12,705 blank** — parsing gap, not a material class |
| `BONE_UNDIFFERENTIATED` | 22,198 | **larger than the collagen class**; 21,131 non-cremated |
| `COLLAGEN` | 20,347 | explicit ultrafiltered / XAD gelatin |
| `SHELL_CONTROL` | 9,561 | |
| `SEDIMENT` | 2,296 | |
| `BONE_APATITE` | **192** | ~70 genuinely bone |

Apatite's rarity reflects that the method is disliked, difficult and controversial —
a fact about practice, not a harmonisation failure. It cannot supply the
specimen-level control the design would ideally use.

---

## 5. Errors introduced during revision (corrected, but be aware)

Two fabrications entered the text during this session and were later caught and fixed.
Flagged here because similar claims may have propagated to derived documents:

1. **`BONE_APATITE (N = 2,946)`** was written into §3.1 with an argument built on
   specimen-level apatite pairing. The true figure is 192 in the master and **zero**
   in the cohort. Removed and replaced with an explicit statement of absence.
2. **"Controls were matched one-to-one on site identity and stratigraphic horizon"**
   was invented to rationalise the equal $N = 18{,}101$. False: only 388 of 4,695
   collagen sites carry a control. Corrected to "size-matched, not site-matched".

A third, softer case: `±8 kJ/mol` for breakpoint uncertainty was asserted before
being computed. The real range is 113.7–133.4, and the honest interval is ±25–40.

---

## 6. Tooling notes

- **`build_manuscript_html_pdf.py`** used the legacy `--headless` Edge flag, which
  silently no-ops (exit 0, no file) when a browser instance is running. Now uses
  `--headless=new` with an isolated profile and **verifies the output mtime changed**,
  raising if not. Do not revert.
- **`plot_comprehensive_deconvolution_4panel.py`** now falls back to the repo data
  copy if `D:` is unavailable and mirrors output to three directories.
- Figure 2 Panel D benchmarks are split into filled (inside envelope) vs hollow
  (occluded matrix); Panel B adds the Ortner intact-bone line. Three stale callouts
  were removed ("Australian Paradox / polymer-in-a-box", "sit directly on the
  boundary", "Validates field $E_a$").

---

## 7. Recommended next steps, in priority order

1. **Rebuild the primary cohort on the reclassified data** (§4.3) and refit from
   there, rather than patching current numbers. This is the single highest-value
   action: it fixes the warm-bin gap and the leverage problem simultaneously.
   Requires re-running `run_thermal_integration.py` over the extended cohort.
2. **Re-state the claim as a bound**, not a point estimate. "The field barrier is
   well below the powdered-bone laboratory value, with independent laboratory
   support at 132.1 on intact bone" is defensible. "133.4 ± 3.6" is not.
3. **Regenerate Figure 2 Panel A** on the reclassified cohort — the collapse it
   currently shows is partly a classification artefact.
4. **Resolve the two `TODO(authors)` markers** in the manuscript: ref 18 carries
   both the camel and the Ellesmere giant bear (Rybczynski 2013 is the camel paper
   only; the bear needs its own citation), and the Sima $\bar{T} = +7\,^\circ$C
   assumption needs Krapp2021 (800 ka, covers the site) or clumped-isotope work.
5. **Run the NaOCl experiment** of §5.3.2 if any of the discrepant material is
   accessible. It is the decisive test of the central mechanistic claim.
6. **Consider whether the paper's framing should change** — the audit arguably
   makes a stronger methodological paper (why boundary fits on radiocarbon
   compilations fail, and what it takes to do them properly) than a kinetics paper.

---

## 8. Related memory

- [[ea-collagen-thermal-age-space-is-circular]] — why thermal-age-space fits cannot test $E_a$
- See also `sensitivity_charcoal_only_report.md`, `probable_collagen_reclassification_report.md`
