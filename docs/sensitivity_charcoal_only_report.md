# Charcoal-Only Sensitivity Analysis — and a Robustness Audit of $E_a$
**Date & Time:** 2026-09-05 (generated from `sensitivity_charcoal_only.py` + leverage/confound tests)

> [!WARNING]
> **The requested sensitivity analysis does not support the manuscript's headline claim.**
> $E_a = 133.4 \pm 3.6\ \mathrm{kJ\,mol^{-1}}$ is not recoverable as a robust
> estimate from this cohort. The result is dominated by a single temperature bin,
> is confounded with sample size, and its site-clustered confidence interval
> includes $173\ \mathrm{kJ\,mol^{-1}}$. Details below.

---

## 1. Cohort composition

| Quantity | N | Sites |
| :--- | ---: | ---: |
| Purified bone collagen | 18,101 | 4,696 |
| Non-collagen organic controls (all) | 18,101 | 8,816 |
| Controls, **strict charcoal only** | 15,154 | 7,436 |
| Collagen at charcoal-witnessed sites | **3,385** | 1,007 |

Retention on restriction to charcoal-witnessed sites: **18.7%**.

> [!IMPORTANT]
> **The analytical cohort contains no bone apatite determinations at all.**
> `BONE_APATITE` has n = 192 in the harmonised master (185,189 records) and none
> survive the coordinate / quality / pairing filters. Of those 192, many are not
> even bone (`pedogenic carbonate`, `geological carbonate`, `calcium carbonate or
> plaster`); genuine bone apatite is roughly 70 records.
> The apatite pairing described in manuscript §3.1 is therefore a **proposed**
> design that this dataset does not implement. §3.1 must be corrected.

---

## 2. The fit rests on eight bins

Bins entering the $E_a$ fit ($n \ge 10$, $13.8 \le \bar{T} \le 25\,^\circ$C):

| $\bar{T}$ bin | n | sites | max $^{14}$C age | $Q_{95}$ |
| ---: | ---: | ---: | ---: | ---: |
| 14 | 318 | 109 | 41,200 | 27,727 |
| 15 | 436 | 95 | 20,260 | 14,120 |
| 16 | 340 | 92 | 13,070 | 10,456 |
| 17 | 99 | 23 | 12,350 | 11,390 |
| 18 | 35 | 15 | 16,795 | 16,153 |
| 19 | 36 | 16 | 12,150 | 11,530 |
| 20 | 85 | 15 | 15,700 | 6,924 |
| 25 | **40** | **5** | 3,000 | 2,707 |

There are **no admissible bins between 21 and 24 °C** (n = 9, 1, 2, 0). The
effective sample size for the headline parameter is 8, not 18,101. Note also
that ages are flat and non-monotonic from 16–20 °C (12,350 / 16,795 / 12,150 /
15,700): the "exponential collapse" is a step between 14 and 16 °C followed by a
plateau, then a single distant anchor.

---

## 3. Robustness tests

### 3.1 Choice of order statistic

| Statistic used as the envelope | $E_a$ (kJ mol$^{-1}$) |
| :--- | ---: |
| bin maximum (as published) | **130.8** |
| bin $Q_{95}$ | 124.5 |
| bin median | 77.5 |

The parameter is not invariant to the summary chosen.

### 3.2 Leverage of the 25 °C bin

| Fit | $E_a$ | bins |
| :--- | ---: | ---: |
| all 8 bins, max | **130.8** | 8 |
| **drop the 25 °C bin**, max | **86.0** | 7 |
| drop the 25 °C bin, $Q_{95}$ | 97.1 | 7 |

Removing one bin — 40 dates from 5 sites — moves $E_a$ by **−44.8 kJ mol$^{-1}$**.
The headline value is carried by that bin's lever arm.

### 3.3 Sample-size confound

Bin sample size falls steeply across the fitted range and is strongly
anti-correlated with temperature:

* corr(log n, $\bar{T}$) = **−0.756**
* corr(log n, log max age) = **+0.561**

A bin maximum is an extreme order statistic and rises with n mechanically.
Because warm bins are smaller, part of the observed decline in maximum age is
sampling, not chemistry. Equal-n subsampling (35 dates per bin, 200 replicates):

| Statistic | $E_a$ | 95% interval |
| :--- | ---: | :--- |
| max, equal-n | **119.9** | [96.4, 140.4] |
| $Q_{95}$, equal-n | **107.7** | [79.2, 134.5] |

### 3.4 Site-clustered bootstrap (B = 400)

| Control set | $E_a$ | SE | 95% CI |
| :--- | ---: | ---: | :--- |
| All organic controls | 130.8 | 26.7 | **[85.1, 191.6]** |
| Strict charcoal only | 88.6 | 52.2 | [−36.7, 159.7] |

**The interval on the full cohort includes 173 kJ mol$^{-1}$.** Determinations
within a site share a thermal history; once that dependence is respected, the
data do not exclude the laboratory value.

### 3.5 Breakpoint sensitivity

| $T_c$ (°C) | $E_a$ all | $E_a$ charcoal | bins |
| ---: | ---: | ---: | ---: |
| 12.0 | 132.4 | 129.2 | 7 |
| 12.5 | 133.4 | 111.8 | 6 |
| 13.0 | 133.4 | 111.8 | 6 |
| 13.5 | 130.8 | 88.6 | 5 |
| 14.0 | 130.8 | 88.6 | 5 |
| 14.5 | 113.7 | 67.9 | 4 |
| 15.0 | 113.7 | 67.9 | 4 |
| 15.5 | 116.7 | — | 3 |
| 16.0 | 116.7 | — | 3 |

Range across plausible $T_c$: **113.7–133.4** (all controls), **67.9–129.2**
(charcoal). The manuscript's stated $\pm 8\ \mathrm{kJ\,mol^{-1}}$ understates this.

### 3.6 The requested charcoal-only result

$E_a$ falls from **130.8** to **88.6 kJ mol$^{-1}$** (−42.1) when restricted to
charcoal-witnessed sites. However, the restriction leaves 5 bins and a bootstrap
SE of 52, so this shift is **not itself well determined**; it should be read as
"the restriction destroys the precision of the estimate", not as "the true value
is 88.6".

---

## 4. What can honestly be said

**Supported:**
1. Bone collagen maximum ages collapse sharply above ~14 °C while plant controls
   do not. This contrast is large, is visible in the raw bins, and is not a
   recovery artefact.
2. The collapse is inconsistent with the *shape* implied by extrapolating
   $E_a = 173$ anchored to temperate radiocarbon limits.
3. The direction of the field estimate is below the laboratory value under every
   variant tested.

**Not supported:**
1. $E_a = 133.4 \pm 3.6$. The precision is roughly an order of magnitude
   overstated; a defensible interval is closer to **$\pm 25$–$40$**.
2. "Incompatible with 173." The site-clustered CI [85, 192] contains it.
3. Any deep-time extrapolation quoted to three significant figures. Propagating
   the real interval to $-10.5\,^\circ$C spans well over an order of magnitude in
   predicted survival.

---

## 5. Recommended next steps

1. **Fix §3.1** — remove the apatite cohort claim (it does not exist in the data).
2. **Recover the 21–24 °C gap.** The single 25 °C bin (5 sites) carries the fit.
   Targeted compilation of warm-climate collagen dates is the highest-value
   addition and would settle the slope.
3. **Model at the determination level, not the bin level** — a censored-regression
   or survival model (Tobit / Cox with the 42 ka ceiling as right-censoring) uses
   all 18,101 dates, avoids the order-statistic and equal-n problems, and yields
   a defensible interval. This is the correct estimator for this question.
4. **Re-state the claim** as a bound: the field barrier is *below* the powdered-bone
   laboratory value, with Ortner's intact-bone 132.1 kJ mol$^{-1}$ as independent
   support — rather than as a point estimate to three significant figures.
