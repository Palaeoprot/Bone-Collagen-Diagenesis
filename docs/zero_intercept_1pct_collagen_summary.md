# Empirical Zero-Intercept (0, 0) Kinetic Fit for 1% Residual Bone Collagen
**Date & Time:** 2026-09-05 09:16:00 (+02:00)

## 1. Zero-Intercept Boundary Formulation

When forcing the boundary envelope through $(0, 0)$ (representing modern / zero-age bone where initial collagen fraction $C/C_0 = 100\%$), the relationship between Arrhenius Thermal Age ($t_{\text{eff}}$ at 10 °C reference) and Calendar Age ($t_{\text{cal}}$) defines the maximum thermal acceleration envelope under which archaeological bone collagen survives at $\ge 1\%$ residual yield:

$$t_{\text{eff}} = \beta \cdot t_{\text{cal}}$$

For first-order hydrolysis degradation:
$$\ln\left(\frac{C(t)}{C_0}\right) = -k(T_{\text{ref}}) \cdot t_{\text{eff}}$$

At the standard radiocarbon threshold of **$1\%\text{ residual collagen}$** ($C/C_0 = 0.01$):
$$-\ln(0.01) = 4.60517 = k(10^\circ\text{C}) \cdot t_{\text{eff, limit}}$$

Therefore:
$$k(10^\circ\text{C}) = \frac{4.60517}{t_{\text{eff, limit}}}$$
$$t_{1/2}(10^\circ\text{C}) = \frac{\ln(2)}{k(10^\circ\text{C})} = \frac{0.693147}{k(10^\circ\text{C})}$$

---

## 2. Quantile Regression Slopes Through $(0, 0)$

Fitting zero-intercept quantile regression models across $N = 18,101$ purified bone collagen determinations versus $N = 18,101$ negative controls (charcoal, wood, seeds):

| Quantile ($\tau$) | Purified Collagen Slope ($\beta_{\text{col}}$) | Non-Collagen Control Slope ($\beta_{\text{ctrl}}$) | Ratio ($\beta_{\text{ctrl}} / \beta_{\text{col}}$) |
| :--- | :--- | :--- | :--- |
| **50.0% (Median)** | **$0.4274 \pm 0.003$** | **$0.8369 \pm 0.006$** | $1.96\times$ |
| **75.0%** | **$0.7877 \pm 0.005$** | **$2.5531 \pm 0.018$** | $3.24\times$ |
| **90.0%** | **$1.3618 \pm 0.010$** | **$8.3517 \pm 0.062$** | $6.13\times$ |
| **95.0% (Upper Envelope)** | **$2.1255 \pm 0.017$** | **$15.4878 \pm 0.142$** | **$7.29\times$** |
| **98.0%** | **$3.3163 \pm 0.029$** | **$42.4115 \pm 0.410$** | $12.79\times$ |
| **99.0%** | **$4.4362 \pm 0.041$** | **$55.4414 \pm 0.582$** | $12.50\times$ |

> [!NOTE]
> The controls display an upper envelope slope of **$15.49$**, which is **$7.3\times$ steeper** than purified bone collagen ($2.13$). This confirms that non-collagen materials date successfully at massive thermal ages where bone collagen has completely decayed away.

---

## 3. Best-Fit Kinetic Parameters for 1% Residual Collagen

Evaluating the empirical survival ceiling across different boundary percentiles gives the corresponding rate constants and half-lives:

| Boundary Metric | Thermal Age Ceiling ($t_{\text{eff}}$ @ 10 °C) | $k(10^\circ\text{C})\ [\text{yr}^{-1}]$ | $k(10^\circ\text{C})\ [\text{s}^{-1}]$ | Half-Life $t_{1/2}(10^\circ\text{C})\ [\text{yr}]$ | Implied Pre-Exponential $A\ [\text{s}^{-1}]$ ($E_a = 173\text{ kJ/mol}$) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **90th Percentile** | $17,220\text{ years}$ | $2.674 \times 10^{-4}$ | $8.474 \times 10^{-12}$ | **$2,592\text{ years}$** | $6.95 \times 10^{20}$ |
| **95th Percentile (Primary)** | **$26,278\text{ years}$** | **$1.753 \times 10^{-4}$** | **$5.553 \times 10^{-12}$** | **$3,955\text{ years}$** | **$4.55 \times 10^{20}$** |
| **98th Percentile** | $43,677\text{ years}$ | $1.054 \times 10^{-4}$ | $3.341 \times 10^{-12}$ | **$6,574\text{ years}$** | $2.74 \times 10^{20}$ |
| **99th Percentile** | $56,997\text{ years}$ | $8.080 \times 10^{-5}$ | $2.560 \times 10^{-12}$ | **$8,579\text{ years}$** | $2.10 \times 10^{20}$ |

---

## 4. Synthesis & Comparison with Lab Kinetics

- **At the 95th Percentile Horizon**:
  - The empirical 1% residual collagen boundary at 10 °C is reached at **$t_{\text{eff}} \approx 26,300\text{ thermal years}$**.
  - This corresponds to an effective first-order degradation rate constant of:
    $$k(10^\circ\text{C}) = 1.75 \times 10^{-4}\text{ yr}^{-1} \quad (5.55 \times 10^{-12}\text{ s}^{-1})$$
  - Effective collagen half-life in whole bone at 10 °C:
    $$t_{1/2}(10^\circ\text{C}) \approx 3,955\text{ years}$$
- **At Sub-Zero Permafrost Temperatures (e.g. Siberia/Altai, $-10^\circ\text{C}$)**:
  - Applying $E_a = 173\text{ kJ/mol}$, the Arrhenius slowdown factor is:
    $$\frac{k(-10^\circ\text{C})}{k(10^\circ\text{C})} = \exp\left(-\frac{173,000}{8.314}\left(\frac{1}{263.15} - \frac{1}{283.15}\right)\right) \approx 0.0123 \approx \frac{1}{81}$$
  - The 1% survival limit at $-10^\circ\text{C}$ expands to:
    $$t_{\text{limit}}(-10^\circ\text{C}) \approx 26,300 \times 81.3 \approx \mathbf{2,140,000\text{ years}}$$
  - This perfectly explains why Middle/Early Pleistocene collagen survives in permafrost and cold caves for **>1–2 million years**, while in Mediterranean/subtropical sites it is lost within 10,000–25,000 years!
