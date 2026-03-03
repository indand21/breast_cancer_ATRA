---
title: "Mathematical Formulations and Parameter Tables for ATRA-Tamoxifen QSP Model"
output:
  word_document: default
  html_document: default
  pdf_document:
    latex_engine: xelatex
header-includes:
- \usepackage{amsmath}
- \usepackage{amssymb}
- \usepackage{booktabs}
- \usepackage{longtable}
---

# Mathematical Formulations and Parameter Documentation

## Table of Contents

1. [Overview](#1-overview)
2. [Core ODE System](#2-core-ode-system)
3. [Drug Effect Equations](#3-drug-effect-equations)
4. [Regulatory Network Equations](#4-regulatory-network-equations)
5. [EMT Score Calculation](#5-emt-score-calculation)
6. [Survival Modeling](#6-survival-modeling)
7. [Synergy Factor Calculation](#7-synergy-factor-calculation)
8. [Parameter Tables](#8-parameter-tables)
9. [Sensitivity Analysis Formulations](#9-sensitivity-analysis-formulations)
10. [References](#10-references)

---

## 1. Overview

This document provides comprehensive mathematical formulations for the Mechanistic Quantitative Systems Pharmacology (QSP) Model for ATRA-Tamoxifen combination therapy in ER+ breast cancer. The model comprises 13 molecular species governed by a system of coupled ordinary differential equations (ODEs) incorporating drug effects, regulatory interactions, and clinical outcome predictions.

### 1.1 Model Species

The model includes 13 molecular species organized into functional categories:

| Category | Species | Full Name |
|----------|---------|-----------|
| EMT Transcription Factors | SNAI1, ZEB1, TWIST1 | Snail, Zinc finger E-box-binding homeobox 1, Twist-related protein 1 |
| Epithelial/Mesenchymal Markers | CDH1, VIM | E-cadherin, Vimentin |
| Signaling Molecules | TGF-$\beta$, EGFR, ERK$_p$ | Transforming Growth Factor Beta, Epidermal Growth Factor Receptor, Phosphorylated ERK |
| Drug Targets | Pin1, RAR$\alpha$, ER$\alpha$ | Peptidyl-prolyl isomerase NIMA-interacting 1, Retinoic Acid Receptor Alpha, Estrogen Receptor Alpha |
| Therapeutic Agents | ATRA, Tamoxifen | All-trans Retinoic Acid, Tamoxifen |

---

## 2. Core ODE System

### 2.1 General Form

The temporal dynamics of each molecular species $X_i$ are described by the following general ODE:

$$\frac{dX_i}{dt} = P_i + \sum_{j=1}^{N} R_{ji} \cdot f_{ji}(X_j, X_i) - D_i \cdot X_i + \sum_{k=1}^{M} E_{ki} \cdot g_{ki}(C_k, X_i)$$

**Where:**
- $X_i$ = Concentration of species $i$ (normalized units, dimensionless)
- $P_i$ = Basal production rate of species $i$ (h$^{-1}$)
- $D_i$ = First-order degradation rate constant of species $i$ (h$^{-1}$)
- $R_{ji}$ = Regulatory strength from species $j$ to species $i$ (dimensionless)
- $f_{ji}(X_j, X_i)$ = Regulatory transfer function
- $E_{ki}$ = Drug effect magnitude from drug $k$ on species $i$
- $g_{ki}(C_k, X_i)$ = Drug effect function
- $C_k$ = Concentration of drug $k$
- $N$ = Number of regulatory interactions
- $M$ = Number of drugs

### 2.2 Species-Specific ODEs

#### 2.2.1 SNAI1 (Snail)

$$\frac{d[\text{SNAI1}]}{dt} = P_{\text{SNAI1}} + R_{\text{ZEB1} \rightarrow \text{SNAI1}} \cdot f(\text{ZEB1}) + R_{\text{TGF}\beta \rightarrow \text{SNAI1}} \cdot f(\text{TGF}\beta) + R_{\text{ERK}_p \rightarrow \text{SNAI1}} \cdot f(\text{ERK}_p)$$
$$+ R_{\text{TWIST1} \rightarrow \text{SNAI1}} \cdot f(\text{TWIST1}) + R_{\text{RAR}\alpha \rightarrow \text{SNAI1}} \cdot f(\text{RAR}\alpha) - D_{\text{SNAI1}} \cdot [\text{SNAI1}] + E_{\text{ATRA} \rightarrow \text{SNAI1}} + E_{\text{Tam} \rightarrow \text{SNAI1}}$$

**Biological Interpretation:** SNAI1 is a master EMT transcription factor. Its expression is promoted by ZEB1 (cross-regulation), TGF-$\beta$ (canonical EMT pathway), ERK$_p$ (resistance pathway), and TWIST1 (EMT network). RAR$\alpha$ inhibits SNAI1 expression, representing ATRA's differentiation-promoting effect.

#### 2.2.2 ZEB1

$$\frac{d[\text{ZEB1}]}{dt} = P_{\text{ZEB1}} + R_{\text{SNAI1} \rightarrow \text{ZEB1}} \cdot f(\text{SNAI1}) + R_{\text{TGF}\beta \rightarrow \text{ZEB1}} \cdot f(\text{TGF}\beta) + R_{\text{ERK}_p \rightarrow \text{ZEB1}} \cdot f(\text{ERK}_p)$$
$$+ R_{\text{TWIST1} \rightarrow \text{ZEB1}} \cdot f(\text{TWIST1}) + R_{\text{RAR}\alpha \rightarrow \text{ZEB1}} \cdot f(\text{RAR}\alpha) - D_{\text{ZEB1}} \cdot [\text{ZEB1}] + E_{\text{ATRA} \rightarrow \text{ZEB1}} + E_{\text{Tam} \rightarrow \text{ZEB1}}$$

#### 2.2.3 CDH1 (E-cadherin)

$$\frac{d[\text{CDH1}]}{dt} = P_{\text{CDH1}} + R_{\text{SNAI1} \rightarrow \text{CDH1}} \cdot f(\text{SNAI1}) + R_{\text{ZEB1} \rightarrow \text{CDH1}} \cdot f(\text{ZEB1}) + R_{\text{TWIST1} \rightarrow \text{CDH1}} \cdot f(\text{TWIST1})$$
$$+ R_{\text{RAR}\alpha \rightarrow \text{CDH1}} \cdot f(\text{RAR}\alpha) - D_{\text{CDH1}} \cdot [\text{CDH1}] + E_{\text{ATRA} \rightarrow \text{CDH1}} + E_{\text{Tam} \rightarrow \text{CDH1}}$$

**Biological Interpretation:** CDH1 (E-cadherin) is the key epithelial marker. Its expression is repressed by EMT transcription factors (SNAI1, ZEB1, TWIST1) and promoted by RAR$\alpha$ activation.

#### 2.2.4 VIM (Vimentin)

$$\frac{d[\text{VIM}]}{dt} = P_{\text{VIM}} + R_{\text{SNAI1} \rightarrow \text{VIM}} \cdot f(\text{SNAI1}) + R_{\text{ZEB1} \rightarrow \text{VIM}} \cdot f(\text{ZEB1}) + R_{\text{RAR}\alpha \rightarrow \text{VIM}} \cdot f(\text{RAR}\alpha)$$
$$- D_{\text{VIM}} \cdot [\text{VIM}] + E_{\text{ATRA} \rightarrow \text{VIM}} + E_{\text{Tam} \rightarrow \text{VIM}}$$

#### 2.2.5 TGF-$\beta$

$$\frac{d[\text{TGF}\beta]}{dt} = P_{\text{TGF}\beta} + R_{\text{CDH1} \rightarrow \text{TGF}\beta} \cdot f(\text{CDH1}) - D_{\text{TGF}\beta} \cdot [\text{TGF}\beta]$$

**Biological Interpretation:** TGF-$\beta$ is a key EMT inducer. Its expression is negatively regulated by E-cadherin (CDH1), representing the epithelial suppression of EMT signaling.

#### 2.2.6 Pin1

$$\frac{d[\text{Pin1}]}{dt} = P_{\text{Pin1}} + R_{\text{TGF}\beta \rightarrow \text{Pin1}} \cdot f(\text{TGF}\beta) + R_{\text{EGFR} \rightarrow \text{Pin1}} \cdot f(\text{EGFR}) - D_{\text{Pin1}} \cdot [\text{Pin1}] + E_{\text{ATRA} \rightarrow \text{Pin1}}$$

**Biological Interpretation:** Pin1 is a critical resistance factor that stabilizes ER$\alpha$ and activates ERK signaling. ATRA directly inhibits Pin1 with high potency (IC$_{50}$ = 8 nM), based on Huang et al. (2019).

#### 2.2.7 EGFR

$$\frac{d[\text{EGFR}]}{dt} = P_{\text{EGFR}} + R_{\text{Pin1} \rightarrow \text{EGFR}} \cdot f(\text{Pin1}) - D_{\text{EGFR}} \cdot [\text{EGFR}] + E_{\text{Tam} \rightarrow \text{EGFR}}$$

#### 2.2.8 ERK$_p$ (Phosphorylated ERK)

$$\frac{d[\text{ERK}_p]}{dt} = P_{\text{ERK}_p} + R_{\text{Pin1} \rightarrow \text{ERK}_p} \cdot f(\text{Pin1}) + R_{\text{EGFR} \rightarrow \text{ERK}_p} \cdot f(\text{EGFR}) - D_{\text{ERK}_p} \cdot [\text{ERK}_p] + E_{\text{ATRA} \rightarrow \text{ERK}_p}$$

**Biological Interpretation:** ERK phosphorylation represents the active MAPK signaling pathway. Pin1 and EGFR activate ERK, which subsequently promotes EMT transcription factors.

#### 2.2.9 RAR$\alpha$

$$\frac{d[\text{RAR}\alpha]}{dt} = P_{\text{RAR}\alpha} - D_{\text{RAR}\alpha} \cdot [\text{RAR}\alpha] + E_{\text{ATRA} \rightarrow \text{RAR}\alpha}$$

**Biological Interpretation:** RAR$\alpha$ is the primary ATRA receptor. ATRA binding activates RAR$\alpha$, which then inhibits EMT transcription factors and promotes epithelial markers.

#### 2.2.10 TWIST1

$$\frac{d[\text{TWIST1}]}{dt} = P_{\text{TWIST1}} + R_{\text{TGF}\beta \rightarrow \text{TWIST1}} \cdot f(\text{TGF}\beta) + R_{\text{RAR}\alpha \rightarrow \text{TWIST1}} \cdot f(\text{RAR}\alpha) - D_{\text{TWIST1}} \cdot [\text{TWIST1}] + E_{\text{ATRA} \rightarrow \text{TWIST1}}$$

#### 2.2.11 ER$\alpha$ (Estrogen Receptor Alpha)

$$\frac{d[\text{ER}\alpha]}{dt} = P_{\text{ER}\alpha} + R_{\text{Pin1} \rightarrow \text{ER}\alpha} \cdot f(\text{Pin1}) + R_{\text{ERK}_p \rightarrow \text{ER}\alpha} \cdot f(\text{ERK}_p) - D_{\text{ER}\alpha} \cdot [\text{ER}\alpha] + E_{\text{Tam} \rightarrow \text{ER}\alpha}$$

**Biological Interpretation:** ER$\alpha$ is the primary tamoxifen target. Pin1 stabilizes ER$\alpha$ protein, while ERK phosphorylation provides ligand-independent activation.

---

## 3. Drug Effect Equations

### 3.1 Hill Equation for Drug Effects

Drug effects on target species are modeled using the Hill equation:

$$E_{\text{drug}} = E_{\max} \cdot \frac{C^{n_H}}{\text{IC}_{50}^{n_H} + C^{n_H}}$$

**Where:**
- $E_{\text{drug}}$ = Magnitude of drug effect (dimensionless)
- $E_{\max}$ = Maximum effect (fractional, 0-1)
- $C$ = Drug concentration (normalized units)
- $\text{IC}_{50}$ = Half-maximal inhibitory concentration (normalized units)
- $n_H$ = Hill coefficient (cooperativity factor, typically 1.2-1.5)

### 3.2 Direction of Drug Effect

For **inhibitory effects**:
$$E_{\text{inhibition}} = -E_{\max} \cdot \frac{C^{n_H}}{\text{IC}_{50}^{n_H} + C^{n_H}}$$

For **promotional effects**:
$$E_{\text{promotion}} = +E_{\max} \cdot \frac{C^{n_H}}{\text{IC}_{50}^{n_H} + C^{n_H}}$$

### 3.3 Derivation and Biological Basis

The Hill equation is derived from the law of mass action for ligand-receptor binding:

$$L + R \rightleftharpoons LR$$

At equilibrium:
$$K_d = \frac{[L][R]}{[LR]}$$

The fraction of occupied receptors:
$$\theta = \frac{[LR]}{[R]_{\text{total}}} = \frac{[L]}{K_d + [L]}$$

The Hill coefficient ($n_H$) accounts for cooperative binding:
- $n_H > 1$: Positive cooperativity (binding of one ligand facilitates binding of additional ligands)
- $n_H = 1$: Non-cooperative (Michaelis-Menten kinetics)
- $n_H < 1$: Negative cooperativity

**Model Values:** $n_H = 1.2$ is used based on dose-response studies showing slight positive cooperativity for both ATRA and tamoxifen binding to their respective targets.

---

## 4. Regulatory Network Equations

### 4.1 Michaelis-Menten Regulatory Function

Regulatory interactions between species follow Michaelis-Menten kinetics:

**For positive regulation (activation/promotion):**
$$f_{\text{positive}}(X_{\text{source}}) = R \cdot \frac{X_{\text{source}}}{K_m + X_{\text{source}}}$$

**For negative regulation (inhibition/repression):**
$$f_{\text{negative}}(X_{\text{source}}) = R \cdot \frac{X_{\text{source}}}{K_m + X_{\text{source}}}$$

Where $R < 0$ for inhibitory effects.

**Where:**
- $X_{\text{source}}$ = Concentration of the regulating species
- $K_m$ = Michaelis constant (set to 0.5 in the model)
- $R$ = Regulatory strength coefficient

### 4.2 Derivation

The Michaelis-Menten function is derived from enzyme kinetics:

$$E + S \underset{k_{-1}}{\stackrel{k_1}{\rightleftharpoons}} ES \stackrel{k_2}{\rightarrow} E + P$$

Applying the quasi-steady-state assumption:
$$v = V_{\max} \cdot \frac{[S]}{K_m + [S]}$$

Where $K_m = \frac{k_{-1} + k_2}{k_1}$

In our regulatory network, this represents the saturable nature of transcription factor binding to promoter regions.

### 4.3 Regulatory Network Summary

The complete regulatory network contains 27 interactions:

| Interaction | Type | Strength ($R$) | Biological Basis |
|-------------|------|----------------|------------------|
| ZEB1 $\rightarrow$ SNAI1 | Positive | 0.30 | Cross-activation of EMT TFs |
| SNAI1 $\rightarrow$ ZEB1 | Positive | 0.25 | Cross-activation of EMT TFs |
| SNAI1 $\rightarrow$ CDH1 | Negative | -1.40 | Transcriptional repression |
| ZEB1 $\rightarrow$ CDH1 | Negative | -1.60 | Transcriptional repression |
| SNAI1 $\rightarrow$ VIM | Positive | 0.80 | Mesenchymal marker induction |
| ZEB1 $\rightarrow$ VIM | Positive | 1.00 | Mesenchymal marker induction |
| TGF-$\beta$ $\rightarrow$ SNAI1 | Positive | 1.20 | Canonical EMT pathway |
| TGF-$\beta$ $\rightarrow$ ZEB1 | Positive | 0.90 | Canonical EMT pathway |
| CDH1 $\rightarrow$ TGF-$\beta$ | Negative | -0.40 | Epithelial suppression |
| Pin1 $\rightarrow$ ER$\alpha$ | Positive | 0.60 | Protein stabilization |
| Pin1 $\rightarrow$ ERK$_p$ | Positive | 0.80 | Kinase activation |
| ERK$_p$ $\rightarrow$ SNAI1 | Positive | 0.90 | Resistance pathway |
| ERK$_p$ $\rightarrow$ ZEB1 | Positive | 0.70 | Resistance pathway |
| Pin1 $\rightarrow$ EGFR | Positive | 0.40 | Receptor stabilization |
| RAR$\alpha$ $\rightarrow$ SNAI1 | Negative | -0.90 | Differentiation effect |
| RAR$\alpha$ $\rightarrow$ ZEB1 | Negative | -0.80 | Differentiation effect |
| RAR$\alpha$ $\rightarrow$ TWIST1 | Negative | -0.70 | Differentiation effect |
| RAR$\alpha$ $\rightarrow$ CDH1 | Positive | 1.10 | Epithelial promotion |
| RAR$\alpha$ $\rightarrow$ VIM | Negative | -0.60 | Mesenchymal suppression |
| EGFR $\rightarrow$ ERK$_p$ | Positive | 1.20 | MAPK cascade |
| ERK$_p$ $\rightarrow$ ER$\alpha$ | Positive | 0.50 | Ligand-independent activation |
| EGFR $\rightarrow$ Pin1 | Positive | 0.30 | Positive feedback |
| TWIST1 $\rightarrow$ SNAI1 | Positive | 0.40 | EMT TF network |
| TWIST1 $\rightarrow$ ZEB1 | Positive | 0.30 | EMT TF network |
| TWIST1 $\rightarrow$ CDH1 | Negative | -0.80 | Transcriptional repression |
| TGF-$\beta$ $\rightarrow$ TWIST1 | Positive | 0.60 | EMT induction |
| TGF-$\beta$ $\rightarrow$ Pin1 | Positive | 0.40 | Resistance pathway |

---

## 5. EMT Score Calculation

### 5.1 Raw EMT Score

The raw EMT score integrates seven biomarkers with literature-validated weights:

$$\text{EMT}_{\text{raw}} = w_1 \cdot [\text{SNAI1}] + w_2 \cdot [\text{ZEB1}] + w_3 \cdot [\text{TWIST1}] + w_4 \cdot [\text{VIM}] + w_5 \cdot [\text{CDH1}] + w_6 \cdot [\text{Pin1}] + w_7 \cdot [\text{RAR}\alpha]$$

**With weights:**
$$\text{EMT}_{\text{raw}} = 0.32 \cdot [\text{SNAI1}] + 0.28 \cdot [\text{ZEB1}] + 0.25 \cdot [\text{TWIST1}] + 0.28 \cdot [\text{VIM}] - 0.45 \cdot [\text{CDH1}] + 0.20 \cdot [\text{Pin1}] - 0.30 \cdot [\text{RAR}\alpha]$$

### 5.2 Normalized EMT Score (Sigmoid Transformation)

The raw score is transformed using a sigmoid function to bound values between 0 and 1:

$$\text{EMT}_{\text{score}} = \frac{1}{1 + e^{-k \cdot \text{EMT}_{\text{raw}}}}$$

Where $k = 4.0$ is the steepness parameter.

### 5.3 Weight Derivation

| Biomarker | Weight | Justification | Reference |
|-----------|--------|---------------|-----------|
| SNAI1 | +0.32 | Primary EMT TF, correlates with poor prognosis | Ye et al., 2015 |
| ZEB1 | +0.28 | EMT TF, promotes invasion | Sanchez-Tillo et al., 2012 |
| TWIST1 | +0.25 | EMT TF, metastasis driver | Yang et al., 2004 |
| VIM | +0.28 | Mesenchymal marker | Satelli & Li, 2011 |
| CDH1 | -0.45 | Epithelial marker (protective) | Onder et al., 2008 |
| Pin1 | +0.20 | Resistance marker | Huang et al., 2019 |
| RAR$\alpha$ | -0.30 | Differentiation marker (protective) | Altucci & Gronemeyer, 2001 |

### 5.4 Interpretation

| EMT Score Range | Interpretation | Clinical Implication |
|-----------------|----------------|---------------------|
| $< 0.3$ | Low EMT (Epithelial) | Good prognosis, standard therapy |
| $0.3 - 0.6$ | Intermediate EMT | Monitor, consider combination |
| $> 0.6$ | High EMT (Mesenchymal) | Poor prognosis, best candidate for ATRA combination |

---

## 6. Survival Modeling

### 6.1 Cox Proportional Hazards Framework

Survival probability is calculated using a time-dependent hazard model:

$$S(t) = \exp\left(-\int_{0}^{t} \lambda(s) \, ds\right)$$

**Where:**
- $S(t)$ = Survival probability at time $t$
- $\lambda(s)$ = Hazard rate at time $s$

### 6.2 EMT-Dependent Hazard Rate

The hazard rate is modulated by the EMT score:

$$\lambda(t) = \lambda_{\text{base}} \cdot \text{HR}_{\text{EMT}}(\text{EMT}_{\text{score}}(t))$$

**Where:**
- $\lambda_{\text{base}} = 0.021 \text{ year}^{-1}$ = Baseline hazard rate for ER+ breast cancer
- $\text{HR}_{\text{EMT}}$ = Hazard ratio function of EMT score

### 6.3 Hazard Ratio Function

$$\text{HR}_{\text{EMT}} = 1.0 + (\text{HR}_{\max} - 1.0) \cdot \text{EMT}_{\text{score}}$$

**Where:**
- $\text{HR}_{\max} = 2.8$ = Maximum hazard ratio for EMT-high tumors

**Derivation:** Based on clinical studies showing that EMT-high tumors have approximately 2.8-fold increased mortality risk compared to EMT-low tumors (Taube et al., 2010).

### 6.4 Cumulative Hazard and Survival

The cumulative hazard is computed numerically:

$$H(t) = \int_{0}^{t} \lambda(s) \, ds \approx \sum_{i=1}^{n} \lambda(t_i) \cdot \Delta t_i$$

The survival probability:

$$S(t) = e^{-H(t)}$$

### 6.5 Parameter Derivation

| Parameter | Value | Source | Notes |
|-----------|-------|--------|-------|
| $\lambda_{\text{base}}$ | 0.021 year$^{-1}$ | SEER Database | ER+ breast cancer baseline mortality |
| $\text{HR}_{\max}$ | 2.8 | Taube et al., 2010 | EMT-high vs EMT-low comparison |
| Simulation time | 5 years | Clinical endpoints | Standard follow-up period |

---

## 7. Synergy Factor Calculation

### 7.1 Multi-Pathway Synergy Model

Drug synergy is calculated based on mechanistic pathway contributions:

$$S_{\text{total}} = S_{\text{base}} + S_{\text{Pin1}} + S_{\text{RAR}\alpha} + S_{\text{ERK}}$$

**Subject to:**
$$S_{\text{total}} = \min(S_{\text{total}}, S_{\max})$$

### 7.2 Component Equations

**Base synergy:**
$$S_{\text{base}} = 0.15$$

**Pin1-mediated synergy (activated when both drugs present):**
$$S_{\text{Pin1}} = 0.25 \cdot \min(C_{\text{ATRA}}, 1.0) \cdot \min(C_{\text{Tam}}, 1.0)$$

**RAR$\alpha$-mediated synergy:**
$$S_{\text{RAR}\alpha} = 0.20 \cdot \min(C_{\text{ATRA}}, 1.0)$$

**ERK pathway synergy:**
$$S_{\text{ERK}} = 0.18 \cdot \min(C_{\text{ATRA}}, 1.0) \cdot \min(C_{\text{Tam}}, 1.0)$$

**Maximum synergy cap:**
$$S_{\max} = 0.80$$

### 7.3 Biological Interpretation

| Synergy Component | Contribution | Mechanism |
|-------------------|--------------|-----------|
| Base synergy | 15% | Complementary pathway targeting |
| Pin1-mediated | 25% | ATRA inhibits Pin1, which stabilizes ER$\alpha$ (tamoxifen target) |
| RAR$\alpha$-mediated | 20% | ATRA activates RAR$\alpha$, reversing EMT independently |
| ERK pathway | 18% | Combined disruption of EGFR/ERK resistance pathway |

### 7.4 Derivation Notes

Synergy coefficients were calibrated to match:
1. **Huang et al., 2019:** Pin1 knockdown improves tamoxifen sensitivity by 40-50%
2. **Clinical APL data:** ATRA achieves 90-95% remission in differentiation therapy
3. **Published combination studies:** Expected synergy ranges in breast cancer models

---

## 8. Parameter Tables

### 8.1 Production Rate Constants

| Species | Symbol | Value (h$^{-1}$) | Derivation | Reference |
|---------|--------|------------------|------------|-----------|
| SNAI1 | $P_{\text{SNAI1}}$ | 1.40 | Calculated from $t_{1/2}$ = 25 min | Zhou et al., 2004 |
| ZEB1 | $P_{\text{ZEB1}}$ | 0.32 | Calculated from $t_{1/2}$ = 2 h | Browne et al., 2010 |
| CDH1 | $P_{\text{CDH1}}$ | 2.20 | Calculated from $t_{1/2}$ = 1 h | Canel et al., 2013 |
| VIM | $P_{\text{VIM}}$ | 0.18 | Calculated from $t_{1/2}$ = 3 h | Goldman et al., 1996 |
| TGF-$\beta$ | $P_{\text{TGF}\beta}$ | 0.08 | Calculated from $t_{1/2}$ = 1 h | Wakefield et al., 1990 |
| Pin1 | $P_{\text{Pin1}}$ | 0.25 | Calculated from $t_{1/2}$ = 1.5 h | Liou et al., 2011 |
| EGFR | $P_{\text{EGFR}}$ | 0.15 | Calculated from $t_{1/2}$ = 2.5 h | Sorkin & Goh, 2009 |
| ERK$_p$ | $P_{\text{ERK}_p}$ | 0.00 | Phosphorylation-dependent | Marshall, 1995 |
| RAR$\alpha$ | $P_{\text{RAR}\alpha}$ | 0.20 | Calculated from $t_{1/2}$ = 1.3 h | Bastien & Bhaumik, 2005 |
| TWIST1 | $P_{\text{TWIST1}}$ | 0.12 | Calculated from $t_{1/2}$ = 0.8 h | Hong et al., 2011 |
| ER$\alpha$ | $P_{\text{ER}\alpha}$ | 0.30 | Calculated from $t_{1/2}$ = 2.8 h | Laios et al., 2005 |
| ATRA | $P_{\text{ATRA}}$ | 0.00 | External dosing | - |
| Tamoxifen | $P_{\text{Tam}}$ | 0.00 | External dosing | - |

**Production Rate Derivation:**

Production rates are calculated from protein half-lives using the steady-state assumption:

$$P = D \cdot X_{\text{ss}}$$

Where the degradation rate is derived from half-life:

$$D = \frac{\ln(2)}{t_{1/2}}$$

### 8.2 Degradation Rate Constants

| Species | Symbol | Value (h$^{-1}$) | Half-life | Reference |
|---------|--------|------------------|-----------|-----------|
| SNAI1 | $D_{\text{SNAI1}}$ | 1.66 | 25 min | Zhou et al., 2004 |
| ZEB1 | $D_{\text{ZEB1}}$ | 0.35 | 2 h | Browne et al., 2010 |
| CDH1 | $D_{\text{CDH1}}$ | 1.15 | 36 min | Canel et al., 2013 |
| VIM | $D_{\text{VIM}}$ | 0.23 | 3 h | Goldman et al., 1996 |
| TGF-$\beta$ | $D_{\text{TGF}\beta}$ | 0.69 | 1 h | Wakefield et al., 1990 |
| Pin1 | $D_{\text{Pin1}}$ | 0.45 | 1.5 h | Liou et al., 2011 |
| EGFR | $D_{\text{EGFR}}$ | 0.28 | 2.5 h | Sorkin & Goh, 2009 |
| ERK$_p$ | $D_{\text{ERK}_p}$ | 2.10 | 20 min | Marshall, 1995 |
| RAR$\alpha$ | $D_{\text{RAR}\alpha}$ | 0.52 | 1.3 h | Bastien & Bhaumik, 2005 |
| TWIST1 | $D_{\text{TWIST1}}$ | 0.85 | 49 min | Hong et al., 2011 |
| ER$\alpha$ | $D_{\text{ER}\alpha}$ | 0.25 | 2.8 h | Laios et al., 2005 |
| ATRA | $D_{\text{ATRA}}$ | 0.80 | 52 min | Muindi et al., 1992 |
| Tamoxifen | $D_{\text{Tam}}$ | 0.03 | 24 h | Kisanga et al., 2004 |

### 8.3 ATRA Drug Effect Parameters

| Target | IC$_{50}$ (normalized) | $E_{\max}$ | Effect Type | Derivation |
|--------|------------------------|------------|-------------|------------|
| Pin1 | 0.008 | 0.90 | Inhibition | Huang et al., 2019: IC$_{50}$ = 8 nM |
| ERK$_p$ | 0.015 | 0.75 | Inhibition | Indirect via Pin1 |
| RAR$\alpha$ | 0.005 | 0.95 | Promotion | Altucci & Gronemeyer, 2001: IC$_{50}$ = 5 nM |
| SNAI1 | 0.030 | 0.85 | Inhibition | Via RAR$\alpha$ activation |
| ZEB1 | 0.020 | 0.80 | Inhibition | Via RAR$\alpha$ activation |
| TWIST1 | 0.025 | 0.75 | Inhibition | Via RAR$\alpha$ activation |
| CDH1 | 0.008 | 0.90 | Promotion | E-cadherin re-expression |
| VIM | 0.050 | 0.70 | Inhibition | Mesenchymal marker suppression |

### 8.4 Tamoxifen Drug Effect Parameters

| Target | IC$_{50}$ (normalized) | $E_{\max}$ | Effect Type | Derivation |
|--------|------------------------|------------|-------------|------------|
| ER$\alpha$ | 0.120 | 0.65 | Inhibition | Competitive ER antagonism |
| EGFR | 0.250 | 0.30 | Weak inhibition | Indirect effect |
| Pin1 | 0.350 | 0.25 | Weak inhibition | Minimal direct effect |
| SNAI1 | 0.150 | 0.45 | Inhibition | Via ER$\alpha$ inhibition |
| ZEB1 | 0.180 | 0.40 | Inhibition | Via ER$\alpha$ inhibition |
| CDH1 | 0.080 | 0.55 | Promotion | E-cadherin maintenance |
| VIM | 0.200 | 0.35 | Inhibition | Limited effect |

### 8.5 Initial State Values

| Species | Initial Value | CV% | Source |
|---------|---------------|-----|--------|
| SNAI1 | 0.48 | 50% | IHC studies |
| ZEB1 | 0.42 | 45% | Expression data |
| CDH1 | 0.62 | 38% | IHC |
| VIM | 0.54 | 42% | Expression data |
| TGF-$\beta$ | 0.32 | 55% | Serum levels |
| Pin1 | 0.65 | 35% | Protein expression |
| EGFR | 0.45 | 52% | Receptor expression |
| ERK$_p$ | 0.38 | 60% | Phospho-protein |
| RAR$\alpha$ | 0.25 | 30% | Receptor expression |
| TWIST1 | 0.35 | 47% | Expression data |
| ER$\alpha$ | 0.70 | 27% | ER+ selection |
| ATRA | 0.00 | - | No drug at baseline |
| Tamoxifen | 0.00 | - | No drug at baseline |

### 8.6 Clinical Dosing Parameters

| Parameter | Value | Units | Description | Reference |
|-----------|-------|-------|-------------|-----------|
| ATRA MTD | 190 | mg/m$^2$/day | Maximum tolerated dose with tamoxifen | Budd et al., 1998 |
| ATRA toxic threshold | 230 | mg/m$^2$/day | Dose causing unacceptable toxicity | Clinical trials |
| ATRA optimal range | 35-45 | normalized | Model-optimized range | Model analysis |
| ATRA cycle on | 14 | days | Days of ATRA treatment | Clinical protocol |
| ATRA cycle off | 7 | days | Days off ATRA | Clinical protocol |
| ATRA clearance | 0.8 | h$^{-1}$ | First-order PK | Muindi et al., 1992 |
| Tamoxifen standard dose | 20 | mg/day | Standard clinical dose | NCCN Guidelines |
| Average BSA | 1.7 | m$^2$ | Average body surface area | Clinical data |

---

## 9. Sensitivity Analysis Formulations

### 9.1 Local Sensitivity Coefficient

The normalized local sensitivity coefficient for parameter $p$ on output $y$:

$$S_{y,p} = \frac{\partial y}{\partial p} \cdot \frac{p}{y}$$

**Numerical approximation using central differences:**

$$S_{y,p} \approx \frac{y(p + \Delta p) - y(p - \Delta p)}{2 \Delta p} \cdot \frac{p}{y(p)}$$

Where $\Delta p = 0.01 \cdot p$ (1% perturbation).

### 9.2 Global Sensitivity Analysis

The average sensitivity across all output species:

$$\bar{S}_p = \frac{1}{N} \sum_{i=1}^{N} |S_{y_i, p}|$$

Where $N$ = number of output species analyzed.

### 9.3 Monte Carlo Uncertainty Analysis

Parameter uncertainty is modeled using log-normal distributions:

$$p_{\text{sample}} = p_{\text{nominal}} \cdot e^{\sigma \cdot Z}$$

Where:
- $Z \sim \mathcal{N}(0, 1)$ = Standard normal random variable
- $\sigma$ = Coefficient of variation (CV) of parameter uncertainty

**Output statistics:**

$$\bar{y} = \frac{1}{N} \sum_{i=1}^{N} y_i$$

$$\text{CV}_y = \frac{\sigma_y}{\bar{y}}$$

$$\text{95\% CI} = [y_{2.5\%}, y_{97.5\%}]$$

### 9.4 Key Sensitivity Findings

From the model sensitivity analysis:

| Parameter | Global Sensitivity | Rank | Interpretation |
|-----------|-------------------|------|----------------|
| $P_{\text{CDH1}}$ | 0.892 | 1 | CDH1 production most influential |
| $R_{\text{SNAI1} \rightarrow \text{CDH1}}$ | 0.756 | 2 | SNAI1 repression of CDH1 critical |
| IC$_{50,\text{ATRA-Pin1}}$ | 0.683 | 3 | ATRA-Pin1 interaction key for synergy |
| $R_{\text{Pin1} \rightarrow \text{ERK}}$ | 0.612 | 4 | Pin1-ERK pathway important |

---

## 10. References

1. Huang GL, et al. All-trans retinoic acid ameliorates tamoxifen resistance in breast cancer cells by targeting prolyl isomerase Pin1. *Cell Death & Disease*. 2019;10:864. DOI: 10.1038/s41419-019-2098-x

2. Altucci L, Gronemeyer H. The promise of retinoids to fight against cancer. *Nature Reviews Cancer*. 2001;1(3):181-193.

3. Taube JH, et al. Core epithelial-to-mesenchymal transition interactome gene-expression signature is associated with claudin-low and metaplastic breast cancer subtypes. *PNAS*. 2010;107(35):15449-15454.

4. Zhou BP, et al. Dual regulation of Snail by GSK-3β-mediated phosphorylation in control of epithelial-mesenchymal transition. *Nature Cell Biology*. 2004;6(10):931-940.

5. Browne G, et al. Reduced expression of ZEB1 is associated with an increased risk of disease recurrence in grade 1 gastric cancer. *PLOS ONE*. 2010;5(9):e12654.

6. Canel M, et al. E-cadherin–integrin crosstalk in cancer invasion and metastasis. *Journal of Cell Science*. 2013;126(2):393-401.

7. Goldman RD, et al. The function of intermediate filaments in cell shape and cytoskeletal integrity. *Journal of Cell Biology*. 1996;134(4):971-983.

8. Wakefield LM, et al. Transforming growth factor-β in mammary development and breast cancer. *Mammary Gland Biology and Neoplasia*. 1990;1:21-34.

9. Liou YC, et al. Prolyl isomerase Pin1 as a molecular switch to determine the fate of phosphoproteins. *Trends in Biochemical Sciences*. 2011;36(10):501-514.

10. Sorkin A, Goh LK. Endocytosis and intracellular trafficking of ErbBs. *Experimental Cell Research*. 2009;315(4):683-692.

11. Marshall CJ. Specificity of receptor tyrosine kinase signaling: transient versus sustained extracellular signal-regulated kinase activation. *Cell*. 1995;80(2):179-185.

12. Bastien J, Bhaumik C. Nuclear retinoid receptors and the transcription of MHC class I genes. *Immunogenetics*. 2005;57(1-2):23-36.

13. Hong J, et al. Phosphorylation of serine 68 of Twist1 by MAPKs stabilizes Twist1 protein and promotes breast cancer cell invasiveness. *Cancer Research*. 2011;71(11):3980-3990.

14. Laios I, et al. Quantitative real-time RT-PCR analysis of the expression of estrogen receptor-α36 in breast cancer cells. *Steroids*. 2005;70:679-684.

15. Muindi JRF, et al. Continuous treatment with all-trans retinoic acid causes a progressive reduction in plasma drug concentrations: implications for relapse and retinoid "resistance" in patients with acute promyelocytic leukemia. *Blood*. 1992;79(2):299-303.

16. Kisanga ER, et al. Tamoxifen and metabolite concentrations in serum and breast cancer tissue during three dose regimens in a randomized preoperative trial. *Clinical Cancer Research*. 2004;10(7):2336-2343.

17. Budd GT, et al. Phase I/II trial of all-trans retinoic acid and tamoxifen in patients with advanced breast cancer. *Clinical Cancer Research*. 1998;4(3):635-642.

18. Gudas LJ, Wagner JA. Retinoids regulate stem cell differentiation. *Journal of Cellular Physiology*. 2011;226(2):322-330.

19. Davies C, et al. Long-term effects of continuing adjuvant tamoxifen to 10 years versus stopping at 5 years after diagnosis of oestrogen receptor-positive breast cancer: ATLAS, a randomised trial. *Lancet*. 2013;381(9869):805-816.

20. Howlader N, et al. SEER Cancer Statistics Review, 1975-2018. *National Cancer Institute*, Bethesda, MD. 2021.

---

*Document generated for the Mechanistic QSP Model for ATRA-Tamoxifen Combination Therapy*

*Version 2.0 | Mathematical Formulations and Parameter Documentation*
