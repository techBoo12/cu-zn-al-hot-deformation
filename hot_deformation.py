"""
hot_deformation.py
==================

Python implementation of the constitutive analysis carried out in the BEng
thesis "Assessment of hot deformability of Copper Zinc Aluminium based alloy
using hyperbolic sine constitutive model" (A. M. Ajufo, supervised by
Prof. K. K. Alaneme, Federal University of Technology, Akure, 2024).

The script

1. reads the peak flow stresses measured in Gleeble 3500 hot-compression
   tests (three temperatures x three strain rates);
2. determines the constants of the Sellars-Tegart hyperbolic-sine
   (Arrhenius-type) constitutive equation, n', beta, alpha, n, Q and A, by
   linear regression on the standard linearised plots;
3. back-calculates the peak stress from the fitted model and reports the
   coefficient of determination (R^2) and the average absolute relative
   error (AARE);
4. saves the fitting plots to ``figures/`` and the numerical results to
   ``results_constants.csv`` and ``results_predictions.csv``.

Usage
-----
    python hot_deformation.py
"""

import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")  # non-interactive backend: write PNG files only
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
R = 8.314  # universal gas constant, J mol^-1 K^-1

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_FILE = os.path.join(HERE, "data", "peak_stress.csv")
FIG_DIR = os.path.join(HERE, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 13,
    "legend.fontsize": 11,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
})

# One colour and marker per test temperature (and per strain rate in Fig. 4).
PALETTE = ["#1f4e79", "#c55a11", "#548235"]
MARKERS = ["o", "s", "^"]
DPI = 200

# Axis-label strings (matplotlib mathtext).
LBL_LN_RATE = r"$\ln\,\dot{\varepsilon}$  ($\dot{\varepsilon}$ in s$^{-1}$)"
LBL_LN_SINH = r"$\ln[\sinh(\alpha\sigma_p)]$"

# ---------------------------------------------------------------------------
# Experimental data
# ---------------------------------------------------------------------------
data = pd.read_csv(DATA_FILE)
temps_C = sorted(data["temperature_C"].unique())       # 250, 350, 450 degC
rates = sorted(data["strain_rate_per_s"].unique())     # 0.1, 1, 5 s^-1

# Peak stress matrix: rows = temperature, columns = strain rate (MPa).
stress = (data.pivot(index="temperature_C", columns="strain_rate_per_s",
                     values="peak_stress_MPa")
          .loc[temps_C, rates].to_numpy())
temps_K = np.array(temps_C) + 273.15
ln_rate = np.log(np.array(rates))

print("Peak stress (MPa); rows: temperature, columns: strain rate")
print(pd.DataFrame(stress, index=[f"{t} °C" for t in temps_C],
                   columns=[f"{r} s^-1" for r in rates]).round(1), "\n")


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def linear_fit(x, y):
    """Ordinary least-squares straight line through (x, y).

    Returns
    -------
    slope, intercept : float
    """
    slope, intercept = np.polyfit(x, y, 1)
    return slope, intercept


def format_axes(ax, xlabel, ylabel, title):
    """Apply a consistent style to a figure's axes."""
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.grid(alpha=0.3)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def save(fig, fname):
    """Save a figure to the figures directory and close it."""
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, fname), dpi=DPI)
    plt.close(fig)


def fit_per_temperature(x_grid, fname, xlabel, title, symbol):
    """Regress ln(strain rate) on x separately at each test temperature.

    Parameters
    ----------
    x_grid : ndarray, shape (n_temperatures, n_rates)
        Abscissa values, one row per temperature.
    fname : str
        Output file name for the plot.
    xlabel, title : str
        Axis label and title of the plot.
    symbol : str
        Mathtext symbol of the slope, used in the legend.

    Returns
    -------
    list of float
        Fitted slope at each temperature.
    """
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    slopes = []
    for i, T in enumerate(temps_C):
        x = x_grid[i]
        m, c = linear_fit(x, ln_rate)
        slopes.append(m)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(x, ln_rate, MARKERS[i], color=PALETTE[i], ms=8,
                label=f"{T} °C  ({symbol} = {m:.3g})")
        ax.plot(xs, m * xs + c, "-", color=PALETTE[i], lw=1.8)
    format_axes(ax, xlabel, LBL_LN_RATE, title)
    ax.legend(frameon=False)
    save(fig, fname)
    return slopes


# ---------------------------------------------------------------------------
# Step 1: power-law exponent n'  (low-stress regime)
#   ln(strain rate) = n' ln(sigma) + const
# ---------------------------------------------------------------------------
n1_each = fit_per_temperature(
    np.log(stress), "1_n1_ln_rate_vs_ln_stress.png",
    r"$\ln\,\sigma_p$  ($\sigma_p$ in MPa)",
    r"(a) Power law: determination of $n'$", r"$n'$")
n1 = np.mean(n1_each)

# ---------------------------------------------------------------------------
# Step 2: exponential-law constant beta  (high-stress regime)
#   ln(strain rate) = beta * sigma + const
# ---------------------------------------------------------------------------
beta_each = fit_per_temperature(
    stress, "2_beta_ln_rate_vs_stress.png",
    r"Peak stress, $\sigma_p$ (MPa)",
    r"(b) Exponential law: determination of $\beta$", r"$\beta$")
beta = np.mean(beta_each)

# ---------------------------------------------------------------------------
# Step 3: stress multiplier alpha = beta / n'  (MPa^-1)
# ---------------------------------------------------------------------------
alpha = beta / n1
ln_sinh = np.log(np.sinh(alpha * stress))

# ---------------------------------------------------------------------------
# Step 4: stress exponent n of the hyperbolic-sine law
#   ln(strain rate) = n ln[sinh(alpha sigma)] + const
# ---------------------------------------------------------------------------
n_each = fit_per_temperature(
    ln_sinh, "3_n_ln_rate_vs_ln_sinh.png", LBL_LN_SINH,
    r"(c) Hyperbolic-sine law: determination of $n$", r"$n$")
n = np.mean(n_each)

# ---------------------------------------------------------------------------
# Step 5: apparent activation energy Q
#   At constant strain rate: ln[sinh(alpha sigma)] = K (1000/T) + const
#   Q = 1000 R n K   (T in kelvin; K averaged over the strain rates)
# ---------------------------------------------------------------------------
inv_T = 1000.0 / temps_K
fig, ax = plt.subplots(figsize=(6.4, 4.8))
K_each = []
for j, rate in enumerate(rates):
    y = ln_sinh[:, j]
    K, c = linear_fit(inv_T, y)
    K_each.append(K)
    xs = np.linspace(inv_T.min(), inv_T.max(), 50)
    ax.plot(inv_T, y, MARKERS[j], color=PALETTE[j], ms=8,
            label=f"{rate:g} s$^{{-1}}$  (slope = {K:.2f})")
    ax.plot(xs, K * xs + c, "-", color=PALETTE[j], lw=1.8)
format_axes(ax, r"$1000/T$ (K$^{-1}$)", LBL_LN_SINH,
            r"(d) Determination of the activation energy $Q$")
ax.legend(frameon=False)
save(fig, "4_Q_ln_sinh_vs_inverse_T.png")
K_mean = np.mean(K_each)
Q = R * n * K_mean * 1000.0  # J mol^-1
Q_kJ = Q / 1000.0

# ---------------------------------------------------------------------------
# Step 6: Zener-Hollomon parameter and material constant A
#   Z = strain_rate * exp(Q / (R T))
#   ln Z = ln A + n ln[sinh(alpha sigma)]  ->  ln A is the intercept
# ---------------------------------------------------------------------------
ln_Z = ln_rate[None, :] + Q / (R * temps_K[:, None])
slope_Z, lnA = linear_fit(ln_sinh.ravel(), ln_Z.ravel())
A = np.exp(lnA)

fig, ax = plt.subplots(figsize=(6.4, 4.8))
for i, T in enumerate(temps_C):
    ax.plot(ln_sinh[i], ln_Z[i], MARKERS[i], color=PALETTE[i], ms=8,
            label=f"{T} °C")
xs = np.linspace(ln_sinh.min(), ln_sinh.max(), 50)
ax.plot(xs, slope_Z * xs + lnA, "-", color="#333333", lw=1.8,
        label=rf"Fit: $\ln Z$ = {slope_Z:.2f}$x$ + {lnA:.2f}")
format_axes(ax, LBL_LN_SINH, r"$\ln Z$  ($Z$ in s$^{-1}$)",
            r"(e) Determination of $\ln A$")
ax.legend(frameon=False)
save(fig, "5_lnA_lnZ_vs_ln_sinh.png")

# ---------------------------------------------------------------------------
# Step 7: model validation
#   sigma = (1/alpha) asinh[(Z/A)^(1/n)]
# ---------------------------------------------------------------------------
Z = np.exp(ln_Z)
predicted = (1.0 / alpha) * np.arcsinh((Z / A) ** (1.0 / n))

measured = stress.ravel()
pred = predicted.ravel()
r_squared = 1 - (np.sum((measured - pred) ** 2)
                 / np.sum((measured - measured.mean()) ** 2))
aare = np.mean(np.abs((measured - pred) / measured)) * 100  # per cent

fig, ax = plt.subplots(figsize=(5.8, 5.6))
lims = [min(measured.min(), pred.min()) * 0.9,
        max(measured.max(), pred.max()) * 1.05]
ax.plot(lims, lims, "--", color="#888888", lw=1.4, label="Line of equality")
for i, T in enumerate(temps_C):
    ax.plot(stress[i], predicted[i], MARKERS[i], color=PALETTE[i], ms=9,
            label=f"{T} °C")
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_aspect("equal")
format_axes(ax, "Measured peak stress (MPa)", "Predicted peak stress (MPa)",
            "(f) Predicted vs measured peak stress")
ax.text(0.97, 0.05, f"$R^2$ = {r_squared:.3f}\nAARE = {aare:.1f} %",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=12,
        bbox=dict(boxstyle="round", fc="white", ec="#888888"))
ax.legend(frameon=False, loc="upper left")
save(fig, "6_predicted_vs_measured.png")

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
print("Slopes at each temperature (250, 350, 450 °C):")
print("  n'  :", np.round(n1_each, 3))
print("  beta:", np.round(beta_each, 5))
print("  n   :", np.round(n_each, 3))
print("Slopes of ln[sinh(alpha*sigma)] vs 1000/T at each strain rate "
      "(0.1, 1, 5 s^-1):", np.round(K_each, 3), "\n")

results = pd.DataFrame(
    {"value": [n1, beta, alpha, n, Q_kJ, lnA, r_squared, aare]},
    index=["n' (power-law exponent)", "beta (MPa^-1)", "alpha (MPa^-1)",
           "n (stress exponent)", "Q (kJ/mol)", "ln A (A in s^-1)",
           "R^2", "AARE (%)"])
results.index.name = "constant"
pd.set_option("display.float_format", lambda v: f"{v:.5g}")
print("Fitted constants")
print(results, "\n")

table = data.copy()
table["ln_Z"] = ln_Z.ravel()
table["predicted_MPa"] = pred
table["error_percent"] = (pred - measured) / measured * 100
print(table.round(2).to_string(index=False), "\n")
print(f"R^2 = {r_squared:.4f}   AARE = {aare:.2f} %")
print(f"Figures saved to: {FIG_DIR}")

results.to_csv(os.path.join(HERE, "results_constants.csv"))
table.to_csv(os.path.join(HERE, "results_predictions.csv"), index=False)
