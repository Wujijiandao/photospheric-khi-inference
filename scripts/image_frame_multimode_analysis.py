#!/usr/bin/env python3
"""Reproduce the image-frame two-mode identifiability and recovery analysis.

Outputs:
  image_frame_pair_conditioning.csv
  image_frame_multimode_recovery.csv
  image_frame_global_scan.csv
  fig08_image_frame_identifiability.pdf
  fig09_image_frame_recovery.pdf
  IMAGE_FRAME_MULTIMODE_SUMMARY.json
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parent


def coefficients(kappa, epsilon):
    q = 1.0 - np.exp(-2.0*kappa)
    A = kappa**2*((epsilon**2/4.0)*q - 1.0)
    B = epsilon*kappa*q
    C = q + kappa**2 - 2.0*kappa - epsilon**2*kappa**2*q/4.0
    return A, B, C


def mode(kappa, epsilon):
    A, B, C = coefficients(kappa, epsilon)
    D = B*B - 4*A*C
    R = -B/(2*A)
    I = np.sqrt(max(0.0, -D))/(2*abs(A))
    return R, I, D


def cutoff(epsilon):
    xs = np.linspace(1e-4, 3.0, 3000)
    ds = np.array([mode(x, epsilon)[2] for x in xs])
    ix = np.where((ds[:-1] < 0) & (ds[1:] >= 0))[0]
    if len(ix) == 0:
        raise RuntimeError("No short-wave cutoff found")
    j = ix[-1]
    return brentq(lambda x: mode(x, epsilon)[2], xs[j], xs[j+1])


def invariants(d, epsilon, k1, k2):
    R1, I1, D1 = mode(k1*d, epsilon)
    R2, I2, D2 = mode(k2*d, epsilon)
    if I1 <= 0 or I2 <= 0 or D1 >= 0 or D2 >= 0:
        return np.array([np.nan, np.nan])
    return np.array([np.log(I1/I2), (R1-R2)/I1])


def jacobian(d, epsilon, k1, k2):
    x = np.log(d)
    hx, he = 1e-5, 1e-5
    f = lambda xx, ee: invariants(np.exp(xx), ee, k1, k2)
    dx = (f(x+hx, epsilon)-f(x-hx, epsilon))/(2*hx)
    de = (f(x, epsilon+he)-f(x, epsilon-he))/(2*he)
    return np.column_stack([dx, de])


def forward(d, U0, epsilon, vc, fractions=(0.30, 0.85)):
    kc = cutoff(epsilon)
    kappas = np.array(fractions)*kc
    k = kappas/d
    lam = 2*np.pi/k
    gamma, vapp = [], []
    for ki, kap in zip(k, kappas):
        R, I, _ = mode(kap, epsilon)
        gamma.append(ki*U0*I)
        vapp.append(vc+U0*R)
    return lam, np.array(gamma), np.array(vapp), kappas


def invert_two(lam, gamma, vapp, x0=(10.0, 0.6)):
    lam = np.asarray(lam, float)
    gamma = np.asarray(gamma, float)
    vapp = np.asarray(vapp, float)
    k = 2*np.pi/lam
    g = gamma/k
    if np.any(g <= 0):
        return None
    obs = np.array([np.log(g[0]/g[1]), (vapp[0]-vapp[1])/g[0]])

    def residual(x):
        th = invariants(x[0], x[1], k[0], k[1])
        if np.any(~np.isfinite(th)):
            return np.array([100.0, 100.0])
        return th-obs

    fit = least_squares(residual, x0, bounds=([2.0, 0.01], [30.0, 0.99]),
                        max_nfev=150, xtol=1e-10, ftol=1e-10, gtol=1e-10)
    if np.linalg.norm(fit.fun) > 1e-3:
        return None
    d, epsilon = fit.x
    RI = [mode(ki*d, epsilon)[:2] for ki in k]
    if any(I <= 1e-10 for R, I in RI):
        return None
    U0 = float(np.mean([gi/(ki*I) for gi, ki, (R, I) in zip(gamma, k, RI)]))
    vc = float(np.mean([vi-U0*R for vi, (R, I) in zip(vapp, RI)]))
    return d, epsilon, U0, vc, float(np.linalg.norm(fit.fun))


def main():
    dtrue, Utrue, etrue, vctrue = 10.0, 1.5, 0.6, 1.2
    fractions = (0.30, 0.85)
    lam0, gam0, vapp0, _ = forward(dtrue, Utrue, etrue, vctrue, fractions)
    kc0 = cutoff(etrue)
    k0 = 2*np.pi/lam0
    J = jacobian(dtrue, etrue, k0[0], k0[1])
    sv = np.linalg.svd(J, compute_uv=False)

    sgrid = np.linspace(0.10, 0.95, 86)
    rows, svmap = [], np.full((len(sgrid), len(sgrid)), np.nan)
    for i, s1 in enumerate(sgrid):
        for j, s2 in enumerate(sgrid):
            if s2 <= s1+0.02:
                continue
            JJ = jacobian(dtrue, etrue, s1*kc0/dtrue, s2*kc0/dtrue)
            if np.all(np.isfinite(JJ)):
                ss = np.linalg.svd(JJ, compute_uv=False)
                svmap[i, j] = ss[-1]
                rows.append([s1, s2, np.linalg.det(JJ), ss[0], ss[-1], ss[0]/ss[-1]])
    pd.DataFrame(rows, columns=["fraction1_cutoff", "fraction2_cutoff", "determinant",
                                "singular_value_max", "singular_value_min", "condition_number"]
                ).to_csv(ROOT/"image_frame_pair_conditioning.csv", index=False)

    plt.figure(figsize=(6.6, 5.1))
    im = plt.imshow(np.ma.masked_invalid(np.log10(svmap)), origin="lower",
                    extent=[sgrid[0], sgrid[-1], sgrid[0], sgrid[-1]], aspect="auto")
    plt.colorbar(im, label=r"$\log_{10}\,\sigma_{\min}(J_{\rm img})$")
    plt.scatter([fractions[1]], [fractions[0]], marker="x", s=70, label="recovery design")
    plt.xlabel(r"Second mode: $\kappa_2/\kappa_c$")
    plt.ylabel(r"First mode: $\kappa_1/\kappa_c$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(ROOT/"fig08_image_frame_identifiability.pdf")
    plt.close()

    rng = np.random.default_rng(20260911)
    nmc = 10000
    mcrows = []
    for sigv in [0.01, 0.03, 0.05]:
        vals, rejected = [], 0
        for _ in range(nmc):
            lam = lam0*(1+rng.normal(0, 0.02, 2))
            gam = gam0*(1+rng.normal(0, 0.05, 2))
            vapp = vapp0*(1+rng.normal(0, sigv, 2))
            if np.any(lam <= 0) or np.any(gam <= 0):
                rejected += 1
                continue
            out = invert_two(lam, gam, vapp)
            if out is None:
                rejected += 1
                continue
            d, e, U0, vc, _ = out
            if not (2.01 < d < 29.99 and 0.011 < e < 0.989):
                rejected += 1
                continue
            vals.append([d, e, U0, vc])
        a = np.asarray(vals)
        row = {"sigma_lambda_rel": 0.02, "sigma_gamma_rel": 0.05,
               "sigma_vapp_rel": sigv, "accepted": len(a), "rejected": rejected}
        for col, name in enumerate(["d_km", "epsilon", "U0_km_s", "vc_km_s"]):
            q16, q50, q84 = np.quantile(a[:, col], [0.16, 0.5, 0.84])
            row[f"{name}_q16"] = q16
            row[f"{name}_median"] = q50
            row[f"{name}_q84"] = q84
        mcrows.append(row)
    mc = pd.DataFrame(mcrows)
    mc.to_csv(ROOT/"image_frame_multimode_recovery.csv", index=False)

    plt.figure(figsize=(6.7, 4.7))
    x = np.arange(3)
    for name, trueval, label, offset in [
        ("d_km", dtrue, r"$d/d_{\rm true}$", -0.12),
        ("epsilon", etrue, r"$\epsilon/\epsilon_{\rm true}$", -0.04),
        ("U0_km_s", Utrue, r"$U_0/U_{0,\rm true}$", 0.04),
        ("vc_km_s", vctrue, r"$v_c/v_{c,\rm true}$", 0.12),
    ]:
        med = mc[f"{name}_median"].to_numpy()/trueval
        lo = (mc[f"{name}_median"]-mc[f"{name}_q16"]).to_numpy()/trueval
        hi = (mc[f"{name}_q84"]-mc[f"{name}_median"]).to_numpy()/trueval
        plt.errorbar(x+offset, med, yerr=np.vstack([lo, hi]), fmt="o", capsize=4, label=label)
    plt.axhline(1.0, linewidth=1)
    plt.xticks(x, ["1%", "3%", "5%"])
    plt.xlabel(r"Relative uncertainty of each $v_{\rm app}$")
    plt.ylabel("Recovered parameter / true parameter")
    plt.legend(fontsize=8, ncol=2)
    plt.tight_layout()
    plt.savefig(ROOT/"fig09_image_frame_recovery.pdf")
    plt.close()

    obs = invariants(dtrue, etrue, k0[0], k0[1])
    scan = []
    for d in np.linspace(2, 30, 141):
        for e in np.linspace(0.01, 0.99, 99):
            th = invariants(d, e, k0[0], k0[1])
            if np.all(np.isfinite(th)):
                scan.append([d, e, float(np.linalg.norm(th-obs))])
    pd.DataFrame(scan, columns=["d_km", "epsilon", "invariant_residual"]).to_csv(
        ROOT/"image_frame_global_scan.csv", index=False)

    summary = {
        "truth": {"d_km": dtrue, "U0_km_s": Utrue, "epsilon": etrue, "vc_km_s": vctrue},
        "mode_fractions_of_cutoff": list(fractions),
        "wavelength_km": lam0.tolist(), "growth_s_inv": gam0.tolist(),
        "vapp_km_s": vapp0.tolist(), "J_img": J.tolist(),
        "det_J_img": float(np.linalg.det(J)), "singular_values_J_img": sv.tolist(),
    }
    (ROOT/"IMAGE_FRAME_MULTIMODE_SUMMARY.json").write_text(json.dumps(summary, indent=2))
    print(mc.to_string(index=False))


if __name__ == "__main__":
    main()
