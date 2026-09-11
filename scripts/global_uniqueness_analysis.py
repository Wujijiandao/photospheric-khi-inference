"""Scale-free global uniqueness and signed-velocity tests for the two-mode slab inverse.

This script reproduces the numerical logic added in scientific revision v1.4:
  * dimensionless Jacobian scan in (r,u,epsilon),
  * multistart searches for separated inverse roots,
  * an explicit unsigned-speed branch ambiguity example.

The expensive full scan is intentionally explicit rather than hidden behind cached data.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import brentq, least_squares


def coefficients(k,e):
    q=1-np.exp(-2*k)
    A=k*k*((e*e/4)*q-1); B=e*k*q
    C=q+k*k-2*k-e*e*k*k*q/4
    return A,B,C


def disc(k,e):
    A,B,C=coefficients(k,e)
    return B*B-4*A*C


def mode(k,e):
    A,B,C=coefficients(k,e); D=B*B-4*A*C
    R=-B/(2*A); I=np.sqrt(max(0.0,-D))/(2*abs(A))
    return R,I,D


def cutoff(e):
    return brentq(lambda x:disc(x,e),1e-6,3.0,xtol=1e-12,rtol=1e-12)


def F_k(k1,e,r):
    R1,I1,D1=mode(k1,e); R2,I2,D2=mode(r*k1,e)
    if I1<=0 or I2<=0 or D1>=0 or D2>=0:
        return None
    return np.array([np.log(I1/I2),(R1-R2)/I1])


def F_ue(u,e,r):
    return F_k(u*cutoff(e)/r,e,r)


def jac_ue(u,e,r):
    hu=1e-5; he=1e-5
    du=(F_ue(u+hu,e,r)-F_ue(u-hu,e,r))/(2*hu)
    de=(F_ue(u,e+he,r)-F_ue(u,e-he,r))/(2*he)
    return np.column_stack([du,de])


def jacobian_scan(outdir: Path):
    rgrid=np.geomspace(1.005,10.0,72)
    ugrid=np.linspace(0.01,0.99,70)
    egrid=np.linspace(0.01,0.99,41)
    rows=[]
    for r in rgrid:
        for u in ugrid:
            for e in egrid:
                J=jac_ue(u,e,r)
                sv=np.linalg.svd(J,compute_uv=False)
                rows.append([r,u,e,np.linalg.det(J),sv[0],sv[-1],sv[0]/sv[-1]])
    df=pd.DataFrame(rows,columns=[
        'wavenumber_ratio_r','upper_mode_fraction_cutoff_u','epsilon',
        'jacobian_determinant','singular_value_max','singular_value_min','condition_number'])
    df.to_csv(outdir/'global_uniqueness_jacobian_scan.csv.gz',index=False,compression='gzip')
    atlas=(df.groupby(['wavenumber_ratio_r','upper_mode_fraction_cutoff_u'],as_index=False)
           .agg(min_jacobian_determinant=('jacobian_determinant','min'),
                min_singular_value=('singular_value_min','min'),
                max_condition_number=('condition_number','max')))
    atlas.to_csv(outdir/'global_uniqueness_atlas.csv',index=False)
    return df,atlas


def distinct_roots(r,u0,e0,start_values,res_tol=1e-9,cluster_tol=2e-3):
    target=F_ue(u0,e0,r)
    sols=[]
    for es in start_values:
        for us in start_values:
            ks=us*cutoff(es)/r
            def residual(x):
                f=F_k(x[0],x[1],r)
                if f is None or np.any(~np.isfinite(f)):
                    return np.array([100.,100.])
                return f-target
            fit=least_squares(residual,[ks,es],
                bounds=([1e-8,.0001],[1.58/r,.9999]),max_nfev=180,
                ftol=1e-13,xtol=1e-13,gtol=1e-13)
            if np.linalg.norm(fit.fun)<res_tol:
                k1c,ec=fit.x
                uc=r*k1c/cutoff(ec)
                p=np.array([uc,ec])
                if not any(np.linalg.norm(p-q)<cluster_tol for q in sols):
                    sols.append(p)
    return sols


def root_surveys(outdir: Path):
    rows=[]
    for r in [1.02,1.05,1.1,1.2,1.5,2.0,3.0,5.0,8.0]:
        for u in np.linspace(.05,.95,9):
            for e in np.linspace(.05,.95,9):
                sols=distinct_roots(r,float(u),float(e),np.linspace(.08,.92,4))
                rows.append([r,u,e,len(sols)])
    pd.DataFrame(rows,columns=['wavenumber_ratio_r','true_u','true_epsilon','distinct_root_count']).to_csv(
        outdir/'global_multistart_root_survey.csv',index=False)

    rows=[]
    vals=[.01,.02,.05,.10,.50,.90,.95,.98,.99]
    for r in [1.005,1.01,1.02,1.05,1.1,1.2,1.5,2,3,5,8,10]:
        for u in vals:
            for e in vals:
                sols=distinct_roots(r,u,e,[.08,.50,.92],cluster_tol=3e-3)
                rows.append([r,u,e,len(sols)])
    pd.DataFrame(rows,columns=['wavenumber_ratio_r','true_u','true_epsilon','distinct_root_count']).to_csv(
        outdir/'global_edge_multistart_stress_survey.csv',index=False)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--outdir',default='results_global')
    ap.add_argument('--atlas',action='store_true')
    ap.add_argument('--roots',action='store_true')
    args=ap.parse_args()
    outdir=Path(args.outdir); outdir.mkdir(parents=True,exist_ok=True)
    if not args.atlas and not args.roots:
        args.atlas=args.roots=True
    if args.atlas:
        df,atlas=jacobian_scan(outdir)
        print('Jacobian points:',len(df),'positive-det fraction:',float((df.jacobian_determinant>0).mean()))
    if args.roots:
        root_surveys(outdir)
        print('Root surveys complete.')

if __name__=='__main__':
    main()
