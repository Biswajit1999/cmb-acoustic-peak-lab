"""Auditable diagnostics for the *phenomenological* CMB toy spectrum.

These tests quantify mismatches and parameter degeneracies. They must not be
used to infer Planck cosmological posteriors; use CLASS/CAMB plus the official
likelihood and covariance for that purpose.
"""
from __future__ import annotations
import json
import math
from pathlib import Path
from typing import Sequence

DEFAULT = {"omegaM": 0.315, "h": 0.674, "ns": 0.965, "lensing": 0.18}
PEAKS = [(220,5740,80), (540,2550,110), (800,2500,130), (1120,1200,170), (1450,650,220)]


def toy_dl(ell: float, params: dict | None = None) -> float:
    """Direct Python reproduction of the JavaScript worker's toy TT D_ell."""
    p = DEFAULT if params is None else params
    if not math.isfinite(ell) or ell <= 0:
        raise ValueError("multipole must be positive and finite")
    if not all(math.isfinite(p[k]) for k in DEFAULT):
        raise ValueError("nonfinite model parameters")
    if p['omegaM'] <= 0 or p['h'] <= 0:
        raise ValueError("omegaM and h must be positive")
    shift=(p['omegaM']/.315)**.18*(.674/p['h'])**.25
    base=850*(ell/80)**.2*math.exp(-ell/1600)
    for center,amp,width in PEAKS:
        base+=amp*math.exp(-.5*((ell-center*shift)/width)**2)
    return base*(ell/220)**(p['ns']-.965)*math.exp(-p['lensing']*ell/4800)


def load_planck_points(path: Path) -> list[tuple[float,float,float]]:
    """Read existing repository's Planck binned TT data without modifying it."""
    d=json.loads(path.read_text(encoding='utf-8'))
    points=[(float(v['x']),float(v['y']),float(v['y_err'])) for v in d['points']]
    if len(points)<3 or any(not all(math.isfinite(x) for x in p) or p[0]<=0 or p[2]<=0 for p in points):
        raise ValueError("invalid Planck points")
    if any(a[0]>=b[0] for a,b in zip(points,points[1:])):
        raise ValueError("multipoles must be increasing")
    if d.get('n_points') is not None and d['n_points'] != len(points):
        raise ValueError("reference point count mismatch")
    return points


def fit_amplitude(points: Sequence[tuple[float,float,float]], params: dict | None = None,
                  model_values: Sequence[float] | None = None) -> dict:
    """Analytic diagonal-weighted amplitude fit with residual chi-square.

    Treats published symmetrized error bars as independent for *diagnostics*
    only. Official Planck likelihood is not diagonal in these bins.
    """
    if len(points)<2:
        raise ValueError("need at least two bins")
    m = list(model_values) if model_values is not None else [toy_dl(x,params) for x,_,_ in points]
    if len(m)!=len(points) or any(not math.isfinite(v) for v in m):
        raise ValueError("invalid model values")
    if any(s<=0 or not math.isfinite(s) for _,_,s in points):
        raise ValueError("invalid uncertainties")
    den=sum((v/s)**2 for v,(_,_,s) in zip(m,points))
    if den<=0:
        raise ValueError("zero model information")
    amp=sum(v*y/s**2 for v,(_,y,s) in zip(m,points))/den
    chi2=sum(((y-amp*v)/s)**2 for v,(_,y,s) in zip(m,points))
    return {"amplitude":amp,"amplitude_sigma_diagonal":1/math.sqrt(den),
            "chi2_diagonal":chi2,"dof":len(points)-1,
            "reduced_chi2_diagonal":chi2/(len(points)-1)}


def _simpson(f, left: float, right: float, steps: int=16) -> float:
    if steps<2 or steps%2:
        raise ValueError("Simpson steps must be positive and even")
    h=(right-left)/steps
    return h/3*(f(left)+f(right)+4*sum(f(left+i*h) for i in range(1,steps,2))+
                2*sum(f(left+i*h) for i in range(2,steps,2)))


def top_hat_bin_diagnostics(points: Sequence[tuple[float,float,float]], params: dict | None=None,
                            steps: int=16) -> list[dict]:
    """Compare centre evaluation with approximate top-hat bin integration.

    Midpoint edges are inferred from bin centres, NOT official Planck windows.
    """
    if len(points)<3 or any(points[i][0]>=points[i+1][0] for i in range(len(points)-1)):
        raise ValueError("need at least three sorted multipole bins")
    x=[p[0] for p in points]
    edges=[x[0]-(x[1]-x[0])/2]+[(a+b)/2 for a,b in zip(x,x[1:])]+[x[-1]+(x[-1]-x[-2])/2]
    if edges[0]<=0:
        raise ValueError("first bin crosses nonphysical ell<=0")
    result=[]
    for i,ell in enumerate(x):
        center=toy_dl(ell,params)
        avg=_simpson(lambda u:toy_dl(u,params),edges[i],edges[i+1],steps)/(edges[i+1]-edges[i])
        result.append({"ell":ell,"left":edges[i],"right":edges[i+1],
                       "center_dl":center,"top_hat_dl":avg,
                       "fractional_difference":(avg-center)/center})
    return result


def sensitivity_cosines(points: Sequence[tuple[float,float,float]], params: dict | None=None,
                        fractional_step: float=1e-4) -> dict:
    """Weighted local finite-difference derivative cosines; +/-1 is degeneracy.

    This is *not* a posterior correlation or a Planck Fisher forecast.
    """
    if not (0<fractional_step<.05):
        raise ValueError("invalid step")
    if len(points)<3:
        raise ValueError("need three or more bins")
    p=dict(DEFAULT if params is None else params)
    keys=list(DEFAULT)
    deriv={}
    for key in keys:
        h=fractional_step*max(abs(p[key]),0.1)
        up=dict(p); down=dict(p)
        up[key]+=h;down[key]-=h
        deriv[key]=[(toy_dl(x,up)-toy_dl(x,down))/(2*h*s)
                    for x,_,s in points]
    norms={k:math.sqrt(sum(v*v for v in deriv[k])) for k in keys}
    if any(v==0 for v in norms.values()):
        raise ValueError("unidentifiable derivative")
    cosine={}
    for i,a in enumerate(keys):
        for b in keys[i+1:]:
            cosine[f"{a}:{b}"]=sum(u*v for u,v in zip(deriv[a],deriv[b]))/(norms[a]*norms[b])
    return {"weighted_derivative_norms":norms,"cosines":cosine,
            "note":"Model omegaM and h enter only through one peak-shift combination; they cannot be separated by this toy spectrum."}


def residual_by_band(points: Sequence[tuple[float,float,float]], amplitude: float,
                     params: dict | None=None) -> list[dict]:
    """Summarize residuals by broad, disjoint multipole intervals."""
    bands=[(0,400),(400,1000),(1000,1600),(1600,2600)]
    out=[]
    for lo,hi in bands:
        part=[(x,y,s) for x,y,s in points if lo<=x<hi]
        if part:
            resid=[(y-amplitude*toy_dl(x,params))/s for x,y,s in part]
            out.append({"ell_range":[lo,hi],"n_bins":len(part),
                        "rms_normalized_residual":math.sqrt(sum(v*v for v in resid)/len(resid)),
                        "mean_normalized_residual":sum(resid)/len(resid)})
    return out