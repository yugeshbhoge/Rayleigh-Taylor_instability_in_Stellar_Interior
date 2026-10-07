#!/usr/bin/env python3
"""
2D finite-volume ideal MHD for magneto-RTI (mu0 = 1, gravity along -z)
Companion code to the PH-819 project report "Rayleigh-Taylor Instability in Stellar Interiors".

Physics
-------
Compressible ideal MHD in the (x, z) plane with uniform gravity g = -g z_hat.
Units: mu_0 = 1, magnetic pressure = B^2/2 and v_A = B/sqrt(rho).
    d_t rho   + div(rho v)                                   = 0
    d_t(rho v) + div(rho v v + (p + B^2/2) I - B B)          = rho g
    d_t E     + div((E + p + B^2/2) v - B (v.B))             = rho v.g
    d_t B     + div(v B - B v) + grad(psi)                   = 0      (GLM)
    d_t psi   + c_h^2 div B                                  = -(c_h/L) psi ...

    E = p/(gamma-1) + rho v^2/2 + B^2/2

Numerics
--------
* Finite-volume, cell-centred, 2 ghost layers
* MUSCL (minmod) reconstruction of primitive variables
* Rusanov / local Lax-Friedrichs flux using the fast-magnetosonic speed
* SSP-RK2 (Heun) time stepping, gravity as a source term
* GLM hyperbolic-parabolic divergence cleaning for div B
* x periodic, z: reflecting, conducting walls with hydrostatic ghost cells

Grid layout
-----------------------------------------------
* Compelte MHD: 256 x 256 grid.
* Smooth tanh layer INTERFACE_CELLS = 16 cells thick.

Linear theory used for comparison (report Eqs. 3.15-3.22, k || B0):
    gamma^2 = At g k - 2 k^2 B0^2 / (rho_t + rho_b)    (mu_0 = 1)
    k_c     = g (rho_t - rho_b) / (2 B0^2),   k_m = k_c / 2

How to run
-----
    python3 rti_mhd_2d.py                       # B0 = 0.2, 256x256, t_end = 3
    python3 rti_mhd_2d.py --b0 0                # pure hydrodynamic RTI
    python3 rti_mhd_2d.py --pr exponential  # B_x(z) = B0 exp((z-z_i)/L_B)
    python3 rti_mhd_2d.py --nx 128 --nz 128 --ic 8 --te 2.5 --md 1,2,3,4 --am 0.01   # quick test
"""

import argparse as ap
import os
import time as tm
from dataclasses import dataclass as dc, field as fl

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

NG = 2
IR, IX, IZ, IE, BX, BZ, PS = range(7)
SW = [0, 2, 1, 3, 5, 4, 6]
RF = 1e-6
PF = 1e-6


@dc
class Cf:
    nx: int = 256
    nz: int = 256
    lx: float = 1.0
    lz: float = 1.0
    gm: float = 5.0 / 3.0
    g: float = 1.0
    rt: float = 2.0
    rb: float = 1.0
    p0: float = 5.0
    b0: float = 0.2
    bz: float = 0.0
    pr: str = "uniform"
    lb: float = 1.0
    ic: int = 16
    cg: int = 16
    md: tuple = (1, 2, 4, 8)
    am: float = 0.005
    sg: float = 0.05
    sd: int = 1
    cfl: float = 0.4
    te: float = 3.0
    pa: float = 0.4
    ns: int = 6
    od: str = "rti_output"
    tg: str = ""
    dx: float = fl(init=False)
    dz: float = fl(init=False)
    zi: float = fl(init=False)

    def __post_init__(self):
        self.dx = self.lx / self.nx
        self.dz = self.lz / self.nz
        self.zi = 0.5 * self.lz


def c2p(U, gm):
    r = np.maximum(U[IR], RF)
    vx = U[IX] / r
    vz = U[IZ] / r
    bx, bz = U[BX], U[BZ]
    p = (gm - 1.0) * (U[IE] - 0.5 * r * (vx**2 + vz**2) - 0.5 * (bx**2 + bz**2))
    return np.array([r, vx, vz, np.maximum(p, PF), bx, bz, U[PS]])


def p2c(W, gm):
    r, vx, vz, p, bx, bz, ps = W
    e = p / (gm - 1.0) + 0.5 * r * (vx**2 + vz**2) + 0.5 * (bx**2 + bz**2)
    return np.array([r, r * vx, r * vz, e, bx, bz, ps])


def cf(r, p, b2, bn2, gm):
    s2 = gm * p / r
    a2 = b2 / r
    n2 = bn2 / r
    d = np.maximum((s2 + a2) ** 2 - 4.0 * s2 * n2, 0.0)
    return np.sqrt(0.5 * (s2 + a2 + np.sqrt(d)))


def mm(a, b):
    return 0.5 * (np.sign(a) + np.sign(b)) * np.minimum(np.abs(a), np.abs(b))


def sl(W, ax):
    a = W - np.roll(W, 1, axis=ax)
    b = np.roll(W, -1, axis=ax) - W
    return mm(a, b)


def fs(W, gm, ch):
    r = np.maximum(W[0], RF)
    vn, vt = W[1], W[2]
    p = np.maximum(W[3], PF)
    bn, bt, ps = W[4], W[5], W[6]
    b2 = bn**2 + bt**2
    pt = p + 0.5 * b2
    e = p / (gm - 1.0) + 0.5 * r * (vn**2 + vt**2) + 0.5 * b2
    vb = vn * bn + vt * bt
    F = np.empty((7,) + r.shape)
    F[0] = r * vn
    F[1] = r * vn * vn + pt - bn * bn
    F[2] = r * vn * vt - bn * bt
    F[3] = (e + pt) * vn - bn * vb
    F[4] = ps
    F[5] = vn * bt - vt * bn
    F[6] = ch**2 * bn
    U = np.array([r, r * vn, r * vt, e, bn, bt, ps])
    a = np.abs(vn) + cf(r, p, b2, bn**2, gm)
    return F, U, a


def rs(L, R, gm, ch):
    FL, UL, aL = fs(L, gm, ch)
    FR, UR, aR = fs(R, gm, ch)
    a = np.maximum(aL, aR)
    F = 0.5 * (FL + FR) - 0.5 * a * (UR - UL)
    for i in (BX, PS):
        F[i] = 0.5 * (FL[i] + FR[i]) - 0.5 * ch * (UR[i] - UL[i])
    return F


class Sv:
    def __init__(self, c):
        self.c = c
        self.U = np.zeros((7, c.nx + 2 * NG, c.nz + 2 * NG))
        self.x = (np.arange(c.nx) + 0.5) * c.dx
        self.z = (np.arange(c.nz) + 0.5) * c.dz
        self.ch = 1.0
        self.t = 0.0
        self.ini()

    def bp(self, z):
        c = self.c
        if c.pr == "uniform":
            return c.b0 * np.ones_like(z)
        if c.pr == "exponential":
            return c.b0 * np.exp((z - c.zi) / c.lb)
        raise ValueError(f"unknown profile {c.pr}")

    def ini(self):
        c = self.c
        X, Z = np.meshgrid(self.x, self.z, indexing="ij")
        dl = c.ic * c.dz / 4.0
        s = (Z - c.zi) / dl
        dr = c.rt - c.rb
        r = c.rb + 0.5 * dr * (1.0 + np.tanh(s))
        lc = np.logaddexp(s, -s) - np.log(2.0)
        pt = c.p0 - c.g * ((c.rb + 0.5 * dr) * (Z - c.zi) + 0.5 * dr * dl * lc)
        bx = self.bp(Z)
        bz = c.bz * np.ones_like(Z)
        p = pt - 0.5 * (bx**2 + bz**2)
        if p.min() <= 0:
            raise RuntimeError("negative gas pressure: lower b0 or raise p0")
        rg = np.random.default_rng(c.sd)
        vz = np.zeros_like(Z)
        for m in c.md:
            vz += c.am * np.cos(2 * np.pi * m * X / c.lx + rg.uniform(0, 2 * np.pi))
        vz *= np.exp(-((Z - c.zi) / c.sg) ** 2)
        W = np.array([r, np.zeros_like(r), vz, p, bx, bz, np.zeros_like(r)])
        self.U[:, NG:-NG, NG:-NG] = p2c(W, c.gm)
        self.bc(self.U)

    def bc(self, U):
        c = self.c
        U[:, :NG, :] = U[:, -2 * NG:-NG, :]
        U[:, -NG:, :] = U[:, NG:2 * NG, :]
        wb = c2p(U[:, :, NG], c.gm)
        wt = c2p(U[:, :, NG + c.nz - 1], c.gm)
        for k in range(1, NG + 1):
            w = wb.copy()
            w[2] = -w[2]
            w[5] = -w[5]
            w[3] = w[3] + w[0] * c.g * k * c.dz
            U[:, :, NG - k] = p2c(w, c.gm)
            w = wt.copy()
            w[2] = -w[2]
            w[5] = -w[5]
            w[3] = np.maximum(w[3] - w[0] * c.g * k * c.dz, PF)
            U[:, :, NG + c.nz - 1 + k] = p2c(w, c.gm)

    def cdt(self):
        c = self.c
        r, vx, vz, p, bx, bz = c2p(self.U[:, NG:-NG, NG:-NG], c.gm)[:6]
        b2 = bx**2 + bz**2
        sx = np.abs(vx) + cf(r, p, b2, bx**2, c.gm)
        sz = np.abs(vz) + cf(r, p, b2, bz**2, c.gm)
        self.ch = max(sx.max(), sz.max())
        return c.cfl / (sx.max() / c.dx + sz.max() / c.dz)

    def swp(self, W, ax, sp):
        c = self.c
        if sp:
            W = W[SW]
        s = sl(W, ax)
        lo = [slice(None)] * 3
        hi = [slice(None)] * 3
        lo[ax] = slice(0, -1)
        hi[ax] = slice(1, None)
        L = (W + 0.5 * s)[tuple(lo)]
        R = (W - 0.5 * s)[tuple(hi)]
        F = rs(L, R, c.gm, self.ch)
        return F[SW] if sp else F

    def rh(self, U):
        c = self.c
        self.bc(U)
        W = c2p(U, c.gm)
        Fx = self.swp(W[:, :, NG:-NG], 1, False)
        Fz = self.swp(W[:, NG:-NG, :], 2, True)
        L = -(Fx[:, NG:NG + c.nx, :] - Fx[:, NG - 1:NG + c.nx - 1, :]) / c.dx
        L -= (Fz[:, :, NG:NG + c.nz] - Fz[:, :, NG - 1:NG + c.nz - 1]) / c.dz
        L[IZ] -= U[IR, NG:-NG, NG:-NG] * c.g
        L[IE] -= U[IZ, NG:-NG, NG:-NG] * c.g
        return L

    def stp(self, dt):
        c = self.c
        U0 = self.U[:, NG:-NG, NG:-NG].copy()
        L0 = self.rh(self.U)
        self.U[:, NG:-NG, NG:-NG] = U0 + dt * L0
        L1 = self.rh(self.U)
        self.U[:, NG:-NG, NG:-NG] = 0.5 * U0 + 0.5 * (self.U[:, NG:-NG, NG:-NG] + dt * L1)
        self.U[PS] *= np.exp(-c.pa * self.ch * dt / min(c.dx, c.dz))
        self.t += dt

    def pm(self):
        return c2p(self.U[:, NG:-NG, NG:-NG], self.c.gm)

    def mx(self):
        c = self.c
        r = self.U[IR, NG:-NG, NG:-NG]
        return np.clip((r - c.rb) / (c.rt - c.rb), 0.0, 1.0)


def dg(s):
    c = s.c
    f = s.mx()
    fb = f.mean(axis=0)
    hm = np.sum(4.0 * fb * (1.0 - fb)) * c.dz
    et = c.lz - f.sum(axis=1) * c.dz - c.zi
    eh = 2.0 * np.abs(np.fft.rfft(et)) / c.nx
    W = s.pm()
    bx, bz = W[4], W[5]
    dv = np.gradient(bx, c.dx, axis=0) + np.gradient(bz, c.dz, axis=1)
    bm = np.sqrt(bx**2 + bz**2).max()
    db = np.abs(dv).max() * min(c.dx, c.dz) / max(bm, 1e-12)
    return dict(t=s.t, hm=hm, am=np.array([eh[m] for m in c.md]), db=db,
                vm=np.abs(W[2]).max())


def bk(a, n):
    nx, nz = a.shape
    return a.reshape(n, nx // n, n, nz // n).mean(axis=(1, 3))


def g2t(k, c):
    return ((c.rt - c.rb) / (c.rt + c.rb) * c.g * k
            - 2.0 * (k * c.b0) ** 2 / (c.rt + c.rb))


def fg(t, a, k, cap=0.1):
    lo, hi = 0.01 / k, cap / k
    i = np.where(a >= lo)[0]
    if len(i) == 0:
        return np.nan
    i0 = i[0]
    o = np.where(a[i0:] > hi)[0]
    i1 = i0 + (o[0] if len(o) else len(a[i0:]))
    if i1 - i0 < 5 or a[i1 - 1] < 2.0 * a[i0]:
        return np.nan
    return np.polyfit(t[i0:i1], np.log(a[i0:i1]), 1)[0]


def pls(sn, c, fn):
    n = len(sn)
    nc = 3
    nr = int(np.ceil(n / nc))
    fig, ax = plt.subplots(nr, nc, figsize=(4.3 * nc, 4.5 * nr), squeeze=False)
    xs = (np.arange(c.nx) + 0.5) * c.dx
    zs = (np.arange(c.nz) + 0.5) * c.dz
    st = max(c.nx // 64, 1)
    for a, (t, r, bx, bz) in zip(ax.ravel(), sn):
        im = a.imshow(r.T, origin="lower", extent=[0, c.lx, 0, c.lz], cmap="RdYlBu_r",
                      vmin=c.rb * 0.95, vmax=c.rt * 1.05)
        if np.abs(bx).max() + np.abs(bz).max() > 0:
            a.streamplot(xs[::st], zs[::st], bx[::st, ::st].T, bz[::st, ::st].T,
                         color="k", linewidth=0.5, density=1.1, arrowsize=0.5)
        a.set_title(f"t = {t:.2f}")
        a.set_xlim(0, c.lx)
        a.set_ylim(0, c.lz)
        a.set_xlabel("x")
        a.set_ylabel("z")
    for a in ax.ravel()[n:]:
        a.axis("off")
    fig.colorbar(im, ax=ax, shrink=0.8, label=r"density $\rho$")
    fig.suptitle(f"B0 = {c.b0}, Bz={c.bz}, modes={c.md}, profile = {c.pr}, {c.nx}x{c.nz} cells")
    fig.savefig(fn, dpi=140, bbox_inches="tight")
    plt.close(fig)


def pld(h, s, fn):
    c = s.c
    t = np.array([d["t"] for d in h])
    A = np.array([d["am"] for d in h])
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))

    for j, m in enumerate(c.md):
        ax[0, 0].semilogy(t, A[:, j] + 1e-12, label=f"m={m}")
    ax[0, 0].set(xlabel="t", ylabel=r"$|\hat\eta_m|$", title="interface mode amplitudes")
    ax[0, 0].legend(fontsize=8)
    ax[0, 0].grid(alpha=0.3)

    ax[0, 1].plot(t, [d["hm"] for d in h])
    ax[0, 1].set(xlabel="t", ylabel=r"$h_{mix}$", title="mixing-layer width")
    ax[0, 1].grid(alpha=0.3)

    kk = np.linspace(0.01, 2 * np.pi * max(c.md) * 1.1, 400)
    g2 = g2t(kk, c)
    ax[1, 0].plot(kk, g2, "k-", label=r"theory ($k\parallel B_0$)")
    ax[1, 0].axhline(0, color="gray", lw=0.5)
    for j, m in enumerate(c.md):
        k = 2 * np.pi * m / c.lx
        gr = fg(t, A[:, j], k)
        if np.isfinite(gr):
            ax[1, 0].plot(k, gr**2, "ro", label="simulation" if j == 0 else None)
    if c.b0 > 0:
        kc = c.g * (c.rt - c.rb) / (2 * c.b0**2)
        ax[1, 0].axvline(kc, color="b", ls="--", label=rf"$k_c$ = {kc:.2f}")
    ax[1, 0].set(xlabel="k", ylabel=r"$\gamma^2$", title=r"$\gamma^2$ vs k",
                 ylim=(min(g2.min(), -1) if c.b0 > 0 else -0.5, g2.max() * 1.3 + 0.5))
    ax[1, 0].legend(fontsize=8)
    ax[1, 0].grid(alpha=0.3)

    ax[1, 1].semilogy(t, [d["db"] + 1e-16 for d in h])
    ax[1, 1].set(xlabel="t", title=r"max|div B| $\Delta$/|B|")
    ax[1, 1].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(fn, dpi=140)
    plt.close(fig)


def plc(s, fn):
    c = s.c
    m = bk(s.mx(), c.cg)
    fig, ax = plt.subplots(figsize=(5, 4.5))
    im = ax.imshow(m.T, origin="lower", extent=[0, c.lx, 0, c.lz], cmap="viridis",
                   vmin=0, vmax=1)
    ax.set(title=f"{c.cg}x{c.cg} mixing fraction (t={s.t:.2f})", xlabel="x", ylabel="z")
    fig.colorbar(im, label="f")
    fig.savefig(fn, dpi=140, bbox_inches="tight")
    plt.close(fig)


def run(c):
    os.makedirs(c.od, exist_ok=True)
    tg = c.tg or f"B{c.b0:g}_{c.bz:g}_{c.md}_{c.pr}_{c.nx}x{c.nz}"
    s = Sv(c)

    if c.b0 > 0:
        kc = c.g * (c.rt - c.rb) / (2 * c.b0**2)
        print(f"k_c = {kc:.3f}  (m_c = {kc / (2 * np.pi):.2f}),  k_max = {kc / 2:.3f}")
    print(f"grid {c.nx}x{c.nz}, interface {c.ic} cells, At = {(c.rt - c.rb) / (c.rt + c.rb):.3f}")

    st = np.linspace(0.0, c.te, c.ns)
    i = 0
    rd = c.te / 400.0
    nr = 0.0
    h, sn = [], []
    t0 = tm.time()
    n = 0
    while s.t < c.te - 1e-12:
        dt = min(s.cdt(), c.te - s.t)
        if s.t >= nr - 1e-12:
            h.append(dg(s))
            nr += rd
        if i < len(st) and s.t >= st[i] - 1e-12:
            W = s.pm()
            sn.append((s.t, W[0].copy(), W[4].copy(), W[5].copy()))
            i += 1
        s.stp(dt)
        n += 1
        if n % 100 == 0:
            if not np.isfinite(s.U).all():
                raise RuntimeError(f"NaN/Inf at step {n}, t = {s.t:.4f}")
            d = h[-1]
            print(f"step {n:6d}  t = {s.t:7.4f}  dt = {dt:.2e}  vz_max = {d['vm']:.3f}  "
                  f"h_mix = {d['hm']:.4f}  divB = {d['db']:.1e}  [{tm.time() - t0:.0f}s]")
    h.append(dg(s))
    W = s.pm()
    sn.append((s.t, W[0].copy(), W[4].copy(), W[5].copy()))

    pls(sn, c, os.path.join(c.od, f"snapshots_{tg}.png"))
    pld(h, s, os.path.join(c.od, f"diagnostics_{tg}.png"))
    plc(s, os.path.join(c.od, f"coarse_map_{tg}.png"))
    np.savez_compressed(os.path.join(c.od, f"final_{tg}.npz"), x=s.x, z=s.z, t=s.t,
                        rho=W[0], vx=W[1], vz=W[2], p=W[3], bx=W[4], bz=W[5])

    t = np.array([d["t"] for d in h])
    A = np.array([d["am"] for d in h])
    print("\nmode   k        gamma_sim   gamma_theory")
    for j, m in enumerate(c.md):
        k = 2 * np.pi * m / c.lx
        gr = fg(t, A[:, j], k)
        g2 = g2t(k, c)
        gt = np.sqrt(g2) if g2 > 0 else float("nan")
        print(f"{m:4d}  {k:7.3f}  {gr:10.3f}  {gt:10.3f}" + ("   (stable)" if g2 <= 0 else ""))
    print(f"\ndone in {tm.time() - t0:.0f}s, outputs in '{c.od}/'")
    return s, h


def prs():
    p = ap.ArgumentParser()
    p.add_argument("--nx", type=int, default=256)
    p.add_argument("--nz", type=int, default=256)
    p.add_argument("--b0", type=float, default=0.2)
    p.add_argument("--bz", type=float, default=0.0)
    p.add_argument("--pr", choices=["uniform", "exponential"], default="uniform")
    p.add_argument("--lb", type=float, default=1.0)
    p.add_argument("--g", type=float, default=1.0)
    p.add_argument("--rt", type=float, default=2.0)
    p.add_argument("--rb", type=float, default=1.0)
    p.add_argument("--ic", type=int, default=16)
    p.add_argument("--cg", type=int, default=16)
    p.add_argument("--md", type=str, default="1,2,4,8")
    p.add_argument("--am", type=float, default=0.005)
    p.add_argument("--te", type=float, default=3.0)
    p.add_argument("--cfl", type=float, default=0.4)
    p.add_argument("--od", type=str, default="rti_output")
    p.add_argument("--tg", type=str, default="")
    a = p.parse_args()
    return Cf(nx=a.nx, nz=a.nz, b0=a.b0, bz=a.bz, pr=a.pr, lb=a.lb, g=a.g, rt=a.rt,
              rb=a.rb, ic=a.ic, cg=a.cg, md=tuple(int(m) for m in a.md.split(",")),
              am=a.am, te=a.te, cfl=a.cfl, od=a.od, tg=a.tg)


if __name__ == "__main__":
    run(prs())