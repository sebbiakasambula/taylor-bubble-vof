"""Dimensionless groups to interFoam inputs for the Taylor bubble. The only place this is done.

Scales: pipe diameter D = 1, liquid density rho_l = 1, gravity g = 1. The velocity scale is then
sqrt(g D) = 1, so the Froude number Fr = U_b / sqrt(g D) is simply the rise velocity U_b.

Definitions, with the density difference drho = rho_l - rho_g driving buoyancy:
    Eo  = drho g D^2 / sigma                       Eotvos number, buoyancy against surface tension
    Nf  = sqrt(drho rho_l g D^3) / mu_l            inverse viscosity number, buoyancy against viscosity
    rho_r = rho_l / rho_g                          density ratio
    mu_r  = mu_l / mu_g                            viscosity ratio, held fixed (see MU_R)

Using drho rather than rho_l changes Eo and Nf by under 1% for rho_r >= 100, but it is the
physically correct choice and matches the usual Taylor-bubble literature. State it in the write-up.
"""
import math

D = 1.0
RHO_L = 1.0
G = 1.0
R = D / 2.0

MU_R = 100.0   # held fixed; air-water is about 55 and oil-gas higher. Not a sweep parameter.


def fluid(Eo, Nf, rho_r, mu_r=MU_R):
    """interFoam transport inputs for one (Eo, Nf, rho_r) point."""
    rho_g = RHO_L / rho_r
    drho = RHO_L - rho_g
    sigma = drho * G * D ** 2 / Eo
    mu_l = math.sqrt(drho * RHO_L * G * D ** 3) / Nf
    mu_g = mu_l / mu_r
    return dict(Eo=Eo, Nf=Nf, rho_r=rho_r, mu_r=mu_r,
                rho_l=RHO_L, rho_g=rho_g, sigma=sigma,
                mu_l=mu_l, mu_g=mu_g,
                nu_l=mu_l / RHO_L, nu_g=mu_g / rho_g)


def morton(f):
    """Morton number Mo = g mu_l^4 drho / (rho_l^2 sigma^3). Depends on fluid only, not on D."""
    drho = f["rho_l"] - f["rho_g"]
    return G * f["mu_l"] ** 4 * drho / (f["rho_l"] ** 2 * f["sigma"] ** 3)


if __name__ == "__main__":
    for Eo, Nf, rr in ((100, 2000, 1000), (10, 10, 100), (500, 500, 500)):
        f = fluid(Eo, Nf, rr)
        print(f"Eo={Eo:4g} Nf={Nf:5g} rho_r={rr:5g} -> sigma={f['sigma']:.5f} "
              f"nu_l={f['nu_l']:.3e} nu_g={f['nu_g']:.3e} rho_g={f['rho_g']:.4f} Mo={morton(f):.2e}")
