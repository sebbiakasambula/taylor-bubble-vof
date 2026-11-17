"""Write a complete interFoam case for one Taylor bubble. Solver dictionaries follow the ESI v2412
damBreak tutorial (an unmodified copy of the tutorial case), and differ only where this
problem differs.

Coordinates: the pipe axis is x, pointing UP, so gravity is (-1 0 0). Radius is y. The wedge is in z.
This is the same 6-vertex collapsed-axis wedge verified in Project 1 (checkMesh clean, and the official
build reproduced Hagen-Poiseuille to 0.009%).

Departures from damBreak, each deliberate:
  top boundary   inletOutlet with inletValue 1, because this boundary sits in LIQUID. damBreak's open
                 boundary sits in air and uses 0; copying it would draw gas into the pipe.
  initial bubble radius set from the estimated developed film thickness, so the bubble starts near its
                 final shape and the start-up transient is short (the Project 1 lesson: start close).
"""
import math
from pathlib import Path

import numpy as np

from mesh_design import radial_cells
from tb_params import R, fluid
from tb_references import FR_DUMITRESCU, film_thickness, fr_viana

HDR = "FoamFile {{ version 2.0; format ascii; class {cls}; object {obj}; }}\n"
TH = math.radians(2.5)
YW, ZW = R * math.cos(TH), R * math.sin(TH)


def block_mesh(L=20.0, nx=480, nr=96, grading=0.25):
    """Axis along x from the closed bottom (x = 0) to the open top (x = L). grading < 1 is wall-fine."""
    return HDR.format(cls="dictionary", obj="blockMeshDict") + f"""
scale 1;
vertices
(
    (0    0     0)          // 0 axis, bottom
    ({L} 0     0)          // 1 axis, top
    ({L} {YW:.9f} -{ZW:.9f})   // 2 wall, top, front
    (0    {YW:.9f} -{ZW:.9f})   // 3 wall, bottom, front
    ({L} {YW:.9f}  {ZW:.9f})   // 4 wall, top, back
    (0    {YW:.9f}  {ZW:.9f})   // 5 wall, bottom, back
);
blocks ( hex (0 1 2 3 0 1 4 5) ({nx} {nr} 1) simpleGrading (1 {grading} 1) );
boundary
(
    bottom {{ type wall;  faces ((0 0 5 3)); }}
    top    {{ type patch; faces ((1 2 4 1)); }}
    wall   {{ type wall;  faces ((3 5 4 2)); }}
    front  {{ type wedge; faces ((0 3 2 1)); }}
    back   {{ type wedge; faces ((0 1 4 5)); }}
);
"""


def capillary_dt(f, dx):
    """Explicit surface-tension time-step limit (Brackbill et al. 1992): dt < sqrt(rho_mean dx^3 / (2 pi sigma)),
    with rho_mean the mean of the two densities and dx the smallest cell size at the interface."""
    rho_mean = 0.5 * (f["rho_l"] + f["rho_g"])
    return math.sqrt(rho_mean * dx ** 3 / (2.0 * math.pi * f["sigma"]))


def control_dict(end_time, n_snap=64, extra_functions="", max_dt=0.05):
    return HDR.format(cls="dictionary", obj="controlDict") + f"""
application     interFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {end_time:.6g};
deltaT          1e-4;
writeControl    adjustable;
writeInterval   {end_time / n_snap:.6g};
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   8;
runTimeModifiable yes;
adjustTimeStep  yes;
maxCo           0.5;
maxAlphaCo      0.2;
maxDeltaT       {max_dt:.6g};

functions
{{
    // liquid volume every 20 steps, for the mass-conservation check (V2)
    liquidVolume
    {{
        type            volFieldValue;
        libs            (fieldFunctionObjects);
        fields          (alpha.water);
        operation       volIntegrate;
        regionType      all;
        writeFields     false;
        writeControl    timeStep;
        writeInterval   20;
    }}
{extra_functions}}}
"""


FV_SCHEMES = HDR.format(cls="dictionary", obj="fvSchemes") + """
ddtSchemes           { default Euler; }
gradSchemes          { default Gauss linear; }
divSchemes
{
    div(rhoPhi,U)   Gauss linearUpwind grad(U);
    div(phi,alpha)  Gauss vanLeer;
    div(phirb,alpha) Gauss linear;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes     { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes        { default corrected; }
"""

FV_SOLUTION = HDR.format(cls="dictionary", obj="fvSolution") + """
solvers
{
    "alpha.water.*"
    {
        nAlphaCorr      2;
        nAlphaSubCycles 1;
        cAlpha          1;
        MULESCorr       yes;
        nLimiterIter    5;
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0;
    }
    "pcorr.*"   { solver PCG; preconditioner DIC; tolerance 1e-5;  relTol 0; }
    p_rgh       { solver PCG; preconditioner DIC; tolerance 1e-07; relTol 0.05; }
    p_rghFinal  { $p_rgh; relTol 0; }
    U           { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-06; relTol 0; }
}
PIMPLE
{
    momentumPredictor   no;
    nOuterCorrectors    1;
    nCorrectors         3;
    nNonOrthogonalCorrectors 0;
}
relaxationFactors { equations { ".*" 1; } }
"""


def transport(f):
    return HDR.format(cls="dictionary", obj="transportProperties") + f"""
phases          (water air);
water {{ transportModel Newtonian; nu {f['nu_l']:.10g}; rho {f['rho_l']:.10g}; }}
air   {{ transportModel Newtonian; nu {f['nu_g']:.10g}; rho {f['rho_g']:.10g}; }}
sigma           {f['sigma']:.10g};
"""


TURBULENCE = HDR.format(cls="dictionary", obj="turbulenceProperties") + "\nsimulationType  laminar;\n"
GRAVITY = HDR.format(cls="uniformDimensionedVectorField", obj="g") + \
    "\ndimensions      [0 1 -2 0 0 0 0];\nvalue           (-1 0 0);\n"

WEDGES = "    front { type wedge; }\n    back  { type wedge; }\n"
ALPHA = HDR.format(cls="volScalarField", obj="alpha.water") + """
dimensions      [0 0 0 0 0 0 0];
internalField   uniform 1;
boundaryField
{
    bottom { type zeroGradient; }
    wall   { type zeroGradient; }
    top    { type inletOutlet; inletValue uniform 1; value uniform 1; }
""" + WEDGES + "}\n"
VELOCITY = HDR.format(cls="volVectorField", obj="U") + """
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 0);
boundaryField
{
    bottom { type noSlip; }
    wall   { type noSlip; }
    top    { type pressureInletOutletVelocity; value uniform (0 0 0); }
""" + WEDGES + "}\n"
P_RGH = HDR.format(cls="volScalarField", obj="p_rgh") + """
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform 0;
boundaryField
{
    bottom { type fixedFluxPressure; value uniform 0; }
    wall   { type fixedFluxPressure; value uniform 0; }
    top    { type totalPressure; p0 uniform 0; }
""" + WEDGES + "}\n"


def set_fields(x0, length, radius):
    """Gas cylinder from x0 to x0 + length with a hemispherical nose; liquid everywhere else."""
    return HDR.format(cls="dictionary", obj="setFieldsDict") + f"""
defaultFieldValues ( volScalarFieldValue alpha.water 1 );
regions
(
    cylinderToCell
    {{
        p1 ({x0:.6g} 0 0);
        p2 ({x0 + length:.6g} 0 0);
        radius {radius:.6g};
        fieldValues ( volScalarFieldValue alpha.water 0 );
    }}
    sphereToCell
    {{
        centre ({x0 + length:.6g} 0 0);
        radius {radius:.6g};
        fieldValues ( volScalarFieldValue alpha.water 0 );
    }}
);
"""


def write_case(dest, Eo, Nf, rho_r, L=20.0, nx=480, nr=96, grading=0.25,
               x0=1.5, bubble_length=3.0, travel=12.0, n_snap=64, dt_factor=1.0):
    """Write the case. The run time is the time for the nose to travel `travel` diameters at 90% of the
    Viana et al. (2003) eq (23) speed, the 90% covering their scatter so a slow bubble still gets a
    terminal window. The outlet check uses the fastest plausible bubble, 1 / 0.9 times Viana, so the nose
    stays 1 D below the outlet even then. maxDeltaT is the explicit capillary limit, eq. (61) of Brackbill et
    al. (1992), on the radial cell where the film interface sits (r = initial bubble radius), never the wall
    cell. Eq. (61) comes from eq. (59), c dt / dx < 1/2, which already carries the safety factor of 1/2, so
    production runs take it as published (dt_factor = 1). The V3 time-step study passes dt_factor = 0.5.
    """
    f = fluid(Eo, Nf, rho_r)
    delta = film_thickness(FR_DUMITRESCU, Nf)
    radius = max(R - delta, 0.25)                       # start near the developed shape
    nose = x0 + bubble_length + radius
    if nose + travel / 0.9 > L - 1.0:
        raise ValueError(f"bubble could reach within 1 D of the outlet: nose {nose:.2f} + fastest travel "
                         f"{travel / 0.9:.2f} > {L - 1:.1f}. Lengthen the pipe or shorten the travel.")
    end_time = travel / (0.9 * fr_viana(Eo, Nf))
    # The film interface sits at r = radius, so that radial cell (not the wall cell) sets the limit. No extra
    # factor on eq. (61), its 1/2 is already built in (see docstring).
    faces = radial_cells(nr, grading)[0]
    i = int(np.searchsorted(faces, radius)) - 1
    max_dt = dt_factor * capillary_dt(f, float(min(faces[i + 1] - faces[i], L / nx)))
    dest = Path(dest)
    for d in ("0", "constant", "system"):
        (dest / d).mkdir(parents=True, exist_ok=True)
    files = {
        "system/blockMeshDict": block_mesh(L, nx, nr, grading),
        "system/controlDict": control_dict(end_time, n_snap, max_dt=max_dt),
        "system/fvSchemes": FV_SCHEMES,
        "system/fvSolution": FV_SOLUTION,
        "system/setFieldsDict": set_fields(x0, bubble_length, radius),
        "constant/transportProperties": transport(f),
        "constant/turbulenceProperties": TURBULENCE,
        "constant/g": GRAVITY,
        "0/alpha.water": ALPHA,
        "0/U": VELOCITY,
        "0/p_rgh": P_RGH,
    }
    for name, text in files.items():
        (dest / name).write_text(text)
    return dict(fluid=f, film_estimate=delta, bubble_radius=radius, end_time=end_time,
                max_dt=max_dt, initial_nose=nose)


GRAVITY_ZERO = HDR.format(cls="uniformDimensionedVectorField", obj="g") + \
    "\ndimensions      [0 1 -2 0 0 0 0];\nvalue           (0 0 0);\n"

# largest velocity magnitude every 10 steps: the spurious-current history for the static-bubble test
MAX_SPEED = """    maxSpeed
    {
        type            fieldMinMax;
        libs            (fieldFunctionObjects);
        fields          (U);
        mode            magnitude;
        writeControl    timeStep;
        writeInterval   10;
    }
"""


def set_fields_sphere(centre, radius):
    return HDR.format(cls="dictionary", obj="setFieldsDict") + f"""
defaultFieldValues ( volScalarFieldValue alpha.water 1 );
regions
(
    sphereToCell
    {{
        centre ({centre:.6g} 0 0);
        radius {radius:.6g};
        fieldValues ( volScalarFieldValue alpha.water 0 );
    }}
);
"""


def write_static_case(dest, Eo, Nf, rho_r, L=4.0, nx=96, nr=96, grading=0.25,
                      radius=0.25, centre=2.0, n_capillary=10.0, n_snap=20):
    """V1: a spherical bubble at rest with gravity off, on the production cell sizes.

    Nothing should move. Any velocity that appears is a spurious current from the surface-tension
    discretisation, and the pressure jump across the interface must equal the Laplace value 2 sigma / R.
    The axial cell matches the production mesh (4 / 96 = 20 / 480 = 0.042 D), because spurious currents
    depend on the cell shape, so testing on any other mesh would not tell me about the production one.
    """
    f = fluid(Eo, Nf, rho_r)
    t_cap = math.sqrt(f["rho_l"] * radius ** 3 / f["sigma"])
    end_time = n_capillary * t_cap
    # With no flow, nothing else limits the step, so without this cap the solver grows it to maxDeltaT and
    # V1 would measure a time-step instability instead of spurious currents. Smallest cell the interface
    # crosses is the radial cell near r = radius; half of eq. (61), static test only (cheap, and the precision test).
    faces = radial_cells(nr, grading)[0]
    near = np.diff(faces)[(faces[:-1] < radius * 1.2)]
    dt_cap = 0.5 * capillary_dt(f, float(min(near.min(), L / nx)))
    dest = Path(dest)
    for d in ("0", "constant", "system"):
        (dest / d).mkdir(parents=True, exist_ok=True)
    files = {
        "system/blockMeshDict": block_mesh(L, nx, nr, grading),
        "system/controlDict": control_dict(end_time, n_snap, extra_functions=MAX_SPEED, max_dt=dt_cap),
        "system/fvSchemes": FV_SCHEMES,
        "system/fvSolution": FV_SOLUTION,
        "system/setFieldsDict": set_fields_sphere(centre, radius),
        "constant/transportProperties": transport(f),
        "constant/turbulenceProperties": TURBULENCE,
        "constant/g": GRAVITY_ZERO,
        "0/alpha.water": ALPHA,
        "0/U": VELOCITY,
        "0/p_rgh": P_RGH,
    }
    for name, text in files.items():
        (dest / name).write_text(text)
    return dict(fluid=f, t_capillary=t_cap, end_time=end_time, max_dt=dt_cap,
                steps_at_cap=end_time / dt_cap,
                laplace_jump=2.0 * f["sigma"] / radius, capillary_velocity=f["sigma"] / f["mu_l"])


if __name__ == "__main__":
    import sys
    info = write_case(sys.argv[1] if len(sys.argv) > 1 else "tb_test_case", Eo=100, Nf=500, rho_r=500)
    print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in info.items() if k != "fluid"})
