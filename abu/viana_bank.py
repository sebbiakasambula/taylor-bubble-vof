"""Measured Taylor bubble Fr near T1-T3 from the appendix of Viana et al. (2003), Tables A1 and A2
(UnivCorrelation2.pdf pp. A-2 to A-12). Typed by hand from the page images, 13 Nov 2026. R in Viana = Nf here."""
import math, sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "scripts"))
from tb_references import fr_viana
g = 9.81
# Table A1 rows (source, fluid, D mm, mu mPa s, rho, sigma mN/m, U mm/s or None, Fr or None), page
A1 = [("Brown 1965 t1 (A-3)", "Glycerol solution", 8, 91, 1219, 65.3, 18.2, None),
      ("Brown 1965 t1 (A-4)", "Marcol", 26.4, 19.42, 865, 29.5, 170, 0.334),
      ("Brown 1965 t1 (A-4)", "Primol", 26.4, 142.3, 881, 30.5, 151, 0.296),
      ("Weber 1986 (A-4)", "100% corn syrup", 37.3, 6120, 1410, 87, None, 0.060),
      ("Weber 1986 (A-4)", "95% corn syrup", 37.3, 1830, 1340, 81, None, 0.158),
      ("Weber 1986 (A-4)", "88% corn syrup", 37.3, 518, 1320, 79, None, 0.273),
      ("Weber 1986 (A-4)", "66% sucrose", 37.3, 194, 1330, 79.1, None, 0.329),
      ("Weber 1986 (A-4)", "60% sucrose", 37.3, 51.1, 1280, 77.5, None, 0.342),
      ("Weber 1986 (A-4)", "88% corn syrup", 13.5, 518, 1320, 79, None, 0.059),
      ("Weber 1986 (A-4)", "66% sucrose", 13.5, 194, 1330, 79.1, None, 0.180),
      ("Weber 1986 (A-4)", "60% sucrose", 13.5, 51.1, 1280, 77.5, None, 0.270),
      ("Weber 1986 (A-4)", "Methanol", 13.5, 0.544, 787, 22.2, None, 0.307),
      ("Weber 1986 (A-4)", "66% sucrose", 22.1, 51.1, 1330, 79.1, None, 0.280),
      ("Weber 1986 (A-4)", "100% corn syrup", 22.1, 6120, 1410, 87, None, 0.025),
      ("Viana 2001 (A-4)", "Purolub 150 oil", 76.2, 424.4, 878, 31.5, 271, None),
      ("Viana this work (A-5)", "Silicon oil", 76.2, 1330, 969, 21.2, 228, None),
      ("Viana this work (A-5)", "Silicon oil", 76.2, 3834, 972, 21.3, 138, None),
      ("Viana this work (A-5)", "60% glycerin", 12.7, 19, 1157, 68.5, 101, None),
      ("Viana this work (A-5)", "75% glycerin", 12.7, 67, 1198, 66.4, 84, None)]
rows = []
for src, fl, D, mu, rho, sig, U, Fr in A1:
    D, mu, sig = D / 1e3, mu / 1e3, sig / 1e3
    Eo = rho * g * D ** 2 / sig; Nf = math.sqrt(D ** 3 * g * rho * rho) / mu
    if Fr is None: Fr = U / 1e3 / math.sqrt(g * D)
    rows.append((src + " computed", fl, Nf, Eo, Fr))
# Table A2 rows (source, fluid, R, Eo, Fr), only R 10 to 5000 and Eo 10 to 450
A2 = [("W&B fig2 (A-7)", "Aqueous ethanol", 1134, 18.94, 0.2754), ("W&B fig2 (A-7)", "Aqueous ethanol", 1642, 31.01, 0.3178),
      ("W&B fig2 (A-8)", "40% sucrose", 632.7, 14.23, 0.2414), ("W&B fig2 (A-8)", "40% sucrose", 1170, 32.3, 0.3226),
      ("W&B fig2 (A-8)", "Ethylene glycol", 147.1, 18.26, 0.2384), ("W&B fig2 (A-8)", "Ethylene glycol", 171.4, 22.38, 0.2672),
      ("W&B fig2 (A-8)", "Ethylene glycol", 312.5, 49.87, 0.326), ("W&B fig2 (A-8)", "Ethylene glycol", 440.4, 78.77, 0.3343),
      ("W&B fig2 (A-8)", "58% sucrose", 97, 16.08, 0.1942), ("W&B fig2 (A-8)", "58% sucrose", 175.7, 35.52, 0.2929),
      ("W&B fig2 (A-8)", "58% sucrose", 270.4, 63.11, 0.327), ("W&B fig2 (A-8)", "58% sucrose", 406.4, 108.6, 0.3319),
      ("W&B fig2 (A-8)", "Tellus oil", 49.28, 21.98, 0.17518), ("W&B fig2 (A-8)", "Tellus oil", 56.84, 26.58, 0.2073),
      ("W&B fig2 (A-8)", "Tellus oil", 103.7, 59.23, 0.2944), ("W&B fig2 (A-8)", "Tellus oil", 157.9, 103.8, 0.3178),
      ("W&B fig2 (A-8)", "Tellus oil", 247.2, 188.7, 0.3285),
      ("W&B fig2 (A-9)", "90% glycerol", 45.44, 40.87, 0.2311), ("W&B fig2 (A-9)", "90% glycerol", 70.18, 72.94, 0.2832),
      ("W&B fig2 (A-9)", "90% glycerol", 105.8, 126.1, 0.3153),
      ("W&B fig2 (A-9)", "95% glycerol", 21.69, 41.8, 0.1435), ("W&B fig2 (A-9)", "95% glycerol", 33.27, 73.94, 0.2204),
      ("W&B fig2 (A-9)", "95% glycerol", 50.34, 128.4, 0.2671),
      ("W&B fig2 (A-9)", "Voluta oil", 26.64, 112.1, 0.2039), ("W&B fig2 (A-9)", "Voluta oil", 41, 199.2, 0.2564),
      ("W&B fig2 (A-9)", "Voluta oil", 72.27, 424.3, 0.2949),
      ("W&B fig2 (A-9)", "Glycerol", 23.47, 134.4, 0.1849), ("W&B fig2 (A-9)", "Glycerol", 40.40, 277.2, 0.2516),
      ("Zukoski fig3 (A-10)", "Water-glycerin", 20.1, 99.59, 0.1373), ("Zukoski fig3 (A-10)", "Water-glycerin", 148, 77.46, 0.326),
      ("Zukoski fig3 (A-10)", "Water-glycerin", 46.98, 34.86, 0.2137), ("Zukoski fig3 (A-10)", "Water-glycerin", 99.868, 34.29, 0.2808),
      ("Zukoski fig3 (A-10)", "Water-glycerin", 233.4, 25.51, 0.3104), ("Zukoski fig3 (A-10)", "Water-glycerin", 81.21, 19.73, 0.2199),
      ("Zukoski fig3 (A-10)", "Water-glycerin", 182.3, 17.72, 0.2574), ("Zukoski fig3 (A-10)", "Water-glycerin", 238.8, 14.94, 0.2418),
      ("Wallis 1969 (A-11)", "50% glycerin", 777.8, 17.88, 0.3149), ("Wallis 1969 (A-11)", "50% glycerin", 1070, 27.37, 0.3427),
      ("Wallis 1969 (A-11)", "50% glycerin", 4094, 163.8, 0.375), ("Wallis 1969 (A-11)", "70% glycerin", 1254, 173.3, 0.3719),
      ("Wallis 1969 (A-11)", "70% glycerin", 256.8, 20.91, 0.2968), ("Wallis 1969 (A-11)", "89% glycerin", 179.6, 174.4, 0.3027),
      ("Wallis 1969 (A-11)", "89% glycerin", 80.88, 60.22, 0.245), ("Wallis 1969 (A-11)", "89% glycerin", 34.45, 19.3, 0.1371),
      ("Wallis 1969 (A-11)", "95% glycerin", 89.04, 235, 0.2873), ("Wallis 1969 (A-11)", "95% glycerin", 35.14, 68.05, 0.1855),
      ("Wallis 1969 (A-11)", "95% glycerin", 22.94, 38.53, 0.1379)]
rows += [(s, f, R, E, F) for s, f, R, E, F in A2]
targets = {"T1": (1000, 300), "T2": (50, 100), "T3": (500, 20)}
for k, (N0, E0) in targets.items():
    print(f"\n{k} (Nf {N0}, Eo {E0}), Viana at target {fr_viana(E0, N0):.4f}")
    d = sorted(rows, key=lambda r: math.hypot(math.log(r[2] / N0), math.log(r[3] / E0)))
    for s, f, N, E, F in d[:7]:
        v = fr_viana(E, N)
        print(f"  {s:28s} {f:18s} Nf {N:8.1f} ({100*(N/N0-1):+5.0f}%) Eo {E:6.1f} ({100*(E/E0-1):+5.0f}%) Fr meas {F:.4f} Viana {v:.4f} meas/Viana {100*(F/v-1):+5.1f}%")
