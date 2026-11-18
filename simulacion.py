"""
Simulación del modelo Leslie-Gower con caza cooperativa y efecto Allee
Dinámica depredador-presa Lobo (Canis lupus) - Bisonte (Bison bison)

Modelo:
    dx/dt = r x (1 - x/K)(x - m) - (q + a y) x y      (bisontes)
    dy/dt = s y (1 - y/(x + c))                         (lobos)

Uso:
    pip install -r requirements.txt
    python simulacion.py

Genera las figuras en la carpeta figuras/ e imprime los equilibrios y su estabilidad.
Los parámetros son ilustrativos (densidades relativas, escaladas respecto a K).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------- parámetros
BASE = dict(
    r=1.0,    # tasa intrínseca de crecimiento del bisonte
    K=1.0,    # capacidad de carga del bisonte
    m=0.15,   # umbral del efecto Allee
    q=0.1,    # tasa base de ataque sin cooperación
    a=0.3,    # coeficiente de caza cooperativa
    s=0.3,    # tasa de crecimiento del lobo
    c=0.05,   # alimento alternativo del lobo
)
X0, Y0 = 0.8, 0.1   # condición inicial (bisontes, lobos)
OUT = "figuras"

BISONTE, LOBO, GRIS = "#A44A22", "#2F5D7C", "#5B6670"


# ------------------------------------------------------------------- modelo
def modelo(t, z, p):
    x, y = max(z[0], 0.0), max(z[1], 0.0)
    dx = p["r"] * x * (1 - x / p["K"]) * (x - p["m"]) - (p["q"] + p["a"] * y) * x * y
    dy = p["s"] * y * (1 - y / (x + p["c"]))
    return [dx, dy]


def simular(p, x0=X0, y0=Y0, T=300, n=3000):
    t = np.linspace(0, T, n)
    sol = solve_ivp(modelo, [0, T], [x0, y0], args=(p,), t_eval=t, rtol=1e-8, atol=1e-10)
    return sol.t, np.clip(sol.y, 0, None)


def rk4(p, x0=X0, y0=Y0, T=300, dt=0.05):
    """Runge-Kutta de 4.º orden escrito a mano (el mismo que usa la página web)."""
    n = int(T / dt)
    xs, ys = np.empty(n + 1), np.empty(n + 1)
    x, y = x0, y0
    xs[0], ys[0] = x, y
    f = lambda x, y: modelo(0, [x, y], p)
    for i in range(1, n + 1):
        k1 = f(x, y)
        k2 = f(x + dt / 2 * k1[0], y + dt / 2 * k1[1])
        k3 = f(x + dt / 2 * k2[0], y + dt / 2 * k2[1])
        k4 = f(x + dt * k3[0], y + dt * k3[1])
        x = max(x + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]), 0)
        y = max(y + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]), 0)
        xs[i], ys[i] = x, y
    return np.linspace(0, T, n + 1), xs, ys


# -------------------------------------------------------------- equilibrios
def isoclina_bisonte(x, p):
    """Rama no trivial: a y^2 + q y - r(1 - x/K)(x - m) = 0."""
    g = p["r"] * (1 - x / p["K"]) * (x - p["m"])
    if g < 0:
        return np.nan
    if p["a"] < 1e-12:
        return g / p["q"]
    return (-p["q"] + np.sqrt(p["q"] ** 2 + 4 * p["a"] * g)) / (2 * p["a"])


def jacobiana(x, y, p, e=1e-6):
    f0 = np.array(modelo(0, [x, y], p))
    fx = np.array(modelo(0, [x + e, y], p))
    fy = np.array(modelo(0, [x, y + e], p))
    return np.column_stack([(fx - f0) / e, (fy - f0) / e])


def equilibrios(p):
    E = [(0.0, 0.0, "Extinción total"), (0.0, p["c"], "Sin bisontes (lobos con alimento alternativo)")]
    h = lambda x: isoclina_bisonte(x, p) - (x + p["c"])
    xs = np.linspace(p["m"] + 1e-9, p["K"] - 1e-9, 4001)
    hs = np.array([h(x) for x in xs])
    for i in range(len(xs) - 1):
        if np.isfinite(hs[i]) and np.isfinite(hs[i + 1]) and hs[i] * hs[i + 1] <= 0:
            lo, hi = xs[i], xs[i + 1]
            for _ in range(60):
                mid = (lo + hi) / 2
                if h(lo) * h(mid) <= 0:
                    hi = mid
                else:
                    lo = mid
            xe = (lo + hi) / 2
            E.append((xe, xe + p["c"], "Coexistencia"))
    out = []
    for x, y, nombre in E:
        J = jacobiana(x, y, p)
        tr, det = np.trace(J), np.linalg.det(J)
        tipo = "Estable" if (det > 0 and tr < 0) else ("Punto silla" if det < 0 else "Inestable")
        out.append((nombre, x, y, tr, det, tipo))
    return out


# ------------------------------------------------------------------ figuras
def estilo():
    plt.rcParams.update({"font.size": 13, "axes.spines.top": False, "axes.spines.right": False})


def fig_escenarios():
    fig, axs = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
    for ax, (a, titulo) in zip(axs, [(0, "Sin cooperación (a = 0)"),
                                     (0.3, "Cooperación moderada (a = 0,3)"),
                                     (0.8, "Cooperación excesiva (a = 0,8)")]):
        t, (x, y) = simular({**BASE, "a": a})
        ax.axhline(BASE["m"], color=BISONTE, ls=":", lw=1.2, label="Umbral m")
        ax.plot(t, x, color=BISONTE, lw=2.5, label="Bisontes")
        ax.plot(t, y, color=LOBO, lw=2.5, ls="--", label="Lobos")
        ax.set_title(titulo, loc="left")
        ax.set_xlabel("Tiempo")
        ax.set_ylim(0, 1.05)
        ax.grid(alpha=0.3)
    axs[0].set_ylabel("Densidad relativa")
    axs[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(f"{OUT}/escenarios.png", dpi=150)


def fig_fases(p=BASE):
    fig, ax = plt.subplots(figsize=(9, 7))
    X, Y = np.meshgrid(np.linspace(0.01, 1.1, 22), np.linspace(0.01, 1.1, 22))
    U = np.vectorize(lambda x, y: modelo(0, [x, y], p)[0])(X, Y)
    V = np.vectorize(lambda x, y: modelo(0, [x, y], p)[1])(X, Y)
    L = np.hypot(U, V) + 1e-12
    ax.quiver(X, Y, U / L, V / L, color="#BDB6A6", scale=40, width=0.0025)
    xs = np.linspace(p["m"], p["K"], 400)
    ax.plot(xs, [isoclina_bisonte(x, p) for x in xs], color=BISONTE, ls="--", lw=2, label="Isoclina del bisonte")
    xx = np.linspace(0, 1.1, 10)
    ax.plot(xx, xx + p["c"], color=LOBO, ls="--", lw=2, label="Isoclina del lobo")
    for x0, y0 in [(0.9, 0.05), (1.0, 0.9), (0.4, 1.0), (0.3, 0.05), (0.12, 0.3), (0.25, 0.6)]:
        t, (x, y) = simular(p, x0, y0, T=400)
        col = "#1C2328" if x[-1] > 0.01 else "#8E8E8E"
        ax.plot(x, y, color=col, lw=1.6)
        ax.plot(x0, y0, "o", color=col, ms=5)
    for nombre, x, y, tr, det, tipo in equilibrios(p):
        ax.plot(x, y, "o", ms=12, mfc=BISONTE if tipo == "Estable" else "white", mec=BISONTE, mew=2, zorder=5)
    ax.set_xlim(0, 1.1)
    ax.set_ylim(0, 1.1)
    ax.set_xlabel("Bisontes (x)")
    ax.set_ylabel("Lobos (y)")
    ax.set_title(f"Plano de fases (a = {p['a']})", loc="left")
    ax.legend(frameon=False, loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{OUT}/plano_fases.png", dpi=150)


def fig_bifurcacion():
    A = np.linspace(0, 1.2, 121)
    fx, fy = [], []
    for a in A:
        t, (x, y) = simular({**BASE, "a": a}, T=600)
        fx.append(x[-1])
        fy.append(y[-1])
    fx = np.array(fx)
    a_crit = A[np.argmax(fx < 0.01)] if np.any(fx < 0.01) else None
    fig, ax = plt.subplots(figsize=(14, 5.5))
    ax.axhline(BASE["m"], color=BISONTE, ls=":", lw=1.2, label="Umbral m")
    ax.plot(A, fx, color=BISONTE, lw=3, label="Bisontes")
    ax.plot(A, fy, color=LOBO, lw=3, ls="--", label="Lobos")
    if a_crit is not None:
        ax.axvline(a_crit, color=GRIS, ls=":", lw=1.2)
        ax.text(a_crit + 0.01, 0.9, f"a crítico ≈ {a_crit:.2f}", color=GRIS)
    ax.set_xlabel("Coeficiente de cooperación a")
    ax.set_ylabel("Densidad final")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{OUT}/bifurcacion.png", dpi=150)
    return a_crit


def fig_caza():
    """Curvas reconstruidas con los valores de MacNulty et al. (OR por lobo y saturación)."""
    logit = lambda p: np.log(p / (1 - p))
    inv = lambda z: 1 / (1 + np.exp(-z))
    n = np.linspace(1, 20, 300)
    ataque = inv(logit(0.15) + np.log(1.67) * (np.minimum(n, 4) - 1))
    captura = inv(logit(0.01) + np.log(1.40) * (np.minimum(n, 11) - 1))
    fig, ax = plt.subplots(figsize=(14, 5.5))
    ax.axvspan(8.3, 12.9, color=LOBO, alpha=0.1, lw=0, label="Manada Mollie's (IC 95 %)")
    ax.plot(n, ataque, color="#C08A2E", lw=3, label="Iniciar el ataque (OR = 1,67, satura en 4)")
    ax.plot(n, captura, color=BISONTE, lw=3, label="Capturar al bisonte (OR = 1,40, satura en 11)")
    ax.set_xlabel("Lobos que participan en la cacería")
    ax.set_ylabel("Probabilidad")
    ax.set_xlim(1, 20)
    ax.set_ylim(0, 0.6)
    ax.legend(frameon=False, loc="lower right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{OUT}/exito_caza.png", dpi=150)


# --------------------------------------------------------------------- main
if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    estilo()

    print("Equilibrios con los parámetros base:", BASE)
    print(f"{'Equilibrio':48s} {'x*':>7s} {'y*':>7s} {'traza':>8s} {'det':>8s}  tipo")
    for nombre, x, y, tr, det, tipo in equilibrios(BASE):
        print(f"{nombre:48s} {x:7.3f} {y:7.3f} {tr:8.3f} {det:8.3f}  {tipo}")

    print()
    print("Estado final (t = 300) por escenario:")
    for a in (0, 0.3, 0.8):
        _, xs, ys = rk4({**BASE, "a": a})
        print(f"  a = {a:<4}  bisontes = {xs[-1]:.3f}  lobos = {ys[-1]:.3f}")

    fig_escenarios()
    fig_fases()
    a_crit = fig_bifurcacion()
    fig_caza()
    if a_crit is not None:
        print()
        print(f"Valor crítico de cooperación: a ≈ {a_crit:.2f}")
    print(f"Figuras guardadas en ./{OUT}/")
