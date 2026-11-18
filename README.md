# Dinámica depredador–presa Lobo–Bisonte

Proyecto final de Modelación Matemática. Simulación del modelo de **Leslie-Gower con caza cooperativa y efecto Allee**, aplicado al lobo (*Canis lupus*) y al bisonte (*Bison bison*) en el Parque Nacional de Yellowstone.

**Simulación en vivo:** [abrir la simulación](https://dylan123412.github.io/ProyectoModelacionMatematica/)

## El modelo

```
dx/dt = r·x·(1 − x/K)·(x − m) − (q + a·y)·x·y     (bisontes)
dy/dt = s·y·(1 − y/(x + c))                       (lobos)
```

| Símbolo | Significado | Valor base |
| --- | --- | --- |
| x | Densidad de bisontes (presa) | x₀ = 0,8 |
| y | Densidad de lobos (depredador) | y₀ = 0,1 |
| r | Tasa de crecimiento del bisonte | 1 |
| K | Capacidad de carga del bisonte | 1 |
| m | Umbral del efecto Allee | 0,15 |
| q | Tasa base de ataque sin cooperación | 0,1 |
| a | Coeficiente de caza cooperativa | 0,3 |
| s | Tasa de crecimiento del lobo | 0,3 |
| c | Alimento alternativo del lobo | 0,05 |

Las densidades son relativas (escaladas respecto a K). Los parámetros son ilustrativos: reproducen el comportamiento descrito en la literatura, no se midieron en campo.

## Resultados principales

| Cooperación a | Bisontes (final) | Lobos (final) | Resultado |
| --- | --- | --- | --- |
| 0 | 0,872 | 0,922 | Coexistencia |
| 0,3 | 0,577 | 0,627 | Coexistencia con menos bisontes |
| 0,8 | 0 | 0,05 | Colapso del bisonte |

A partir de **a ≈ 0,57** el bisonte cae por debajo del umbral m y colapsa de forma repentina.

## Archivos

| Archivo | Qué es |
| --- | --- |
| `index.html` | Simulación interactiva (deslizadores, gráfica en el tiempo, plano de fases, equilibrios, diagrama de bifurcación y datos de campo) |
| `simulacion.py` | Script de Python que resuelve el modelo, calcula los equilibrios y genera las figuras |
| `requirements.txt` | Librerías de Python necesarias |

## Cómo ejecutarlo

**Simulación web:** abre `index.html` en el navegador, o usa el enlace de GitHub Pages de arriba.

**Python:**

```bash
pip install -r requirements.txt
python simulacion.py
```

El script imprime los puntos de equilibrio con su estabilidad, el estado final de cada escenario y el valor crítico de cooperación, y guarda las figuras en `figuras/`.

## Métodos

- Integración numérica con Runge-Kutta de 4.º orden (Δt = 0,05) y `solve_ivp` de SciPy.
- Equilibrios: cruce de isoclinas, encontrado por bisección.
- Estabilidad: traza y determinante de la matriz jacobiana.
- Bifurcación: densidad final para 121 valores de a entre 0 y 1,2.

## Referencias

- MacNulty, D. R., Tallian, A., Stahler, D. R. y Smith, D. W. (2014). *Influence of group size on the success of wolves hunting bison.*
- Teixeira Alves, M. y Hilker, F. M. (2017). Hunting cooperation and Allee effects in predators. *Journal of Theoretical Biology.*
- González-Olivares, E. y Rojas-Palma, A. Modelo de Leslie-Gower con efecto Allee y caza cooperativa.
- Lotka (1925), Volterra (1926), Allee (1931), Leslie y Gower (1960), Holling (1959).
