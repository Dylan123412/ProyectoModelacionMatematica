# Dinámica depredador–presa Lobo–Bisonte

Proyecto final de Modelación Matemática. Simulación del modelo de **Leslie-Gower con caza cooperativa y efecto Allee**, aplicado al lobo (*Canis lupus*) y al bisonte (*Bison bison*) en el Parque Nacional de Yellowstone.

**Simulación en vivo:** [abrir la simulación](https://dylan123412.github.io/ProyectoModelacionMatematica/)

## El modelo

```
dx/dt = r·x·(1 − x/K)·(x − m) − (q + a·y)·x·y     (bisontes)
dy/dt = s·y·(1 − y/(x + c))                       (lobos)
```

| Símbolo | Significado | Valor base | De dónde sale |
| --- | --- | --- | --- |
| x | Densidad de bisontes (presa), en fracciones de K | x₀ = 0,8 | — |
| y | Densidad de lobos (depredador) | y₀ = 0,1 | — |
| K | Capacidad de carga del bisonte | 1 (= 6.200 bisontes) | Plumb et al. (2009) |
| r | Tasa de crecimiento del bisonte | 1 por año | Crecimiento de 14,5 % anual con ~4.900 bisontes (Geremia, 2022): r·(1 − x)·(x − m) = ln 1,145 da r ≈ 1,01 |
| s | Tasa de crecimiento del lobo | 0,3 por año | De 21 lobos (1995) a 171 (2004): ≈ 0,23/año ya con competencia (Yellowstone Wolf Project) |
| m | Umbral del efecto Allee | 0,15 | Sin medición directa; se prueba en la sensibilidad |
| q | Tasa base de ataque sin cooperación | 0,1 | Pequeño: un lobo solo acierta 1 % (MacNulty et al., 2014) |
| a | Coeficiente de caza cooperativa | 0 a 1,2 | Parámetro estudiado |
| c | Alimento alternativo del lobo | 0,05 | Pequeño: en invierno el alce deja el Valle de Pelican |

El tiempo se mide en años. Como en González-Olivares y Rojas-Palma (2021) la ecuación del lobo es s·y·(1 − y/(n·x + c)); aquí n = 1 midiendo y en unidades de n·K lobos.

## Resultados principales

| Cooperación a | Bisontes (final) | Lobos (final) | Resultado |
| --- | --- | --- | --- |
| 0 | 0,872 | 0,922 | Coexistencia |
| 0,3 | 0,577 | 0,627 | Coexistencia con menos bisontes |
| 0,8 | 0 | 0,05 | Colapso del bisonte |

A partir de **a ≈ 0,57** el bisonte cae por debajo del umbral m y colapsa de forma repentina.

### Análisis de sensibilidad

Valor crítico de a al variar cada parámetro entre 50 % y 150 % de su valor base (los demás fijos):

| Parámetro | 50 % | 75 % | Base | 125 % | 150 % |
| --- | --- | --- | --- | --- | --- |
| r | 0,16 | 0,36 | 0,57 | 0,75 | 0,91 |
| m | 1,18 | 0,80 | 0,57 | 0,41 | 0,30 |
| q | 0,70 | 0,63 | 0,57 | 0,50 | 0,44 |
| c | 0,67 | 0,61 | 0,57 | 0,52 | 0,48 |
| s | 0,49 | 0,54 | 0,57 | 0,57 | 0,58 |

El colapso existe con todos los valores probados (incluso con m = 0,01, donde exige a ≈ 3,1). Los parámetros solo cambian cuánta cooperación hace falta; lo que más pesa es r y m.

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

- MacNulty, D. R., Tallian, A., Stahler, D. R. y Smith, D. W. (2014). Influence of group size on the success of wolves hunting bison. *PLoS ONE*, 9(11), e112884. https://doi.org/10.1371/journal.pone.0112884
- González-Olivares, E. y Rojas-Palma, A. (2021). Influencia del efecto Allee en las presas y de la colaboración entre los depredadores en un modelo de depredación del tipo Leslie-Gower. *Revista de Modelamiento Matemático de Sistemas Biológicos*, 1(2), 21–29.
- Teixeira Alves, M. y Hilker, F. M. (2017). Hunting cooperation and Allee effects in predators. *Journal of Theoretical Biology*, 419, 13–22.
- Plumb, G. E., White, P. J., Coughenour, M. B. y Wallen, R. L. (2009). Carrying capacity, migration, and dispersal in Yellowstone bison. *Biological Conservation*, 142(11), 2377–2387.
- Geremia, C. (2022). *Status Report on the Yellowstone Bison Population*. National Park Service.
- Yellowstone Wolf Project, informes anuales (National Park Service).
- Lotka (1925), Volterra (1926), Allee (1931), Leslie y Gower (1960), Holling (1959).
