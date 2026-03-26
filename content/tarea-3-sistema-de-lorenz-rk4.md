---
title: "Inciso 1: Sistema de Lorenz con RK4"
---

# Inciso 1: Sistema de Lorenz con RK4

## Planteamiento
Resolver numericamente el sistema de Lorenz en la notacion del curso, usando las variables `W`, `T1` y `T2`, mediante un esquema de Runge-Kutta de cuarto orden (`RK4`). El inciso pide mantener `Pr = 10`, `b = 8/3`, `dt = 0.005`, integrar en `0 <= t <= 50` y partir de la condicion inicial `y_0 = (0, 0.5, 0.5)`.

El sistema integrado fue:

`dW/dt = Pr (T1 - W)`

`dT1/dt = -W T2 + r W - T1`

`dT2/dt = W T1 - b T2`

Estas ecuaciones son equivalentes al sistema de Lorenz clasico bajo el renombre `W <-> x`, `T1 <-> y` y `T2 <-> z`, pero en este capitulo se conservan los simbolos del enunciado para que la resolucion coincida exactamente con la guia.

## Enfoque de trabajo
La implementacion numerica es necesaria, pero no suficiente. El objetivo principal del inciso es analizar criticamente como cambia el comportamiento del sistema cuando `r` pasa de un regimen estable a uno transicional y, finalmente, a un regimen caotico. Por eso, cada simulacion se interpreta a partir de dos salidas complementarias:

- la evolucion temporal de `T1(t)` y `T2(t)`, que permite identificar estabilizacion u oscilaciones persistentes;
- la trayectoria en el espacio de fases `W-T1-T2`, que permite reconocer puntos fijos, cambios de lobo y estructura atractora.

## Objetivos del inciso
- Implementar `RK4` de forma explicita y reproducible.
- Resolver el sistema para `r = 2`, `10`, `24`, `25` y `30` en el intervalo `0 <= t <= 50`.
- Graficar la evolucion temporal de `T1(t)` y `T2(t)` para `r = 2`, `10`, `24` y `25`.
- Graficar la trayectoria en el espacio de fases `W, T1, T2` para `r = 10`, `24`, `25` y `30`.
- Mostrar la sensibilidad a condiciones iniciales para `r = 30` comparando `T1(t)` para `y_0 = (0, 0.5, 0.5)` y `y_0 = (0, 0.5, 0.50001)`.

## Entregables
- Notebook reproducible con la implementacion de `RK4` y las simulaciones pedidas.
- Documento de resultados y discusion para los distintos valores de `r`.
- Figuras de `T1(t)` y `T2(t)` y de espacio de fases `W-T1-T2` exportadas para el book.

## Resolucion y notebook asociado
- [Resolucion Inciso 1: Sistema de Lorenz con RK4](resultados-discusion-lorenz-rk4.md)
- [Analisis de sensibilidad: Sistema de Lorenz con RK4](analisis-sensibilidad-lorenz-rk4.md)
- [Notebook asociado: Sistema de Lorenz con RK4](../notebooks/Guia_practica_Lorenz_RK4.ipynb)

## Referencias base (opcional)
- Rutinas MATLAB de apoyo en `Tarea2ModAmb/dydt.m`, `Tarea2ModAmb/lorenz.m` y `Tarea2ModAmb/rkstep.m`.
- Lorenz, E. N. (1963). Deterministic nonperiodic flow. *Journal of the Atmospheric Sciences*, 20(2), 130-141.
