---
title: "Inciso 1: Sistema de Lorenz con RK4"
---

# Inciso 1: Sistema de Lorenz con RK4

## Planteamiento
Resolver numericamente el sistema de Lorenz mediante un esquema de Runge-Kutta de cuarto orden (`RK4`) y analizar como cambia la dinamica cuando varia el parametro `r`, manteniendo `Pr = 10`, `b = 8/3`, `dt = 0.005` y la condicion inicial estandar `(x0, y0, z0) = (1, 1, 1)`.

El sistema integrado fue:

`dx/dt = Pr (y - x)`

`dy/dt = x (r - z) - y`

`dz/dt = xy - b z`

## Objetivos del inciso
- Implementar `RK4` de forma explicita y reproducible.
- Resolver el sistema para `r = 2`, `10`, `24`, `25` y `30` en el intervalo `0 <= t <= 50`.
- Comparar la evolucion temporal de `x(t)` y la trayectoria en el espacio de fases.
- Mostrar la sensibilidad a condiciones iniciales para `r = 30`.

## Entregables
- Notebook reproducible con la implementacion de `RK4` y las simulaciones pedidas.
- Documento de resultados y discusion para los distintos valores de `r`.
- Figuras de series temporales y espacio de fases exportadas para el book.

## Resolucion y notebook asociado
- [Resolucion Inciso 1: Sistema de Lorenz con RK4](resultados-discusion-lorenz-rk4.md)
- [Notebook asociado](../notebooks/Guia_practica_Lorenz_RK4.ipynb)

## Referencias base (opcional)
- Rutinas MATLAB de apoyo en `Tarea2ModAmb/dydt.m`, `Tarea2ModAmb/lorenz.m` y `Tarea2ModAmb/rkstep.m`.
- Lorenz, E. N. (1963). Deterministic nonperiodic flow. *Journal of the Atmospheric Sciences*, 20(2), 130-141.
