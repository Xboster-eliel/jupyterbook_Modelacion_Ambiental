---
title: "Inciso 3: Modelos lineales vs procesos"
---

# Inciso 3: Modelos lineales vs procesos

## Planteamiento
Modelar y comparar modelos estadisticos y modelos basados en procesos para predecir el caudal en cuencas colombianas de tamano mediano a pequeno. Para la modelacion se usa el `75%` inicial de cada serie y el `25%` final para evaluacion.

En este inciso se resolvio el ejercicio para dos cuencas:

- `APARTADO [12017060]`
- `PUENTE CARRETERA - AUT [21137030]`

En ambos casos se implementaron:

1. un modelo lineal `ARX(3,1)` como referencia estadistica principal,
2. una familia `ARIMAX`,
3. una familia `SARIMAX`,
4. un modelo fisico parsimonioso basado en el balance de agua de largo plazo.

## Preguntas guia
- Que tan competitivos son los modelos lineales frente a un modelo fisico simple de balance hidrologico.
- Como cambia el desempeno predictivo entre cuencas con comportamiento hidrologico distinto.
- Que se gana en interpretabilidad hidrologica y que se pierde en habilidad predictiva al usar un modelo basado en procesos muy parsimonioso.

## Entregables
- Notebook reproducible con el procesamiento y analisis de ambas cuencas.
- Documento de resultados y discusion consolidado.
- Comparacion transversal entre cuencas y entre familias de modelos.

## Resolucion y notebook asociado
- [Resolucion Inciso 3: Modelos lineales vs procesos](resultados-discusion-modelos-linealesvsprocesos)
- [Notebook asociado](guia-practica-modelos-linealesvsprocesos)

## Referencias base (opcional)
- NASA POWER API mensual por punto para temperatura del aire a 2 m.
- CHIRPS para precipitacion diaria.
- IDEAM para caudal medio diario.
- Thornthwaite para PET mensual.
