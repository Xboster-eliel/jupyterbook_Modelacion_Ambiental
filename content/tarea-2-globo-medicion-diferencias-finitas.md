---
title: "Inciso 1: Globo de medicion de variables meteorologicas"
---

# Inciso 1: Globo de medicion de variables meteorologicas

## Planteamiento
Un globo de medicion de variables meteorologicas asciende verticalmente hasta alcanzar su nivel de boyancia neutra. Cerca de ese punto, el movimiento se modela como una oscilacion suave alrededor de la posicion de equilibrio.

La ecuacion diferencial ordinaria que representa este proceso es:

$$
m\,\frac{d^2 y}{dt^2} = -\beta\,\frac{dy}{dt} - \gamma\,y,
$$

donde `y(t)` es el desplazamiento vertical respecto al nivel de equilibrio, `y = 0` representa la posicion de equilibrio, `m` es la masa del globo, `\beta` representa la resistencia del aire y `\gamma` es la constante lineal de boyancia restauradora. En forma equivalente, tambien puede escribirse como `m y'' + \beta y' + \gamma y = 0`.

Las condiciones iniciales del inciso son:

$$
y(0) = 0, \qquad y'(0) = v_0.
$$

En esta tarea revisaremos la solucion numerica de la ecuacion y, mas importante aun, el analisis critico de los resultados que entrega la discretizacion.

## Actividades
### a. Derivacion del esquema
Usando diferencias centrales para la segunda derivada y diferencia hacia adelante para la primera derivada en el punto `t_n`, derive el esquema de diferencias finitas para obtener `y_(n+1)`.

Considere la aproximacion centrada de segundo orden:

$$
y''(t_n) \approx \frac{y_{n+1} - 2y_n + y_{n-1}}{h^2}.
$$

### b. Comparacion con Euler hacia adelante
Explique en que se diferencia este esquema respecto del metodo de Euler hacia adelante. En particular, discuta que cambia en la forma de usar las condiciones iniciales y por que este detalle modifica la implementacion numerica.

### c. Inicio de la integracion y exploracion numerica
Explique como se inicia la integracion numerica. Plantee su estrategia para calcular el primer paso y luego integre tantos casos como considere utiles para explorar como esta discretizacion resuelve el problema.

## Objetivos del inciso
- Derivar de forma consistente un esquema explicito de dos pasos para una EDO amortiguada de segundo orden.
- Comparar la estructura del metodo con Euler hacia adelante y discutir el papel de las condiciones iniciales.
- Evaluar criticamente el efecto del paso temporal y de los parametros `m`, `\beta` y `\gamma` sobre la solucion numerica.

## Entregables
- Derivacion algebraica clara del esquema de diferencias finitas.
- Discusion comparativa entre el esquema obtenido y el metodo de Euler hacia adelante.
- Analisis critico de resultados numericos para distintos casos de integracion.

## Resolucion asociada
- [Resolucion Inciso 1: Globo de medicion de variables meteorologicas](resultados-discusion-globo-medicion.md)
- [Notebook asociado](../notebooks/Guia_practica_Globo_medicion.ipynb)
