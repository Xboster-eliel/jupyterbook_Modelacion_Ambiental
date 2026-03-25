---
title: "Resolucion Inciso 1: Globo de medicion de variables meteorologicas"
---

# Resolucion Inciso 1: Globo de medicion de variables meteorologicas

## Modelo continuo
El movimiento vertical del globo cerca de su nivel de equilibrio se representa mediante la ecuacion:

$$
m\,\ddot y + b\,\dot y + k\,y = 0,
$$

donde `m` es la masa, `b` el coeficiente de amortiguamiento por resistencia del aire, `k` la constante lineal de boyancia y `y(t)` el desplazamiento vertical respecto al equilibrio. Las condiciones iniciales son:

$$
y(0) = y_0, \qquad \dot y(0) = v_0.
$$

## Derivacion del esquema de diferencias
Se discretiza el tiempo con un paso uniforme `h = Delta t` y se aproxima la segunda derivada con diferencias centrales:

$$
\ddot y(t_n) \approx \frac{y_{n+1} - 2y_n + y_{n-1}}{h^2},
$$

mientras que la primera derivada se aproxima con diferencia hacia adelante:

$$
\dot y(t_n) \approx \frac{y_{n+1} - y_n}{h}.
$$

Al sustituir ambas expresiones en la ecuacion diferencial, se obtiene:

$$
m\frac{y_{n+1} - 2y_n + y_{n-1}}{h^2} + b\frac{y_{n+1} - y_n}{h} + ky_n = 0.
$$

Multiplicando por `h^2` y reordenando terminos:

$$
(m + bh) y_{n+1} = (2m + bh - kh^2) y_n - m y_{n-1}.
$$

Por tanto, el esquema explicito para avanzar la solucion es:

$$
y_{n+1} = \frac{2m + bh - kh^2}{m + bh} y_n - \frac{m}{m + bh} y_{n-1}.
$$

Este resultado muestra que el metodo depende de dos estados previos, `y_n` y `y_(n-1)`, por lo que se trata de un esquema de dos pasos.

## Diferencia con Euler hacia adelante
La diferencia central con Euler hacia adelante no es solo algebraica, sino tambien estructural.

- Euler hacia adelante es un metodo de un paso: para avanzar solo necesita el valor actual de la solucion y la pendiente evaluada en ese mismo instante.
- En una ecuacion de segundo orden, Euler suele aplicarse despues de reescribir el problema como un sistema de primer orden para `y` y `v = dy/dt`.
- El esquema derivado aqui no necesita introducir explicitamente una variable de velocidad en cada paso, pero si exige disponer de dos posiciones consecutivas.
- Por esa razon, las condiciones iniciales no bastan por si solas para iterar directamente la formula. Ademas de `y_0`, hace falta construir un valor inicial adicional `y_1`.

Desde el punto de vista numerico, este detalle cambia la implementacion: Euler arranca de forma inmediata con `(y_0, v_0)`, mientras que el esquema de diferencias requiere una etapa de arranque antes de usar la recurrencia principal.

## Inicio de la integracion numerica
La integracion se inicia con `y_0 = y(0)` y con una aproximacion para `y_1`. La opcion mas simple consiste en usar Euler hacia adelante sobre la posicion:

$$
y_1 = y_0 + v_0 h.
$$

Sin embargo, una aproximacion mas consistente con la dinamica del problema usa Taylor de segundo orden:

$$
y_1 = y_0 + v_0 h + \frac{1}{2} a_0 h^2,
$$

donde la aceleracion inicial se obtiene de la ecuacion original:

$$
a_0 = \ddot y(0) = -\frac{b v_0 + k y_0}{m}.
$$

Con `y_0` y `y_1` ya disponibles, la solucion se genera de manera recursiva para `n = 1, 2, 3, ...` usando el esquema derivado.

## Analisis critico de la discretizacion
El interes principal de este ejercicio no es solo obtener una trayectoria numerica, sino evaluar si la discretizacion representa correctamente la fisica del sistema.

### Precision esperada
La aproximacion de la segunda derivada es de orden `O(h^2)`, mientras que la derivada primera se aproxima con una formula hacia adelante de orden `O(h)`. En consecuencia, el esquema completo no hereda automaticamente un orden alto en todos sus terminos, y la eleccion de `h` sigue siendo decisiva para controlar el error.

### Papel del paso temporal
Si `h` es demasiado grande, la solucion puede distorsionar la frecuencia de oscilacion, exagerar o atenuar artificialmente el amortiguamiento e incluso producir inestabilidades numericas. Si `h` es suficientemente pequeno, el metodo captura mejor la oscilacion amortiguada y el decaimiento progresivo de la amplitud.

### Interpretacion fisica
Cuando `b > 0`, la amplitud debe disminuir con el tiempo porque el sistema pierde energia por resistencia del aire. Una solucion numerica razonable debe reproducir ese decaimiento sin introducir crecimiento espurio. Asimismo, la rapidez de la oscilacion debe estar ligada a la escala `sqrt(k/m)`, de modo que cambios en `k` o `m` modifiquen la frecuencia observada de manera fisicamente coherente.

### Comparacion conceptual con Euler
Para un mismo `h`, este esquema suele representar mejor la curvatura temporal de la posicion que un Euler directo sobre el sistema equivalente de primer orden. Sin embargo, esa ventaja potencial exige un arranque bien construido y una revision critica de la estabilidad. En otras palabras, un metodo aparentemente mas sofisticado no garantiza por si solo mejores resultados si se usa con un paso inadecuado o con una inicializacion pobre.

## Conclusiones

1. El esquema de diferencias finitas derivado para el globo es explicito y de dos pasos, porque calcula `y_(n+1)` a partir de `y_n` y `y_(n-1)`.
2. La diferencia principal frente a Euler hacia adelante esta en la forma de iniciar la integracion: aqui no basta con `y_0` y `v_0`, sino que se debe estimar tambien `y_1`.
3. La calidad de la solucion depende de manera fuerte del paso temporal `h`, tanto por precision como por estabilidad.
4. El analisis critico debe verificar si la amplitud decrece de forma fisicamente razonable y si la frecuencia numerica concuerda con la dinamica esperada del sistema amortiguado.
5. El mejor uso de esta discretizacion no consiste solo en ejecutar la recurrencia, sino en comparar casos y discutir cuando el metodo representa bien o mal el movimiento del globo.

## Referencias base
- Ortega, R. (2014). *Diferenciacion numerica: aplicaciones computacionales*. Universidad de Santiago de Chile. https://mecanica-usach.mine.nu/media/uploads/L06_DiferenciacionNumerica.pdf
- Cordero, P., y Soto, R. (2011). *Ecuaciones diferenciales ordinarias: integracion numerica*. Universidad de A Coruna. http://caminos.udc.es/info/asignaturas/grado_itop/221/images/Imagenes_complementarios/Edos_teoria.pdf
