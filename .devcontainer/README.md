# Para ejecucion se tiene los siguientes comandos:

# Compilar
make

# Ejecutar benchmarks
make benchmark

# Ejecutar an análisis
make analysis

# Limpiar
make clean


Por defecto se tienen lo valores:

Número de nodos: 10000
Pasos de simulación: 1000
Paso de tiempo: 0.01
Tamaño de red: 100x100 (para redes 2D)
Coeficiente de difusión: D = 0,1
Coeficiente de amortiguación: γ = 0,01
Amplitud inicial: A0 = 1,0
Fuente externa: S = 0,1 sin(ωt)

Para probar en una dimension se deja un ejemplo de la prueba a realizar:


--------


Esto debe ser modificado en el archivo "main.cpp".