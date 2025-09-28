# Usar imagen base oficial de Ubuntu
FROM ubuntu:22.04

# Instalar compiladores gcc y g++, make y librería OpenMP
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    libomp-dev \
    make \
    && rm -rf /var/lib/apt/lists/*

# Establecer directorio de trabajo
WORKDIR /app

# Copiar todos los archivos de código fuente
COPY . .

# Compilar el código con g++ y OpenMP (Si tienes varios archivos .cpp)
RUN g++ main.cpp Network.cpp Node.cpp WavePropagator.cpp -o main -fopenmp

# Comando para ejecutar al iniciar el contenedor
CMD ["./main"]