import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import imageio
import os

archivo = "wave_evolution.dat"     

# TIPO DE ENCABEZADO ESPERADO 1D:
# Simulacion 1D 10000 pasos=1000 dt=0.01
# D=0.1 gamma=0.01 amp0=1
# Fuente: 0.1 sin(6.28319t)
# Formato: step time x amp

# TIPO DE ENCABEZADO ESPERADO 2D:
# Simulacion 2D 100x100 pasos=1000 dt=0.01
# D=0.1 gamma=0.01 amp0=1
# Fuente: 0.1 sin(6.28319t)
# Formato: step time x y amp

# Detección y lectura del encabezado
with open(archivo) as f:
    encabezado = [next(f).lower().strip() for _ in range(5)]
    formato_line = [line for line in encabezado if "formato:" in line][0]

    is_2d = ("x y amp" in formato_line)
    is_1d = ("x amp" in formato_line)
    # Puedes agregar lógica para otros formatos si lo necesitas

matrices = []

with open(archivo) as f:
    for _ in range(5): next(f)  # Saltar siempre 5 líneas de encabezado
    if is_2d:
        nx, ny = 100, 100       # Ajusta según tus parámetros reales (puedes extraer del encabezado si lo deseas)
        matriz = np.zeros((ny, nx))
        last_step = 0
        for linea in f:
            partes = linea.strip().split()
            if not partes or len(partes) < 5: continue
            step, _, x, y, amp = map(float, partes)
            x, y = int(x), int(y)
            if step != last_step:
                matrices.append(matriz)
                matriz = np.zeros((ny, nx))
                last_step = step
            matriz[y, x] = amp
        matrices.append(matriz)

    elif is_1d:
        nx = 10000              # Ajusta aquí si tu tamaño cambia
        amplitudes = np.zeros(nx)
        last_step = 0
        for linea in f:
            partes = linea.strip().split()
            if not partes or len(partes) < 4: continue
            step, _, x, amp = map(float, partes)
            x = int(x)
            if step != last_step:
                matrices.append(amplitudes.copy())
                amplitudes[:] = 0.0
                last_step = step
            amplitudes[x] = amp
        matrices.append(amplitudes.copy())

imagenes = []

# Generación de animaciones
if is_2d:
    X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
    for idx, matriz in enumerate(matrices[::10]):
        fig = plt.figure(figsize=(7,6))
        ax = fig.add_subplot(111, projection='3d')
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')
        Z = np.maximum(matriz, 0)
        surf = ax.plot_surface(X, Y, Z, cmap='viridis', linewidth=0, antialiased=False)
        ax.set_zlim(0, 1)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Amplitud")
        ax.set_title(f"Step {idx*10}")
        ax.view_init(elev=40, azim=225)
        plt.tight_layout()
        plt.savefig(f'_tmp3d_{idx}.png', bbox_inches='tight', pad_inches=0, facecolor=fig.get_facecolor())
        plt.close(fig)
        imagenes.append(imageio.imread(f'_tmp3d_{idx}.png'))
        os.remove(f'_tmp3d_{idx}.png')
    imageio.mimsave('ondas_3d_base0.gif', imagenes, fps=30)
    print("Animación 2D guardada como ondas_3d_base0.gif")

elif is_1d:
    X = np.arange(nx)
    for idx, amplitudes in enumerate(matrices[::10]):
        fig = plt.figure(figsize=(9, 5))
        plt.plot(X, amplitudes, color='royalblue')
        plt.xlim(0, nx-1)
        plt.ylim(0, 1)
        plt.xlabel("Posición")
        plt.ylabel("Amplitud")
        plt.title(f"Step {idx*10}")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f'_tmp1d_{idx}.png', bbox_inches="tight", pad_inches=0, facecolor='white')
        plt.close(fig)
        imagenes.append(imageio.imread(f'_tmp1d_{idx}.png'))
        os.remove(f'_tmp1d_{idx}.png')
    imageio.mimsave('onda_1d_evolucion.gif', imagenes, fps=30)
    print("Animación 1D guardada como onda_1d_evolucion.gif")

else:
    print("Formato no reconocido en el encabezado del archivo.")
