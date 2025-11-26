import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import imageio

nx, ny = 100, 100
archivo = "wave_evolution.dat"
matrices = []

with open(archivo) as f:
    for _ in range(5): next(f)
    matriz = np.zeros((ny, nx))
    last_step = 0
    for linea in f:
        partes = linea.strip().split()
        if not partes or len(partes) < 5:
            continue
        step, _, x, y, amp = map(float, partes)
        x, y = int(x), int(y)
        if step != last_step:
            matrices.append(matriz)
            matriz = np.zeros((ny, nx))
            last_step = step
        matriz[y, x] = amp
    matrices.append(matriz)

imagenes = []
X, Y = np.meshgrid(np.arange(nx), np.arange(ny))

for idx, matriz in enumerate(matrices[::10]):
    fig = plt.figure(figsize=(7,6))
    ax = fig.add_subplot(111, projection='3d')
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    Z = np.maximum(matriz, 0)    # Trunca negativos a cero
    surf = ax.plot_surface(X, Y, Z, cmap='viridis', linewidth=0, antialiased=False)
    ax.set_zlim(0, 1)            # Base en cero, amplitud máxima 1 (ajusta según tus datos)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Amplitud")
    ax.set_title(f"Step {idx*10}")
    ax.view_init(elev=40, azim=225)
    plt.tight_layout()
    plt.savefig(f'_tmp3d_{idx}.png', bbox_inches='tight', pad_inches=0, facecolor=fig.get_facecolor())
    plt.close(fig)
    imagenes.append(imageio.imread(f'_tmp3d_{idx}.png'))

imageio.mimsave('ondas_3d_base0.gif', imagenes, fps=30)
print("Animación 3D con fondo blanco y base en cero guardada como ondas_3d_base0.gif")
