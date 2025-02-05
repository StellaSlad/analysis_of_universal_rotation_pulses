import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_transfer_efficiency(name: str, offsets: np.ndarray):
    """
    Plot the transfer efficiency results of the analysis.

    Parameters:
    name (str): The name of the file.
    offsets (ndarray): Array of offset values.
    """
    line_styles = ['-', '--', '-.']
    colors = [
        [0.32, 0.9, 0.8],
        [0.1, 0.6, 0.52],
        'k'
    ]
    labels = [
        ['x→x', 'x→y', 'x→z'],
        ['y→x', 'y→y', 'y→z'],
        ['z→x', 'z→y', 'z→z']
    ]

    M = np.loadtxt(f'figures/offset/{name}_results.txt')

    fig, axs = plt.subplots(1, 3, figsize=(17.32, 6.14))  # 44 cm x 7.8 cm

    for i in range(3):
        axs[i].plot(offsets / 1000, M[:, 3 * i], line_styles[0], color=colors[0], linewidth=0.8)
        axs[i].plot(offsets / 1000, M[:, 3 * i + 1], line_styles[1], color=colors[1], linewidth=0.8)
        axs[i].plot(offsets / 1000, M[:, 3 * i + 2], line_styles[2], color=colors[2], linewidth=0.8)
        axs[i].set_xlabel('Offset [kHz]', fontsize=18, fontweight='bold')
        axs[i].set_ylabel('Transfer efficiency', fontsize=18, fontweight='bold')
        axs[i].legend(labels[i], loc='eastoutside', fontsize=18)
        miny = [min(M[:, 3 * i]), min(M[:, 3 * i + 1]), min(M[:, 3 * i + 2])]
        maxy = [max(M[:, 3 * i]), max(M[:, 3 * i + 1]), max(M[:, 3 * i + 2])]
        axs[i].set_ylim([min(miny) - 0.1, max(maxy) + 0.1])
        axs[i].tick_params(axis='both', which='major', labelsize=14, direction='out')

    plt.tight_layout()
    plt.savefig(f'figures/offset/te_plot_{name}.jpg')
    plt.savefig(f'figures/offset/te_plot_{name}.pdf')
    plt.savefig(f'figures/offset/te_plot_{name}.png')
    plt.show()

def plot_rotation_axis(name: str, offsets: np.ndarray, rx: list, ry: list, rz: list, Lxx: list, Lyy: list, Lzz: list):
    """
    Plot the rotation axis and save the plot to a file.

    Parameters:
    name (str): The name of the file.
    offsets (np.ndarray): Array of offset values.
    rx (list): Rotation vector x components.
    ry (list): Rotation vector y components.
    rz (list): Rotation vector z components.
    Lxx (list): Directional cosines x components.
    Lyy (list): Directional cosines y components.
    Lzz (list): Directional cosines z components.

    Returns:
    None
    """
    line_styles = ['-', '--', '-.']
    colors = [
        [0.32, 0.9, 0.8],
        [0.1, 0.6, 0.52],
        'k'
    ]
    labels = ['r_x', 'r_y', 'r_z']

    fig = plt.figure(figsize=(10.5, 16))

    ax1 = fig.add_subplot(211)
    ax1.plot(offsets / 1000, rx, line_styles[0], color=colors[0], linewidth=1.2, label=labels[0])
    ax1.plot(offsets / 1000, ry, line_styles[1], color=colors[1], linewidth=1.2, label=labels[1])
    ax1.plot(offsets / 1000, rz, line_styles[2], color=colors[2], linewidth=1.2, label=labels[2])
    ax1.set_xlabel('Offset [kHz]', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Rotation vector', fontsize=14, fontweight='bold')
    ax1.legend(loc='center left', fontsize=14)

    ax2 = fig.add_subplot(212, projection='3d')
    for i in range(len(Lxx)):
        ax2.quiver(0, 0, 0, Lxx[i], Lyy[i], Lzz[i], color=[0.12, 0.7, 0.6], linewidth=0.5, arrow_length_ratio=0.1)

    ax2.set_xlabel('x', fontsize=14, fontweight='bold')
    ax2.set_ylabel('y', fontsize=14, fontweight='bold')
    ax2.set_zlabel('z', fontsize=14, fontweight='bold')
    ax2.set_xlim([-1, 1])
    ax2.set_ylim([-1, 1])
    ax2.set_zlim([-1, 1])
    ax2.grid(True, which='minor')

    plt.savefig(f'figures/rotation_axis/rot_axis_{name}.pdf')
    plt.show()

def plot_transfer_efficiency2(name, offsets):
    """
    Plot the results of the analysis.

    Parameters:
    name (str): The name of the file.
    offsets (ndarray): Array of offset values.
    """
    line_styles = ['-', '--', '-.']
    colors = [
        [0.32, 0.9, 0.8],
        [0.1, 0.6, 0.52],
        'k'
    ]
    labels = [
        ['x→x', 'x→y', 'x→z'],
        ['y→x', 'y→y', 'y→z'],
        ['z→x', 'z→y', 'z→z']
    ]

    M = np.loadtxt(f'figures/offset/{name}_results.txt')

    fig, axs = plt.subplots(1, 3, figsize=(17.32, 6.14))  # 44 cm x 7.8 cm

    for i in range(3):
        axs[i].plot(offsets / 1000, M[:, 3 * i], line_styles[0], color=colors[0], linewidth=0.8)
        axs[i].plot(offsets / 1000, M[:, 3 * i + 1], line_styles[1], color=colors[1], linewidth=0.8)
        axs[i].plot(offsets / 1000, M[:, 3 * i + 2], line_styles[2], color=colors[2], linewidth=0.8)
        axs[i].set_xlabel('Offset [kHz]', fontsize=18, fontweight='bold')
        axs[i].set_ylabel('Transfer efficiency', fontsize=18, fontweight='bold')
        axs[i].legend(labels[i], loc='right', fontsize=18)
        miny = [min(M[:, 3 * i]), min(M[:, 3 * i + 1]), min(M[:, 3 * i + 2])]
        maxy = [max(M[:, 3 * i]), max(M[:, 3 * i + 1]), max(M[:, 3 * i + 2])]
        axs[i].set_ylim([min(miny) - 0.1, max(maxy) + 0.1])
        axs[i].tick_params(axis='both', which='major', labelsize=14, direction='out')

    plt.tight_layout()
    plt.savefig(f'figures/offset/M_plot_{name}.jpg')
    plt.savefig(f'figures/offset/M_plot_{name}.pdf')
    plt.savefig(f'figures/offset/M_plot_{name}.png')
    plt.show()