import os
import numpy as np
import matplotlib.pyplot as plt

import lib.calculations as c
import lib.basis as basis

# Parameters
TIMESTEP = 2.0e-6
IDEAL_AXIS = np.array([1, 0, 0])  # ideal axis
ANGLE = 90  # ideal angle
OFFSRANGE = 10000  # in Hz
N_OFFSETS = 101  # offsrange/100 is good
U_DESIRED = np.array([[np.sqrt(2)/2, np.sqrt(2)*1j/2], [-np.sqrt(2)/2, np.sqrt(2)*1j/2]])  # for 90x pulse
OFFSETS = np.linspace(-OFFSRANGE/2, OFFSRANGE/2, N_OFFSETS)

def ensure_directories():
    """
    Ensure that necessary directories exist.
    """
    dirs = ['figures', 'figures/rotation_axis', 'figures/offset']
    for directory in dirs:
        if not os.path.exists(directory):
            os.makedirs(directory)

def get_pulse_filenames():
    """
    Get all filenames that start with 'pulse'.
    """
    return [f for f in os.listdir() if f.startswith('pulse')]

def plot_results(name, offsets):
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

    data = np.loadtxt(f'figures/offset/{name}_results.txt')

    fig, axes = plt.subplots(1, 3, figsize=(44, 7.8))
    plt.subplots_adjust(wspace=0.3)

    for i, ax in enumerate(axes):
        for j in range(3):
            ax.plot(offsets / 1000, data[:, i * 3 + j], line_styles[j], color=colors[j], linewidth=0.8, label=labels[i][j])
        ax.set_xlabel('Offset [kHz]', fontsize=18, fontweight='bold')
        ax.set_ylim([np.min(data[:, i * 3]) - 0.1, np.max(data[:, i * 3 + 2]) + 0.1])
        ax.set_xticks([-60, -40, -20, 0, 20, 40, 60])
        ax.legend(loc='center left', fontsize=18)
        if i == 0:
            ax.set_ylabel('Transfer efficiency', fontsize=18, fontweight='bold')

    plt.savefig(f'figures/offset/M_plot_{name}.pdf')

def main():
    """
    Main function to execute the analysis.
    """
    ensure_directories()
    pulse_files = get_pulse_filenames()

    for file_name in pulse_files:
        c.calculate_effective_propagator(file_name, OFFSRANGE, N_OFFSETS)
        c.calculate_rotation_axis(file_name, OFFSETS, N_OFFSETS)
        c.final_states(file_name, N_OFFSETS)
        plot_results(file_name, OFFSETS)
        c.calculate_quality_factor(file_name, N_OFFSETS, U_DESIRED)

if __name__ == "__main__":
    main()