import os
import numpy as np
import matplotlib.pyplot as plt

import lib.calculations as c
import lib.plot_functions as pf

# Parameters
DIRECTORY = 'pulses/'
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
    dirs = ['txt', 'images']
    for directory in dirs:
        if not os.path.exists(directory):
            os.makedirs(directory)

def get_pulse_filenames():
    """
    Get all filenames that start with 'pulse'.
    """
    return [f for f in os.listdir() if f.startswith('pulse')]

def main():
    """
    Main function to execute the analysis.
    """
    ensure_directories()
    os.chdir(DIRECTORY)
    pulse_files = get_pulse_filenames()

    for file_name in pulse_files:
        c.calculate_effective_propagator(file_name, OFFSRANGE, N_OFFSETS)
        Lxx, Lyy, Lzz, rx, ry, rz = c.calculate_rotation_axis(file_name, OFFSETS, N_OFFSETS)
        pf.plot_rotation_axis(file_name, OFFSETS, rx, ry, rz, Lxx, Lyy, Lzz)
        c.calculate_final_states(file_name, N_OFFSETS)
        pf.plot_transfer_efficiency(file_name, OFFSETS)
        c.calculate_quality_factor(file_name, N_OFFSETS, U_DESIRED)

if __name__ == "__main__":
    main()