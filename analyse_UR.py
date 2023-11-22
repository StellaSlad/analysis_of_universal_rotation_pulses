# Analyse UR
# 
# This python script aims to provide a thorough analysis of broadband UR pulses.
# More details..

import os
import numpy as np
import scipy.linalg as scla
import pandas as pd
import matplotlib.pyplot as plt
import calculations as c

#_____________________________________________
# Parameters
# timestep = 2.0e-6
ideal_axis = np.array([1, 0, 0])  # ideal axis
angle = 90  # ideal angle
offsrange = 10000  # in Hz
n_offsets = 101  # offsrange/100 is good
# Desired target propagator
# U_F = np.array([[1, 0], [0, 1j]])  # for 180x pulses
U_desired = np.array([[np.sqrt(2)/2, np.sqrt(2)*1j/2], [-np.sqrt(2)/2, np.sqrt(2)*1j/2]])  # for 90x pulse
offsets = np.linspace(-offsrange/2, offsrange/2, n_offsets)

#________________
# Functions

def plot_results(name, offsets):
    # Parameters - should be commented out
    # offsrange = 120000  # in Hz
    # n_offsets = 501  # offsrange/100 is good
    # offsets = np.linspace(-offsrange/2, offsrange/2, n_offsets)

    # Preparation
    l = ['-', '--', '-.']  # line style
    li = [str(x) + 'k' for x in l]
    M = np.loadtxt(f'figures/offset/{name}_results.txt')

    # Plotting
    fig = plt.figure(figsize=(44, 7.8))
    plt.subplots_adjust(wspace=0.3)  # Adjust the width between subplots

    ax1 = plt.subplot(1, 3, 1)
    ax1.plot(offsets / 1000, M[:, 0], li[0], color=[0.32, 0.9, 0.8], linewidth=0.8, label='x→x')
    ax1.plot(offsets / 1000, M[:, 1], li[1], color=[0.1, 0.6, 0.52], linewidth=0.8, label='x→y')
    ax1.plot(offsets / 1000, M[:, 2], li[2], color='k', linewidth=0.8, label='x→z')
    ax1.set_xlabel('Offset [kHz]', fontsize=18, fontweight='bold')
    ax1.set_ylabel('Transfer efficiency', fontsize=18, fontweight='bold')
    ax1.legend(loc='center left', fontsize=18)
    ax1.set_ylim([np.min(M[:, 0]) - 0.1, np.max(M[:, 2]) + 0.1])
    ax1.set_xticks([-60, -40, -20, 0, 20, 40, 60])

    ax2 = plt.subplot(1, 3, 2)
    ax2.plot(offsets / 1000, M[:, 3], li[0], color=[0.32, 0.9, 0.8], linewidth=0.8, label='y→x')
    ax2.plot(offsets / 1000, M[:, 4], li[1], color=[0.1, 0.6, 0.52], linewidth=0.8, label='y→y')
    ax2.plot(offsets / 1000, M[:, 5], li[2], color='k', linewidth=0.8, label='y→z')
    ax2.set_xlabel('Offset [kHz]', fontsize=18, fontweight='bold')
    ax2.set_ylim([np.min(M[:, 3]) - 0.1, np.max(M[:, 5]) + 0.1])
    ax2.set_xticks([-60, -40, -20, 0, 20, 40, 60])

    ax3 = plt.subplot(1, 3, 3)
    ax3.plot(offsets / 1000, M[:, 6], li[0], color=[0.32, 0.9, 0.8], linewidth=0.8, label='z→x')
    ax3.plot(offsets / 1000, M[:, 7], li[1], color=[0.1, 0.6, 0.52], linewidth=0.8, label='z→y')
    ax3.plot(offsets / 1000, M[:, 8], li[2], color='k', linewidth=0.8, label='z→z')
    ax3.set_xlabel('Offset [kHz]', fontsize=18, fontweight='bold')
    ax3.set_ylim([np.min(M[:, 6]) - 0.1, np.max(M[:, 8]) + 0.1])
    ax3.set_xticks([-60, -40, -20, 0, 20, 40, 60])

    plt.savefig(f'figures/offset/M_plot_{name}.pdf')

#_____________________________________________
# Main

# Create necessary subdirectories
if not os.path.exists('figures'):
    os.mkdir('figures')
if not os.path.exists('figures/rotation_axis'):
    os.mkdir('figures/rotation_axis')
if not os.path.exists('figures/offset'):
    os.mkdir('figures/offset')

# Get all file names
S = [f for f in os.listdir() if f.startswith('pulse')]

# Execute different functions to analyze the pulse
for file_name in S:
    c.calculate_effective_propagator(file_name, offsrange, n_offsets)
    c.calculate_rotation_axis(file_name, offsets, n_offsets)
    c.final_states(file_name, n_offsets)
    plot_results(file_name, offsets)
    c.calculate_quality_factor(file_name, n_offsets, U_desired)