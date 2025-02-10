import numpy as np
import pandas as pd
import scipy.linalg as scla
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from lib.basis import basis

def read_pulse_file(name: str) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Read the shaped pulse file.

    Parameters:
    name (str): The name of the file.

    Returns:
    tuple: x_amp, y_amp, timestep
    """
    pulse = pd.read_table(name, skiprows=4, header=None, names=['Var1', 'Var2', 'Var3'], sep=r'\s+')
    x_amp = pulse['Var1'].astype(float).values
    y_amp = pulse['Var2'].astype(float).values
    timestep = float(pulse['Var3'][0])
    return x_amp, y_amp, timestep

def calculate_time_propagator(timestep: float, Hevo: np.ndarray) -> np.ndarray:
    """
    Calculate the time propagator.

    Parameters:
    timestep (float): The time step.
    Hevo (np.ndarray): The Hamiltonian.

    Returns:
    np.ndarray: The time propagator.
    """
    return scla.expm(-1j * timestep * Hevo)

def calculate_effective_propagator(name: str, offsrange: float, n_offsets: int):
    """
    Calculate the effective propagator for a given shaped pulse file over a range of offsets.

    This function performs the following steps:
    1. Initializes the spin system.
    2. Reads the shaped pulse file to obtain the pulse amplitudes and timestep.
    3. Creates an offset grid based on the specified range and number of offsets.
    4. Initializes the initial state and preallocates space for the results.
    5. Loops over the offset grid to calculate the effective propagator at each offset.
    6. For each offset, loops over the time grid to propagate the state using the time propagator.
    7. Concatenates the results and saves them to a file.

    Parameters:
    name (str): The name of the file.
    offsrange (float): The offset range.
    n_offsets (int): The number of offsets.

    Returns:
    None
    """
    nspins = 1
    _, ix, iy, iz, _, _, _, _, _ = basis(nspins)

    x_amp, y_amp, timestep = read_pulse_file(name)
    offsets = np.linspace(-offsrange / 2, offsrange / 2, n_offsets)
    rhoinit = iz[:, :, 0]
    nsteps = len(x_amp)
    U_eff = np.empty((2, 0), dtype=np.complex128)  # Initialize as empty array

    for n in range(n_offsets):
        rho = rhoinit
        H_chemical_shift = 2 * np.pi * offsets[n] * iz[:, :, 0]
        U_eff_dummy = np.eye(2)

        for k in range(nsteps):
            H_rf = 2 * np.pi * (x_amp[k] * ix[:, :, 0] + y_amp[k] * iy[:, :, 0])
            H_total = H_chemical_shift + H_rf
            U_total = calculate_time_propagator(timestep, H_total)
            rho = np.dot(U_total, np.dot(rho, np.conj(U_total).T))
            U_eff_dummy = np.dot(U_eff_dummy, U_total)

        U_eff = np.concatenate((U_eff, U_eff_dummy), axis=1)
        # print(f"U_eff shape after offset {n}: {U_eff.shape}")  # Debug statement to check the shape

    # Print the final dimensions of U_eff
    print(f"Final dimensions of U_eff: {U_eff.shape}")
    np.savetxt(f'../txt/{name}_u.txt', U_eff, delimiter='\t', fmt='%.6f')

def calculate_transfer_function(rho: np.ndarray, operator: np.ndarray) -> float:
    """
    Calculate the transfer function component.

    Parameters:
    rho (np.ndarray): The density matrix.
    operator (np.ndarray): The operator matrix.

    Returns:
    float: The transfer function component.
    """
    return np.real(np.trace(np.dot(np.conj(rho).T, operator)) / np.trace(np.dot(operator, operator)))

def calculate_final_states(name: str, n_offsets: int):
    """
    Calculate the final states for a given shaped pulse file over a range of offsets.

    Parameters:
    name (str): The name of the file.
    n_offsets (int): The number of offsets.

    Returns:
    None
    """
    nspins = 1
    _, ix, iy, iz, _, _, _, _, _ = basis(nspins)

    U = np.loadtxt(f'../txt/{name}_u.txt', dtype=np.complex128)
    
    M_xx, M_xy, M_xz = np.zeros(n_offsets), np.zeros(n_offsets), np.zeros(n_offsets)
    M_yx, M_yy, M_yz = np.zeros(n_offsets), np.zeros(n_offsets), np.zeros(n_offsets)
    M_zx, M_zy, M_zz = np.zeros(n_offsets), np.zeros(n_offsets), np.zeros(n_offsets)

    for i in range(n_offsets):
        U_eff = U[0:2, 2 * i:2 * i + 2]

        rho_x = np.dot(np.dot(U_eff, ix[:, :, 0]), np.conj(U_eff).T)
        rho_y = np.dot(np.dot(U_eff, iy[:, :, 0]), np.conj(U_eff).T)
        rho_z = np.dot(np.dot(U_eff, iz[:, :, 0]), np.conj(U_eff).T)

        M_xx[i] = calculate_transfer_function(rho_x, ix[:, :, 0])
        M_xy[i] = calculate_transfer_function(rho_x, iy[:, :, 0])
        M_xz[i] = calculate_transfer_function(rho_x, iz[:, :, 0])

        M_yx[i] = calculate_transfer_function(rho_y, ix[:, :, 0])
        M_yy[i] = calculate_transfer_function(rho_y, iy[:, :, 0])
        M_yz[i] = calculate_transfer_function(rho_y, iz[:, :, 0])

        M_zx[i] = calculate_transfer_function(rho_z, ix[:, :, 0])
        M_zy[i] = calculate_transfer_function(rho_z, iy[:, :, 0])
        M_zz[i] = calculate_transfer_function(rho_z, iz[:, :, 0])

    M = np.column_stack((M_xx, M_xy, M_xz, M_yx, M_yy, M_yz, M_zx, M_zy, M_zz))
    np.savetxt(f'../txt/{name}_results.txt', M, delimiter='\t', fmt='%.6f')

def calculate_quality_factor(name: str, n_offsets: int, desired_propagator: np.ndarray):
    """
    Calculate the quality factor for a given shaped pulse file over a range of offsets.

    Parameters:
    name (str): The name of the file.
    n_offsets (int): The number of offsets.
    desired_propagator (np.ndarray): The desired propagator.

    Returns:
    None
    """
    U = np.loadtxt(f'../txt/{name}_u.txt', dtype=np.complex128)

    transfer_efficiency = np.zeros(n_offsets + 1)
    for i in range(n_offsets):
        effective_propagator = U[0:2, 2 * i:2 * i + 2]
        transfer_efficiency[i] = np.real(np.trace(np.conj(desired_propagator).T @ effective_propagator))
    
    transfer_efficiency[n_offsets] = np.sum(transfer_efficiency) / n_offsets
    np.savetxt(f'../txt/quality_factor_{name}.txt', transfer_efficiency, delimiter='\t', fmt='%.6f')

def calculate_rotation_axis(name: str, offsets: np.ndarray, n_offsets: int):
    """
    Calculate the rotation axis for a given shaped pulse file over a range of offsets.

    Parameters:
    name (str): The name of the file.
    offsets (np.ndarray): Array of offset values.
    n_offsets (int): The number of offsets.

    Returns:
    tuple: Lists of Lxx, Lyy, Lzz, rx, ry, rz.
    """
    U = np.loadtxt(f'../txt/{name}_u.txt', dtype=np.complex128)

    a = U[0, 0::2]
    b = U[0, 1::2]

    A = np.imag(b)
    B = np.real(b)
    C = np.imag(a)
    D = np.real(a)

    omega_half = np.arccos(D)
    S = np.sin(omega_half)

    Lxx_new = A / S
    Lyy_new = B / S
    Lzz_new = C / S

    rx_new = 2 * omega_half * Lxx_new
    ry_new = 2 * omega_half * Lyy_new
    rz_new = 2 * omega_half * Lzz_new

    Lxx, Lyy, Lzz, rx, ry, rz = Lxx_new.tolist(), Lyy_new.tolist(), Lzz_new.tolist(), rx_new.tolist(), ry_new.tolist(), rz_new.tolist()
    axis = np.column_stack((Lxx, Lyy, Lzz))
    np.savetxt(f'../txt/{name}_axis.txt', axis, delimiter='\t', fmt='%.6f')

    return Lxx, Lyy, Lzz, rx, ry, rz


# def plot_rotation_axis(name: str, offsets: np.ndarray, rx: list, ry: list, rz: list, Lxx: list, Lyy: list, Lzz: list):
#     """
#     Plot the rotation axis and save the plot to a file.

#     Parameters:
#     name (str): The name of the file.
#     offsets (np.ndarray): Array of offset values.
#     rx (list): Rotation vector x components.
#     ry (list): Rotation vector y components.
#     rz (list): Rotation vector z components.
#     Lxx (list): Directional cosines x components.
#     Lyy (list): Directional cosines y components.
#     Lzz (list): Directional cosines z components.

#     Returns:
#     None
#     """
#     line_styles = ['-', '--', '-.']
#     colors = [
#         [0.32, 0.9, 0.8],
#         [0.1, 0.6, 0.52],
#         'k'
#     ]
#     labels = ['r_x', 'r_y', 'r_z']

#     fig = plt.figure(figsize=(10.5, 16))

#     ax1 = fig.add_subplot(211)
#     ax1.plot(offsets / 1000, rx, line_styles[0], color=colors[0], linewidth=1.2, label=labels[0])
#     ax1.plot(offsets / 1000, ry, line_styles[1], color=colors[1], linewidth=1.2, label=labels[1])
#     ax1.plot(offsets / 1000, rz, line_styles[2], color=colors[2], linewidth=1.2, label=labels[2])
#     ax1.set_xlabel('Offset [kHz]', fontsize=14, fontweight='bold')
#     ax1.set_ylabel('Rotation vector', fontsize=14, fontweight='bold')
#     ax1.legend(loc='center left', fontsize=14)

#     ax2 = fig.add_subplot(212, projection='3d')
#     for i in range(len(Lxx)):
#         ax2.quiver(0, 0, 0, Lxx[i], Lyy[i], Lzz[i], color=[0.12, 0.7, 0.6], linewidth=0.5)

#     ax2.set_xlabel('x', fontsize=14, fontweight='bold')
#     ax2.set_ylabel('y', fontsize=14, fontweight='bold')
#     ax2.set_zlabel('z', fontsize=14, fontweight='bold')
#     ax2.set_xlim([-1, 1])
#     ax2.set_ylim([-1, 1])
#     ax2.set_zlim([-1, 1])
#     ax2.grid(True, which='minor')

#     plt.savefig(f'figures/rotation_axis/axis_{name}.pdf')