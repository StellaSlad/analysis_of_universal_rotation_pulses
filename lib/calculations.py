# Example usage:
# rot_axis('your_file_name', np.linspace(-132000 / 2, 132000 / 2, 551), 551)

import numpy as np
import scipy.linalg as scla
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def basis(nspins):
    # Matrices for single spin 1/2
    mix = 0.5 * np.array([[0, 1], [1, 0]])
    miy = 0.5 * np.array([[0, -1j], [1j, 0]])
    miz = 0.5 * np.array([[1, 0], [0, -1]])
    mip = np.array([[0, 1], [0, 0]])
    mim = np.array([[0, 0], [1, 0]])
    mia = np.array([[1, 0], [0, 0]])
    mib = np.array([[0, 0], [0, 1]])
    ione = np.array([[1, 0], [0, 1]])

    # Initializations
    norder = 2**nspins
    iu = np.zeros((norder, norder, nspins), dtype=complex)
    ix = np.zeros((norder, norder, nspins), dtype=complex)
    iy = np.zeros((norder, norder, nspins), dtype=complex)
    iz = np.zeros((norder, norder, nspins), dtype=complex)
    ip = np.zeros((norder, norder, nspins), dtype=complex)
    im = np.zeros((norder, norder, nspins), dtype=complex)
    ia = np.zeros((norder, norder, nspins), dtype=complex)
    ib = np.zeros((norder, norder, nspins), dtype=complex)

    # Creation of product operator matrices
    for ispins in range(1, nspins + 1):
        dummy_u = ione
        dummy_x = mix
        dummy_y = miy
        dummy_z = miz
        dummy_p = mip
        dummy_m = mim
        dummy_a = mia
        dummy_b = mib

        for j in range(2, nspins + 1):
            if j > ispins:
                dummy_u = np.kron(dummy_u, ione)
                dummy_x = np.kron(dummy_x, ione)
                dummy_y = np.kron(dummy_y, ione)
                dummy_z = np.kron(dummy_z, ione)
                dummy_p = np.kron(dummy_p, ione)
                dummy_m = np.kron(dummy_m, ione)
                dummy_a = np.kron(dummy_a, ione)
                dummy_b = np.kron(dummy_b, ione)
            else:
                dummy_u = np.kron(ione, dummy_u)
                dummy_x = np.kron(ione, dummy_x)
                dummy_y = np.kron(ione, dummy_y)
                dummy_z = np.kron(ione, dummy_z)
                dummy_p = np.kron(ione, dummy_p)
                dummy_m = np.kron(ione, dummy_m)
                dummy_a = np.kron(ione, dummy_a)
                dummy_b = np.kron(ione, dummy_b)

        iu[:, :, ispins - 1] = dummy_u
        ix[:, :, ispins - 1] = dummy_x
        iy[:, :, ispins - 1] = dummy_y
        iz[:, :, ispins - 1] = dummy_z
        ip[:, :, ispins - 1] = dummy_p
        im[:, :, ispins - 1] = dummy_m
        ia[:, :, ispins - 1] = dummy_a
        ib[:, :, ispins - 1] = dummy_b

    return iu, ix, iy, iz, ip, im, ia, ib, norder

def calculate_effective_propagator(name, offsrange, n_offsets):
    # spin system initialization
    nspins = 1
    iu, ix, iy, iz, ip, im, ia, ib, norder = basis(nspins)

    # read shaped pulse
    pulse = pd.read_table(name, skiprows=4, header=None, names=['Var1', 'Var2', 'Var3'], sep='\s+')
    x_amp = pulse['Var1'].astype(float).values
    y_amp = pulse['Var2'].astype(float).values
    timestep = float(pulse['Var3'][0])
    # alternativ: timestep = pulse['Var3'].astype(float).iloc[0]

    # offset grid
    offsets = np.linspace(-offsrange/2, offsrange/2, n_offsets)

    # initial state
    rhoinit = iz[:, :, 0]

    # number of time steps
    nsteps = len(x_amp)

    # preallocate space for results
    U_eff = np.eye(2)

    # loop over the offset grid
    for n in range(n_offsets):
        # initialize rho at t=0
        rho = rhoinit
        # chemical shift hamiltonian
        Hcs = 2 * np.pi * offsets[n] * iz[:, :, 0]

        # initial U_eff is the identity matrix
        U_eff_dummy = np.eye(2)

        # loop over time grid
        for k in range(nsteps):
            # grab the pulse amplitudes and set the control hamiltonian
            Hrf = 2 * np.pi * (x_amp[k] * ix[:, :, 0] + y_amp[k] * iy[:, :, 0])

            # set the total hamiltonian
            Hevo = Hcs + Hrf

            # calculate the time propagator
            uevo = scla.expm(-1j * timestep * Hevo)

            # use the time propagator to propagate the current state
            rho = np.dot(uevo, np.dot(rho, np.conj(uevo).T))

            # save the projections into the dummy variables
            U_eff_dummy = np.dot(U_eff_dummy, uevo)

        # grab the trajectory at the current offset point
        U_eff = np.concatenate((U_eff, U_eff_dummy), axis=1)

    np.savetxt(f'figures/offset/{name}_u.txt', U_eff, delimiter='\t', fmt='%.6f')

def final_states(name, n_offsets):
    nspins = 1
    iu, ix, iy, iz, ip, im, ia, ib, norder = basis(nspins)

    U = np.loadtxt(f'figures/offset/{name}_u.txt',dtype=np.complex128)
    
    Mxx, Mxy, Mxz = np.zeros(n_offsets), np.zeros(n_offsets), np.zeros(n_offsets)
    Myx, Myy, Myz = np.zeros(n_offsets), np.zeros(n_offsets), np.zeros(n_offsets)
    Mzx, Mzy, Mzz = np.zeros(n_offsets), np.zeros(n_offsets), np.zeros(n_offsets)

    for i in range(n_offsets):
        # Extracting the U_eff for the offset i
        U_eff = U[0:2, 2*i:2*i+2]

        # Propagation of the initial state ix
        rhox = np.dot(np.dot(U_eff, ix[:, :, 0]), np.conj(U_eff).T)
        rhoy = np.dot(np.dot(U_eff, iy[:, :, 0]), np.conj(U_eff).T)
        rhoz = np.dot(np.dot(U_eff, iz[:, :, 0]), np.conj(U_eff).T)

        # Normalized final magnetization components
        Mxx[i] = np.real(np.trace(np.dot(np.conj(rhox).T, ix[:, :, 0])) / np.trace(np.dot(ix[:, :, 0], ix[:, :, 0])))
        Mxy[i] = np.real(np.trace(np.dot(np.conj(rhox).T, iy[:, :, 0])) / np.trace(np.dot(iy[:, :, 0], iy[:, :, 0])))
        Mxz[i] = np.real(np.trace(np.dot(np.conj(rhox).T, iz[:, :, 0])) / np.trace(np.dot(iz[:, :, 0], iz[:, :, 0])))

        Myx[i] = np.real(np.trace(np.dot(np.conj(rhoy).T, ix[:, :, 0])) / np.trace(np.dot(ix[:, :, 0], ix[:, :, 0])))
        Myy[i] = np.real(np.trace(np.dot(np.conj(rhoy).T, iy[:, :, 0])) / np.trace(np.dot(iy[:, :, 0], iy[:, :, 0])))
        Myz[i] = np.real(np.trace(np.dot(np.conj(rhoy).T, iz[:, :, 0])) / np.trace(np.dot(iz[:, :, 0], iz[:, :, 0])))

        Mzx[i] = np.real(np.trace(np.dot(np.conj(rhoz).T, ix[:, :, 0])) / np.trace(np.dot(ix[:, :, 0], ix[:, :, 0])))
        Mzy[i] = np.real(np.trace(np.dot(np.conj(rhoz).T, iy[:, :, 0])) / np.trace(np.dot(iy[:, :, 0], iy[:, :, 0])))
        Mzz[i] = np.real(np.trace(np.dot(np.conj(rhoz).T, iz[:, :, 0])) / np.trace(np.dot(iz[:, :, 0], iz[:, :, 0])))

    M = np.column_stack((Mxx, Mxy, Mxz, Myx, Myy, Myz, Mzx, Mzy, Mzz))
    np.savetxt(f'figures/offset/{name}_results.txt', M, delimiter='\t', fmt='%.6f')

def calculate_quality_factor(name, n_offsets, desired_propagator):
    
    U = np.loadtxt(f'figures/offset/{name}_u.txt',dtype=np.complex128)

    transfer_efficiency = np.zeros(n_offsets + 1)
    for i in range(n_offsets):
        # Calculating the elements a and b
        effective_propagator = U[0:2, 2*i:2*i+2]
        # Extracting real and imaginary parts A, B, C, D
        transfer_efficiency[i] = np.trace(np.conj(desired_propagator).T @ effective_propagator)
    
    transfer_efficiency[n_offsets] = np.sum(transfer_efficiency) / n_offsets
    np.savetxt(f'figures/offset/{name}_te.txt', transfer_efficiency, delimiter='\t', fmt='%.6f')

def calculate_rotation_axis(name, offsets, n_offsets):
    # Read the file containing unitary rotation for different offsets
    U = np.loadtxt(f'figures/offset/{name}_u.txt',dtype=np.complex128)

    # Calculate elements a and b
    a = U[0, 0::2]
    b = U[0, 1::2]

    # Extracting real and imaginary parts A, B, C, D
    A = np.imag(b)
    B = np.real(b)
    C = np.imag(a)
    D = np.real(a)

    # Calculate half of the rotation angle omega (or omega/2)
    # Then, directional cosines are calculated - Lxx, Lyy, Lzz
    omega_half = np.arccos(D)
    S = np.sin(omega_half)

    # Directional cosines
    Lxx_new = A / S
    Lyy_new = B / S
    Lzz_new = C / S

    # Rotation vector
    rx_new = 2 * omega_half * Lxx_new
    ry_new = 2 * omega_half * Lyy_new
    rz_new = 2 * omega_half * Lzz_new

    # Save everything to lists
    Lxx, Lyy, Lzz, rx, ry, rz = [], [], [], [], [], []

    for i in range(n_offsets):
        Lxx.append(Lxx_new[i])
        Lyy.append(Lyy_new[i])
        Lzz.append(Lzz_new[i])
        rx.append(rx_new[i])
        ry.append(ry_new[i])
        rz.append(rz_new[i])

    axis = np.column_stack((Lxx, Lyy, Lzz))
    np.savetxt(f'figures/rotation_axis/{name}_axis.txt', axis, delimiter='\t', fmt='%.6f')

    # Plotting
    l = ['-', '--', '-.']  # line style
    li = [str(x) + 'k' for x in l]

    fig = plt.figure(figsize=(10.5, 16))

    ax1 = fig.add_subplot(211)
    ax1.plot(offsets / 1000, rx, li[0], color=[0.32, 0.9, 0.8], linewidth=1.2, label='r_x')
    ax1.plot(offsets / 1000, ry, li[1], color=[0.1, 0.6, 0.52], linewidth=1.2, label='r_y')
    ax1.plot(offsets / 1000, rz, li[2], color='k', linewidth=1.2, label='r_z')
    ax1.set_xlabel('Offset [kHz]', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Rotation vector', fontsize=14, fontweight='bold')
    ax1.legend(loc='center left', fontsize=14)

    ax2 = fig.add_subplot(212, projection='3d')
    for i in range(n_offsets):
        ax2.quiver(0, 0, 0, Lxx[i], Lyy[i], Lzz[i], color=[0.12, 0.7, 0.6], linewidth=0.5)

    ax2.set_xlabel('x', fontsize=14, fontweight='bold')
    ax2.set_ylabel('y', fontsize=14, fontweight='bold')
    ax2.set_zlabel('z', fontsize=14, fontweight='bold')
    ax2.set_xlim([-1, 1])
    ax2.set_ylim([-1, 1])
    ax2.set_zlim([-1, 1])
    ax2.grid(True, which='minor')
    
    plt.savefig(f'figures/rotation_axis/axis_{name}.pdf')