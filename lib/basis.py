import numpy as np

def create_single_spin_matrices():
    """
    Create matrices for a single spin 1/2.

    Returns:
    dict: Dictionary containing single spin matrices.
    """
    return {
        "mix": 0.5 * np.array([[0, 1], [1, 0]]),
        "miy": 0.5 * np.array([[0, -1j], [1j, 0]]),
        "miz": 0.5 * np.array([[1, 0], [0, -1]]),
        "mip": np.array([[0, 1], [0, 0]]),
        "mim": np.array([[0, 0], [1, 0]]),
        "mia": np.array([[1, 0], [0, 0]]),
        "mib": np.array([[0, 0], [0, 1]]),
        "ione": np.array([[1, 0], [0, 1]])
    }

def initialize_matrices(norder, nspins):
    """
    Initialize basis matrices.

    Parameters:
    norder (int): Order of the system.
    nspins (int): Number of spins.

    Returns:
    dict: Dictionary containing initialized matrices.
    """
    shape = (norder, norder, nspins)
    return {
        "iu": np.zeros(shape, dtype=complex),
        "ix": np.zeros(shape, dtype=complex),
        "iy": np.zeros(shape, dtype=complex),
        "iz": np.zeros(shape, dtype=complex),
        "ip": np.zeros(shape, dtype=complex),
        "im": np.zeros(shape, dtype=complex),
        "ia": np.zeros(shape, dtype=complex),
        "ib": np.zeros(shape, dtype=complex)
    }

def create_product_operator_matrices(nspins, matrices):
    """
    Create product operator matrices for each spin.

    Parameters:
    nspins (int): Number of spins.
    matrices (dict): Dictionary containing single spin matrices.

    Returns:
    dict: Dictionary containing product operator matrices for each spin.
    """
    ione = matrices["ione"]
    result = initialize_matrices(2**nspins, nspins)

    for ispins in range(nspins):
        dummy_matrices = {
            "dummy_u": ione, 
            "dummy_x": matrices["mix"], 
            "dummy_y": matrices["miy"], 
            "dummy_z": matrices["miz"],
            "dummy_p": matrices["mip"], 
            "dummy_m": matrices["mim"], 
            "dummy_a": matrices["mia"], 
            "dummy_b": matrices["mib"]
        }

        for j in range(1, nspins):
            for key, value in dummy_matrices.items():
                dummy_matrices[key] = np.kron(value, ione) if j >= ispins else np.kron(ione, value)

        for key, value in dummy_matrices.items():
            result[key.replace("dummy_", "")][:, :, ispins] = value

    return result

def basis(nspins):
    """
    Generates basis matrices for a given number of spins.

    Parameters:
    nspins (int): The number of spins.

    Returns:
    tuple: Basis matrices and the order of the system.

    Example usage: iu, ix, iy, iz, ip, im, ia, ib, norder = basis(2)
    """
    matrices = create_single_spin_matrices()
    product_matrices = create_product_operator_matrices(nspins, matrices)
    norder = 2 ** nspins

    return (product_matrices["iu"], product_matrices["ix"], product_matrices["iy"], 
            product_matrices["iz"], product_matrices["ip"], product_matrices["im"], 
            product_matrices["ia"], product_matrices["ib"], norder)