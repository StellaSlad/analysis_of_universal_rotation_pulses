import numpy as np

def create_single_spin_matrices() -> dict[str, np.ndarray]:
    """
    Create matrices for a single spin 1/2.

    Returns:
    dict: Dictionary containing single spin matrices.
    """
    return {
        "m_ix": 0.5 * np.array([[0, 1], [1, 0]]),
        "m_iy": 0.5 * np.array([[0, -1j], [1j, 0]]),
        "m_iz": 0.5 * np.array([[1, 0], [0, -1]]),
        "m_ip": np.array([[0, 1], [0, 0]]),
        "m_im": np.array([[0, 0], [1, 0]]),
        "m_ia": np.array([[1, 0], [0, 0]]),
        "m_ib": np.array([[0, 0], [0, 1]]),
        "i_one": np.array([[1, 0], [0, 1]])
    }

def initialize_matrices(norder: int, nspins: int) -> dict[str, np.ndarray]:
    """
    Initialize basis matrices.

    Parameters:
    norder (int): Order of the system.
    nspins (int): Number of spins.

    Returns:
    dict: Dictionary containing initialized matrices.
    """
    if norder <= 0 or nspins <= 0:
        raise ValueError("norder and nspins must be positive integers")

    shape = (norder, norder, nspins)
    return {
        "i_u": np.zeros(shape, dtype=complex),
        "i_x": np.zeros(shape, dtype=complex),
        "i_y": np.zeros(shape, dtype=complex),
        "i_z": np.zeros(shape, dtype=complex),
        "i_p": np.zeros(shape, dtype=complex),
        "i_m": np.zeros(shape, dtype=complex),
        "i_a": np.zeros(shape, dtype=complex),
        "i_b": np.zeros(shape, dtype=complex)
    }

def create_product_operator_matrices(nspins: int, matrices: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """
    Create product operator matrices for each spin.

    Parameters:
    nspins (int): Number of spins.
    matrices (dict): Dictionary containing single spin matrices.

    Returns:
    dict: Dictionary containing product operator matrices for each spin.
    """
    i_one = matrices["i_one"]
    result = initialize_matrices(2**nspins, nspins)

    for ispins in range(nspins):
        dummy_matrices = {
            "dummy_u": i_one, 
            "dummy_x": matrices["m_ix"], 
            "dummy_y": matrices["m_iy"], 
            "dummy_z": matrices["m_iz"],
            "dummy_p": matrices["m_ip"], 
            "dummy_m": matrices["m_im"], 
            "dummy_a": matrices["m_ia"], 
            "dummy_b": matrices["m_ib"]
        }

        for j in range(1, nspins):
            for key, value in dummy_matrices.items():
                dummy_matrices[key] = np.kron(value, i_one) if j >= ispins else np.kron(i_one, value)

        for key, value in dummy_matrices.items():
            result[key.replace("dummy_", "i_")][:, :, ispins] = value

    return result

def basis(nspins: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
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

    return (product_matrices["i_u"], product_matrices["i_x"], product_matrices["i_y"], 
            product_matrices["i_z"], product_matrices["i_p"], product_matrices["i_m"], 
            product_matrices["i_a"], product_matrices["i_b"], norder)