# Analysis of Broadband Universal Rotation Pulses (for NMR Spectroscopy)

This repository contains different functions for the analysis of broadband Universal Rotation (UR) pulses used in NMR spectroscopy. 
This package was designed for a compact analysis of BURBOP pulses. However, it can also be used for the analysis of pulses that belong to another pulse class (e. g. BEBOP, xyBEBOP). You can create graphs showing the Magnetisation Transfer Function of the pulses and the total rotation axis for different offsets. 
The script will analyze all pulses in the DIRECTORY specified. A sample input file showing our pulse format is available.

Plotting the rotation vectors for spin with different offsets is an essential step for determining, whether a pulse performs a Univeral Rotation.
A UR pulse will rotate all the spins about the same axis. Optimised UR pulses show very small deviations from this behaviour.
Here, you can download BURBOP pulses optimised for different use cases: https://www.ioc.kit.edu/luy/111.php#burbop19f

A demo notebook is available [here](https://github.com/StellaSlad/analysis_of_universal_rotation_pulses/blob/main/demo.ipynb).

## Literature on fundamental theoretical concepts, BURBOP pulses and their applications:
1. Optimal Control of Coupled Spin Dynamics: Design of NMR Pulse Sequences by Gradient Ascent Algorithms:
https://www.ch.nat.tum.de/fileadmin/w00bzu/ocnmr/pdf/94_GRAPE_JMR_05_.pdf
2. Broadband 180◦ universal rotation pulses for NMR spectroscopy designed by optimal control:
<br /> https://arxiv.org/pdf/1111.6647
4. Exploring the limits of broadband 90° and 180° universal rotation pulses:
https://www.sciencedirect.com/science/article/abs/pii/S1090780712003126?via%3Dihub
5. Comprehensive and High-Throughput Exploration of Chemical Space Using Broadband 19F NMR-Based Screening:
https://onlinelibrary.wiley.com/doi/10.1002/ange.202002463

# Modules

## main.py

The `main.py` module serves as the entry point for running the quantum simulation and analysis tasks. It sets up the necessary parameters, ensures required directories exist, retrieves pulse filenames, and plots the results of the analysis.

### Key Components

- **Parameters**: Defines the key parameters for the simulation, such as the directory for pulse files, timestep, ideal axis, angle, offset range, number of offsets, desired propagator, and offsets array.

- **Functions**:
  - `ensure_directories()`: Ensures that the necessary directories for storing figures and results exist.
  - `get_pulse_filenames() -> list[str]`: Retrieves all filenames in the current directory that start with 'pulse'.
  - `plot_results(name: str, offsets: np.ndarray)`: Plots the results of the analysis and saves the plots to files.

- **Main Execution**:
  - The `main()` function orchestrates the execution of the simulation and analysis tasks. It ensures directories exist, retrieves pulse filenames, and iterates over each pulse file to perform calculations and plot results.

## basis.py:

The `basis.py` module is responsible for generating the basis matrices for a given number of spins in a quantum system. It provides functions to create single spin matrices and product operator matrices for multiple spins. These matrices are essential for simulating and analyzing quantum systems, particularly in the context of NMR (Nuclear Magnetic Resonance) and quantum computing.

### Functions

- `create_single_spin_matrices() -> dict[str, np.ndarray]`
  - Creates and returns a dictionary containing the matrices for a single spin 1/2 system. These matrices include the Pauli matrices and identity matrix.

- `initialize_matrices(norder: int, nspins: int) -> dict[str, np.ndarray]`
  - Initializes and returns a dictionary containing zero matrices for the specified order and number of spins.

- `create_product_operator_matrices(nspins: int, matrices: dict[str, np.ndarray]) -> dict[str, np.ndarray]`
  - Creates and returns product operator matrices for each spin in a multi-spin system. This function uses the single spin matrices to generate the product operators.

- `basis(nspins: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]`
  - Generates and returns the basis matrices for a given number of spins, along with the order of the system. The basis matrices include the identity matrix and the Pauli matrices for each spin.

### Example Usage

```python
from lib.basis import basis

# Generate basis matrices for a 2-spin system
iu, ix, iy, iz, ip, im, ia, ib, norder = basis(2)
```

## calculations.py

The `calculations.py` module is responsible for performing various calculations related to quantum systems, particularly in the context of NMR (Nuclear Magnetic Resonance) and quantum computing. It provides functions to read pulse files, calculate time propagators, effective propagators, final states, quality factors, and rotation axes.

### Functions

- `read_pulse_file(name: str) -> tuple[np.ndarray, np.ndarray, float]`
  - Reads a shaped pulse file and returns the pulse amplitudes (x and y) and the timestep.

- `calculate_time_propagator(timestep: float, Hevo: np.ndarray) -> np.ndarray`
  - Calculates and returns the time propagator for a given timestep and Hamiltonian.

- `calculate_effective_propagator(name: str, offsrange: float, n_offsets: int)`
  - Calculates the effective propagator for a given shaped pulse file over a range of offsets and saves the results to a file.

- `final_states(name: str, n_offsets: int)`
  - Calculates the final states for a given shaped pulse file over a range of offsets and saves the results to a file.

- `calculate_quality_factor(name: str, n_offsets: int, desired_propagator: np.ndarray)`
  - Calculates the quality factor for a given shaped pulse file over a range of offsets and saves the results to a file.

- `calculate_rotation_axis(name: str, offsets: np.ndarray, n_offsets: int)`
  - Calculates the rotation axis for a given shaped pulse file over a range of offsets and saves the results to a file.

- `plot_rotation_axis(name: str, offsets: np.ndarray, rx: list, ry: list, rz: list, Lxx: list, Lyy: list, Lzz: list)`
  - Plots the rotation axis and saves the plot to a file.

### Example Usage

```python
from lib.calculations import calculate_effective_propagator, final_states, calculate_quality_factor, calculate_rotation_axis

# Calculate the effective propagator
calculate_effective_propagator('pulse_file.txt', 10000, 101)

# Calculate the final states
final_states('pulse_file.txt', 101)

# Calculate the quality factor
desired_propagator = np.array([[np.sqrt(2)/2, np.sqrt(2)*1j/2], [-np.sqrt(2)/2, np.sqrt(2)*1j/2]])
calculate_quality_factor('pulse_file.txt', 101, desired_propagator)

# Calculate the rotation axis
offsets = np.linspace(-5000, 5000, 101)
calculate_rotation_axis('pulse_file.txt', offsets, 101)
```

### Planned improvements:

- add function that can read .bruker input files
- add a function in the beginning that automatically recognizes files with other endings and puts them into the subdirectory
"other_files" or creates output message "Detected a file format that cannot be processed, file name: ... . Please move this file to a different
directory".
