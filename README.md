# Analysis of Broadband Universal Rotation Pulses (for NMR Spectroscopy)

This repo contains different functions for the analysis of broadband Universal Rotation pulses used in NMR spectroscopy. 
This package was designed for a compact analysis of UR-pulses. However, it can also be used for the analysis of pulses that belong to another pulse class (e. g. BEBOP, xyBEBOP).

You can create graphs showing the Transfer Function of the pulses, or the total rotation axis for different offsets. 
The script will analyze all pulses in the DIRECTORY specified.

A sample input file showing our pulse format is available. 
It is important not to store any other file types in the same directory!


main.py:
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

lib/basis.py:

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

lib/calculations.py:

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

  
### Example plots:

![Universal Rotation Pulse](https://github.com/StellaSlad/analysis_of_univeral_rotation_pulses/blob/main/images/UR.jpg)
This is a universal rotation pulse, because the rotation axes are almost identical for different offsets.

Comment: need to find the original figure for a not-UR pulse (so far, only have a small figute from a word file..)
(figure discription: This is not a universal rotation pulse, because the rotation axes are different for different offsets.)

### Planned improvements:

- add function that can read .bruker input files
- add a function in the beginning that automatically recognizes files with other endings and puts them into the subdirectory
"other_files" or creates output message "Detected a file format that cannot be processed, file name: ... . Please move this file to a different
directory".
