# MaxCutPlanarAndToroidal

Upper bounds and approximated solutions for GSet toroidal Max-Cut problems.

## Overview

This project provides algorithms and tools for computing upper bounds and finding approximated solutions to the Maximum Cut (Max-Cut) problem on graphs embedded in toroidal and planar topologies. The implementation focuses on GSet benchmark toroidal problems: ##11,12,13,32,33,34,48,49,50,57,62,65,66,67,72,77,81, offering both theoretical bounds and practical heuristic solutions.

## Project Structure

### Main Scripts

- **`main_maxcut_bound.py`** - Calculates upper bounds for toroidal Max-Cut problems
  - Loads GSet problems from benchmark datasets
  - Computes topological properties (cycle structure, average distance)
  - Generates upper bound estimates

- **`main_maxcut_solve.py`** - Finds approximated solutions to toroidal Max-Cut problems
  - Implements topological and direct painting solving strategies
  - Generates and outputs solution masks and groupings in hex format
  - Includes problem analysis and visualization capabilities
  - Saves results to the `sols/` directory

- **`main_maxcut_verify.py`** - Verifies and validates Max-Cut solutions

### Core Libraries

- **`toroid_maxcut_lib.py`** - Main library for toroidal Max-Cut algorithms
  - Upper bound computation
  - Topological solving methods
  - Direct painting techniques

- **`toroid_faces_lib.py`** - Face-based algorithms for toroidal graphs
- **`toroid_cycles_lib.py`** - Cycle-based algorithms for toroidal graphs

- **`planar_maxcut_lib.py`** - Algorithms optimized for planar graphs
- **`graph_lib.py`** - General graph utilities and operations

### Utility Modules

- **`file_lib.py`** - File I/O operations for loading/saving GSet problems and solutions
- **`params_config.py`** - Configuration parameters for problem selection
- **`params_const.py`** - Constant definitions

### Directories

- **`data/`** - Contains GSet problem instances and benchmark datasets
- **`sols/`** - Output directory for computed solutions

## Installation

### Requirements

- Python 3.7+

### Dependencies

Install required packages using pip:

```bash
pip install -r requirements.txt
```

Key dependencies include:
- **numpy** - Numerical computations
- **scipy** - Scientific computing utilities
- **networkx** - Graph algorithms and manipulation
- **matplotlib** - Visualization and plotting
- **jupyter** - Interactive notebooks (optional)

## Usage

### Computing Upper Bounds

To calculate an upper bound for a specific GSet problem:

```bash
python main_maxcut_bound.py
```

Configure the problem number in `params_config.py` by setting `gsetnum` to the desired GSet instance.

### Finding Approximated Solutions

To solve a Max-Cut problem and generate an approximated solution:

```bash
python main_maxcut_solve.py
```

The script will:
1. Load the GSet problem specified in `params_config.py`
2. Analyze graph properties (number of nodes, edges, cycle structure)
3. Apply topological and painting-based algorithms
4. Output solutions to the `sols/` directory in hexadecimal format
5. Log results to `log.txt`

### Verifying Solutions

To verify the correctness of computed solutions:

```bash
python main_maxcut_verify.py
```

## Configuration

Edit `params_config.py` to select which GSet problem instance to work with:

```python
gsetnum = 12  # Change to desired GSet toroidal instance (e.g.: 11,12,13,32,33,34,48,49,50,57,62,65,66,67,72,77,81)
```

## Output Format

### Solution Files

Solutions are saved in the `sols/` directory:

- **`gset_<N>_maxcut_masks.txt`** - Binary masks representing the cut partition
- **`gset_<N>_maxcut_groups_hex.txt`** - Hexadecimal representation of vertex groupings

### Logs

Results are appended to `log.txt` with the following format:
```
GSet #<N>, #Nodes=<V>, #Edges=<E>, D=<avg_distance>, MC=<max_cut_value>, time: <duration> sec.
```

## Algorithm Details

### Toroidal Topology

The implementation handles graphs embedded on a toroid, which allows for:
- Wrap-around edges in both horizontal and vertical directions
- Face-based cycle detection
- Topological constraint satisfaction

### Solving Strategies

1. **Topological Method** - Exploits the toroidal structure to partition vertices
2. **Direct Painting** - Uses graph coloring and heuristic vertex assignment

### Upper Bound Computation

Computes theoretical upper bounds on the maximum cut size using topological properties and graph characteristics.

## Performance Considerations

- Graph analysis includes cycle size computation for topological properties
- Solutions are computed using heuristic methods optimized for large sparse graphs
- Results are cached where applicable to improve runtime performance

## References

This project works with GSet benchmark problems, a standard dataset for evaluating Max-Cut algorithms.

## License

[Specify your license here, if applicable]

## Author

[Your name/organization]

## Contributing

Contributions are welcome! Please ensure that:
- Code follows the existing style
- Solutions are validated with the verification script
- Results are logged appropriately

## Support

For issues, questions, or contributions, please open an issue or contact the project maintainer.
