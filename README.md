# MaxCutPlanarAndToroidal

Upper bounds and approximated solutions for GSet toroidal Max-Cut problems.

## Overview

This project provides algorithms and tools for computing upper bounds and finding approximated solutions to the Maximum Cut (Max-Cut) problem on graphs embedded in toroidal and planar topologies. It includes implementations for analyzing graph properties, computing theoretical upper bounds, generating solutions using topological methods, and verifying solution quality on GSet benchmark problems.

## Project Structure

### Main Scripts

- **`main_maxcut_bound.py`** - Calculates upper bounds for toroidal Max-Cut problems
  - Loads GSet problems from benchmark datasets
  - Computes topological properties (cycle structure, average cycle distance)
  - Generates and logs upper bound estimates with timing information

- **`main_maxcut_solve.py`** - Finds approximated solutions to toroidal Max-Cut problems
  - Implements topological solving strategies
  - Generates and outputs solution masks and groupings in hexadecimal format
  - Includes problem analysis and graph visualization capabilities
  - Saves results to the `sols/` directory
  - Logs detailed metrics including odd vertices and active short paths counts

- **`main_maxcut_verify.py`** - Verifies and validates Max-Cut solutions
  - Loads computed solutions from hex format
  - Calculates actual cut values for each solution
  - Detects and filters duplicate/complementary solutions
  - Reports solution statistics per GSet instance

### Core Libraries

- **`toroid_maxcut_lib.py`** - Main library for toroidal Max-Cut algorithms
  - Upper bound computation (`get_toroid_max_cut_bound`)
  - Topological solving methods (`get_toroid_max_cut_topological`)
  - Direct painting techniques

- **`toroid_faces_lib.py`** - Face-based algorithms for toroidal graphs
- **`toroid_cycles_lib.py`** - Cycle-based algorithms for toroidal graphs
- **`toroid_lp_lib.py`** - Liers-Pardella (LP) algorithm adapted to toroidal graphs (`get_toroid_max_cut_lp`)
  - Maximum-weight perfect matching on the expanded (K4) dual of the toroidal embedding gives an upper bound
  - If the matched edges are not a cut, it is repaired and improved by re-solving the matching with the edges of two non-contractible cycles fixed
  - The cut is proven optimal when it reaches the bound, otherwise it is a heuristic solution

- **`planar_maxcut_lib.py`** - Algorithms optimized for planar graphs
- **`graph_lib.py`** - General graph utilities and operations

### Utility Modules

- **`file_lib.py`** - File I/O operations for loading/saving GSet problems and solutions
- **`params_config.py`** - Configuration parameters for problem selection and algorithm settings
- **`params_const.py`** - Constant definitions
- **`params_manager.py`** - Parameter management interface

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
- **numpy** (2.3.4+) - Numerical computations
- **scipy** (1.16.2+) - Scientific computing utilities
- **networkx** (3.5+) - Graph algorithms and manipulation
- **matplotlib** (3.10.7+) - Visualization and plotting
- **jupyter** - Interactive notebooks (optional)

## Configuration

Edit `params_config.py` to configure the solver:

```python
# Toroidal graphs: [11, 12, 13, 32, 33, 34, 48, 49, 50, 57, 62, 65, 66, 67, 72, 77, 81]
gsetnum = [11]              # List of GSet problem numbers to solve (comment out to run all)
use_multi_proc = 1          # Enable multiprocessing (0 or 1)
ncycles = 18                # Number of cycles for algorithm
MinPathThreshhold = 5       # Threshold for path filtering (or 1E6 for full matrix)
alg_num = 3                 # Sum of the algorithms to run: 1 - Algorithm#1, 2 - Algorithm#2, 4 - LP (e.g. 5 - Algorithm#1 and LP, 7 - all)
lp_max_iters = 50           # LP: max number of improvement iterations
lp_max_stall = 10           # LP: stop after this number of iterations without improvement
```

### Supported Toroidal Instances

The following GSet instances are toroidal:
`[11, 12, 13, 32, 33, 34, 48, 49, 50, 57, 62, 65, 66, 67, 72, 77, 81]`

## Usage

### Computing Upper Bounds

To calculate an upper bound for GSet problems specified in `params_config.py`:

```bash
python main_maxcut_bound.py
```

**Output:**
- Console displays: `GSet #<N>, #Nodes=<V>, D=<distance>`
- `log.txt` records: `<timestamp> GSet=<N> MaxCutBound=<value> time=<duration>sec`

### Finding Approximated Solutions

To solve a Max-Cut problem and generate an approximated solution:

```bash
python main_maxcut_solve.py
```

**Process:**
1. Loads the GSet problems specified in `params_config.py`
2. Analyzes graph properties (vertices, edges, cycle structure, average distance)
3. Applies topological algorithms based on `alg_num` setting
4. Outputs solutions to the `sols/` directory in hexadecimal format
5. Logs detailed results to `log.txt`

**Console Output Example:**
```
GSet #11, #Nodes=800, D=20
max cut: 2891
masks array shape: (1, 800)
groups array shape: (1, 200)
GSetNum=11, #Nodes=800, #Edges=4694, D=20, #OddVertices=0, #ActShortPaths=1542, MC=2891, time=45 sec.
```

### Verifying Solutions

To verify the correctness of computed solutions:

```bash
python main_maxcut_verify.py
```

**Process:**
1. Reads hexadecimal solution files from `sols/` directory
2. Converts hex format back to binary groupings
3. Calculates actual cut values for each solution
4. Identifies and removes duplicate or complementary solutions
5. Reports statistics for each GSet instance

**Console Output Example:**
```
GSet #11 #Nodes=  800 D= 20 [ 1/ 3] [2891, 2856, 2745]
```

## Output Format

### Solution Files

Solutions are saved in the `sols/` directory:

- **`gset_<N>_maxcut_masks.txt`** - Binary masks in hexadecimal format representing the cut partition
- **`gset_<N>_maxcut_groups_hex.txt`** - Hexadecimal representation of vertex groupings for solutions

**Format Details:**
- Each solution line is represented as a hexadecimal string
- Each hex character encodes 4 bits (2 vertices in the grouping format)
- Allows compact storage and efficient I/O operations

### Logs

Results are appended to `log.txt`:

**Bound Computation Format:**
```
<timestamp> GSet=<N> MaxCutBound=<value> time=<duration>sec
```

**Solution Finding Format:**
```
<timestamp>: GSetNum=<N>, #Nodes=<V>, #Edges=<E>, D=<distance>, #OddVertices=<ov>, #ActShortPaths=<asp>, MC=<max_cut_value>, time=<duration> sec.
```

When LP is enabled (`alg_num` includes 4), `main_maxcut_bound.py` also computes the LP bound, which does not depend on `MinPathThreshhold`. It is saved to `bounds/gset<N>_max_cut_bound_lp.txt` and logged as:
```
<timestamp> GSet=<N> Alg=LP MaxCutBound=<value> time=<duration>sec
```

In `main_maxcut_solve.py`, when LP is enabled it is run after Algorithms #1/#2, timed separately, and logged as:
```
<timestamp>: GSetNum=<N>, #Nodes=<V>, #Edges=<E>, D=<distance>, Alg=LP, MCBound=<bound>, #Iters=<iterations>, MC=<max_cut_value>, time=<duration> sec.
```
Its solutions are saved to `sols/gset_<N>_maxcut_masks_lp.txt` and `sols/gset_<N>_maxcut_groups_hex_lp.txt`, and a `run times` line comparing the algorithms is printed for each problem.

## Algorithm Details

### Toroidal Topology

The implementation handles graphs embedded on a toroid, which enables:
- Wrap-around edges in both horizontal and vertical directions
- Face-based cycle detection for topological constraints
- Exploitation of toroidal structure for efficient partitioning

### Solving Strategies

1. **Topological Method** - Exploits the toroidal embedding structure to partition vertices based on cycles and faces
2. **Direct Painting** - Uses graph coloring heuristics and vertex assignment strategies optimized for toroidal graphs

### Key Metrics

- **D (Average Distance)** - Mean cycle distance computed from connected vertex patterns
- **#OddVertices** - Count of vertices with odd degree in the induced subgraph
- **#ActShortPaths** - Number of active/significant shortest paths used in solution computation
- **MC (Max Cut Value)** - The actual cut size achieved by the solution

### Upper Bound Computation

Computes theoretical upper bounds on the maximum cut size using:
- Topological properties of the toroidal embedding
- Graph characteristics (density, degree distribution)
- Cycle structure and path constraints

## Performance Considerations

- Graph analysis includes cycle size computation for characterizing topological properties
- Solutions are computed using heuristic methods optimized for large sparse graphs
- Multiprocessing support available via `use_multi_proc` setting for handling multiple instances
- Results are cached where applicable to improve runtime performance
- Hexadecimal encoding provides efficient memory usage for storing large binary solutions

## Supported GSet Instances

This project specifically targets toroidal GSet instances:
- **Group 1:** 11, 12, 13
- **Group 2:** 32, 33, 34
- **Group 3:** 48, 49, 50
- **Group 4:** 57, 62, 65, 66, 67, 72, 77, 81

## References

This project works with [GSet benchmark problems](http://web.stanford.edu/~yyye/yyye/Gset/), a standard dataset for evaluating Max-Cut algorithms.

## License

[Specify your license here, if applicable]

## Author

[Your name/organization]

## Contributing

Contributions are welcome! Please ensure that:
- Code follows the existing style and conventions
- Solutions are validated with the verification script
- Results are logged appropriately in `log.txt`
- Configuration changes are documented in `params_config.py`

## Support

For issues, questions, or contributions, please open an issue or contact the project maintainer.
