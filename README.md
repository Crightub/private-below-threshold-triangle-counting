# Publishing Below-Threshold Triangle Counts under Local Weight Differential Privacy

This repository contains the reference implementation for the paper:

**Publishing Below-Threshold Triangle Counts under Local Weight Differential Privacy**

The code implements algorithms for privately releasing triangle count statistics under local weight differential privacy, including efficient sensitivity computation and optimized assignment routines. 
It also includes benchmarking code to reproduce the experimental results reported in the paper.

## Repository Structure
```
.
├── benchmark/ # Output for benchmark results
├── data/ # Input datasets 
├── include/ # Public C++ headers
├── src/ # Core implementation  
├── test/ # Unit and correctness tests
├── CMakeLists.txt
└── README.md
```
## Requirements

The code is written in C++ and requires the following dependencies:

- **CMake** (≥ 3.10)
- **C++20-compatible compiler** (GCC or Clang recommended)
- **Boost**: https://www.boost.org
- **nlohmann/json**: https://github.com/nlohmann/json
- **OpenMP**: https://www.openmp.org

On most systems, the dependencies can be installed via the system package manager.

## Build Instructions

We recommend an out-of-source build.
From the root directory:

```bash
mkdir build
cd build
cmake ..
make benchmark
```
## Executing Benchmarks

The benchmark executable supports two execution modes: a single-iteration benchmark and a multi-iteration benchmark.
The mode is determined by the number of command-line arguments.

### Single-Iteration Benchmark
Runs the benchmark once for a fixed set of parameters.

```bash
  ./benchmark <graph_name> <instances> <weight_eps> <count_eps> <lambda>
```

Arguments:
- graph_name: Name of the input graph (see supported graphs below). 
- instances: Number of independent benchmark instances to execute.
- weight_eps: Privacy parameter ε for releasing weight vectors.
- count_eps: Privacy parameter ε for releasing local counts.
- lambda: Threshold for triangles

Example:
```bash
  ./benchmark telecom-278 10 1.0 1.0 4
```

### Multi-Iteration Benchmark
Evaluates performance and accuracy over multiple parameter values.

```bash
./benchmark <graph_name> <instances_per_size> <iterations_count> <start_value> <step_size> <param> <weight_eps> <count_eps> <lambda>
```

Arguments:
- graph_name: Name of the input graph (see supported graphs below).
- instances_per_size: Number of benchmark instances executed per parameter setting.
- iterations_count: Number of iterations (parameter settings) to evaluate.
- start_value: Initial value of the iterated parameter.
- step_size: Increment applied to the parameter after each iteration.
- param: Name of the parameter to be varied.
- weight_eps: Privacy parameter ε for releasing weight vectors.
- count_eps: Privacy parameter ε for releasing local counts.
- lambda: Threshold for triangles

### Parameter Options
The following parameters are valid for the multi-iteration benchmark:
'weight_eps', 'count_eps', 'instance_size', 'lambda'

### Supported Graphs

Supported Graphs
The benchmark currently supports the following built-in graph datasets:
```aiignore
telecom-278
```
Telecommunications network graph with 278 nodes.
```aiignore
gmwcs   
```
Graph derived from the GMWCS benchmark dataset.
Graphs are loaded internally by name from the files stored in /data.

### AI Acknowledgement
GPT 5.2 was used while writing this code.