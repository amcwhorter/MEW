# MEW: Marked Edge Walk

This repository contains the code to replicate the results from "The Marked Edge Walk: A Novel MCMC Algorithm for Sampling of Graph Partitions" by Atticus McWhorter and Daryl DeFord.

## Overview

MEW (Marked Edge Walk) is an MCMC algorithm designed for sampling graph partitions. This implementation provides tools to run the algorithm on dual graphs, with specific examples using New Hampshire and Texas data.

## Repository Structure

### Main Files

- **`batched_wrapper_NH.jl`** - Wrapper for running multiple chains on the New Hampshire dual graph
- **`batched_wrapper_TX.jl`** - Wrapper for running multiple chains on the Texas dual graph
- **`main_NH.jl`** - Executable functions for running MEW on New Hampshire data
- **`main_TX.jl`** - Executable functions for running MEW on Texas data
- **`beano2.3.jl`** - Helper functions used by all Julia files above
- **`push_reader2.jl`** - Code to read observables from partitions
- **`ks1.jl`** - Code to calculate KS statistics from observables
- **`make_all_figures.jl`** - Regenerates Figures 3-8 from the data in `data/`

### Configuration Files

- **`Project.toml`** - Julia project dependencies
- **`Manifest.toml`** - Exact versions of dependencies

### Data Folders

- **`NH/`** - Shape files and dual graph data for New Hampshire
- **`Texas/`** - Shape files and dual graph data for Texas
- **`data/`** - Compact processed data backing each figure (see `make_all_figures.jl`)
- **`figures/`** - Output of `make_all_figures.jl`

## Getting Started

### Prerequisites

- Julia (version specified in Project.toml)
- Required Julia packages (automatically installed via Project.toml)

### Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/amcwhorter/MEW.git
   cd MEW
   ```

2. Activate the Julia environment and install dependencies:
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```

## Usage

### Running Single Chains

To run MEW on New Hampshire data:
```julia
include("main_NH.jl")
```

To run MEW on Texas data:
```julia
include("main_TX.jl")
```

### Running Multiple Chains

To run multiple chains on New Hampshire:
```julia
include("batched_wrapper_NH.jl")
```

To run multiple chains on Texas:
```julia
include("batched_wrapper_TX.jl")
```

### Recreating Paper Figures

Figures 3-8 can be regenerated directly from the data in `data/`:
```bash
julia make_all_figures.jl
```
This writes each figure as a PNG to `figures/`.

For figures not covered by `make_all_figures.jl`, after running chains, observables can be calculated with:
```julia
include("push_reader2.jl")
```
and KS statistics can be calculated with: 
```julia
include("ks1.jl")
```
