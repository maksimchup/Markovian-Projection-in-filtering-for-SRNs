# Markovian-Projection-in-filtering-for-SRNs

This repository contains an implementation of the Markov projection methods for the filtering problem for Stochastic Reaction Networks (SRNs). Two methods for dimensionality reduction are developed: standard Markovian Projection (MP) and Filtered Markov Projection (FMP).

A detailed description of the methods is provided in the paper:
> Hammouda, C. B., Chupin, M., Münker, S., & Tempone, R. (2025). Filtered Markovian Projection: Dimensionality Reduction in Filtering for Stochastic Reaction Networks. 
[arXiv preprint arXiv:2502.07918.]([http://www.example.com](https://arxiv.org/abs/2502.07918))

## Running simulations

Use the `main_bistable_network.m` or `main_linear_cascade.m` files to run the simulations.

## Code structure
The main functions and classes are listed below

- `@SRN` class stores the parameters of the stochastic process and provides functions to sample the trajectories using SSA and Tau-leap methods
- `@PF` class implements Particle Filter. Used at the first step of projection methods.
- `@FFSP` class implements the Filtered Finite State Projection method. Used to obtain the reference solution and at the second step of projection methods.

- `markovian_projection_extrap` function implements the standard MP and returns SRN of lower dimensionality
- `generate_observations` function generates input data for the filtering problem



## Potential improvements
- Add parallel implementation of the FFSP and PF (see `@FFSP/A.m` and `@PF/evolve.m`)
- Use more efficient extrapolation methods in `extrapolate_a_bar_`
