# Sparse-Affine-Projection-Algorithm
This repository includes the codes for the following paper [SAPA:Sparse Affine Projection Algorithm in ADMM-LP decoding of LDPC codes](https://ieeexplore.ieee.org/document/9817674/citations#citations), presented in CWIT '22, 
and part of my master's thesis, [Approximate and Randomized ADMM-LP Decoding Using the Geometric Information of the Parity Polytope](https://tspace.library.utoronto.ca/bitstream/1807/125679/1/Asadzadeh_Amirreza_202211_MAS_thesis.pdf).

## Background
This project is focused on decoding low-density parity-check (LDPC) codes using an alternating direction method of multipliers (ADMM) framework
for solving the linear-programming (LP) decoding problem. This algorithm was initially introduced in the paper [Decomposition Methods for Large Scale LP Decoding](https://ieeexplore.ieee.org/abstract/document/6595057) and was significantly improved in the papers
[the ADMM Penalized Decoder for LDPC Codes](https://ieeexplore.ieee.org/abstract/document/7456284) and [Hardware Based Projection onto the Parity Polytope and Probability Simplex](https://ieeexplore.ieee.org/abstract/document/7421292).

ADMM-LP decoding algorithm iteratively applies message passing decoding on the Tanner graph of LDPC codes, while stroing the residual information in the Lagrange multipliers.
As part of this algorithm, in each iteration, there exists $M$ projections onto a geometric object known as parity polytope, for all $M$ check nodes of the graph.
Such projections can be implemented in practice via a water-filling process, which includes sorting and thresholding operations.
The average time complexity of projecting a vector of size $d$ onto its parity polytope of dimension $d$ is $O(d \log(d))$, due to the sorting operation.
It is known that this projection step is the complexity bottleneck of ADMM-LP decoding, and simplifying this step results in the complexity reduction of the overall decoding algorithm.

## Contribution
In this project, we approximate the parity polytope projection step in order to simplify it.
To approximate such projection in dimension $d$, instead of projecting a vector onto the convex hull of its closest d vertices on the polytope, we project it onto the affine hull of the closest $\chi$ vertices.
The hyperparameter $\chi$ can be tuned and set to a value significantly smaller than the dimension, $d$.
Our experiments validate that our sparse affine projection achieves a linear time-complexity, i.e. $O(d)$, while the overall decoding performance of the ADMM-LP algorithm will not be affected given enough number of iterations.



