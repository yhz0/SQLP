# SD Algorithm

A multi-epigraph version of two-stage regularized SD algorithm.

## Difference from the original SD algorithm

1. This version supports weighted scenario, and multiple weighted epigraph variables. Each epigraph variable maintains its own cut pool. This can be used in conjunction with variance reduction methods.
2. The incumbent cut is separately maintained. The incumbent cut is re-generated each time, and the old one is discarded immediately. In original SD, when the incumbent changes, the cut is still in the pool, but scaled down and then (probably) dropped at a later time when it is slack.
3. The master program sets LOWER BOUND for epigraph variables. In other words, the convex piecewise linear function approximation on the epigraph always includes an implied "horizontal" lower bound piece. When evaluating the cuts, all LOWER BOUNDS on epigraphs are also enforced. In original SD, this piece is not recorded.

## Current Issues (TODO List)

1. (Performance issue) Improve data structure for dual vertices. 
Currently the dual vertices are stored with Set{Vector{Float64}}. This
compares two vectors exactly and creates unnecessary copies.
    - One way is to use norm and cosine similarity to define is_equal(), then
    hash them with the norm.
    - Alternative is to use local sensitive hashing.

2. (Performance Issue) Profiling needed. Probably the subproblem solving is slow because
the problem is passed to the solver again and again. Should reuse instance.

3. Need to implement stopping criteria.

4. Validation phase should use number of samples to target a specified level of accuracy
in the estimates.

5. Implement override scenario weight in add_scenario! to support importance sampling. (easy to do)

6. Implement Harsha's [random cost coefficients](https://doi.org/10.1287/ijoc.2019.0929).
This may need data structure changes, and redesign tests.

7. Implement SMPS sampling methods (antithetic, stratified).

8. For importance sampling, need to override the epigraph's total_weight when adding new samples.