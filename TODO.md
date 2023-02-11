## Bug/Incomplete Code

## Enhancement

1. (Performance issue) Improve data structure for dual vertices. 
Currently the dual vertices are stored with Set{Vector{Float64}}. This
compares two vectors exactly and creates unnecessary copies.
    - One way is to use norm and cosine similarity to define is_equal(), then
    hash them with the norm.
    - Alternative is to use local sensitive hashing.

3. (Performance Issue) Profiling needed. Probably the subproblem solving is slow because
the problem is passed to the solver again and again. Should reuse instance.

4. Need to implement stopping criteria.

5. Validation phase should use number of samples to target a specified level of accuracy
in the estimates.

6. Implement override scenario weight in add_scenario! to support importance sampling. (easy to do)

7. Implement Harsha's [random cost coefficients](https://doi.org/10.1287/ijoc.2019.0929).
This may need data structure changes, and redesign tests.

8. Implement SMPS sampling methods (antithetic, stratified).

9. For importance sampling, need to override the epigraph's total_weight when adding new samples.