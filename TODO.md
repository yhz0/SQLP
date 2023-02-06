1. (Performance issue) Improve data structure for dual vertices. 
Currently the dual vertices are stored with Set{Vector{Float64}}. This
compares two vectors exactly and creates unnecessary copies.
    - One way is to use norm and cosine similarity to define is_equal(), then
    hash them with the norm.
    - Alternative is to use local sensitive hashing.

2. Implement crash methods to initialize an appropriate starting solution.

3. (Performance Issue) Profiling needed. Probably the subproblem solving is slow because
the problem is passed to the solver again and again. Should reuse instance.

4. Need to implement stopping criteria.

5. Validation phase should use number of samples to target a specified level of accuracy
in the estimates.