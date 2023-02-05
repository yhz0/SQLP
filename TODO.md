1. (Performance issue) Improve data structure for dual vertices. New types should be defined with
overridden equal sign.
    - Currently the dual vertices are stored with Set{Vector{Float64}}. This
    compares two vectors exactly and creates unnecessary copies.
    - One way is to use norm and cosine similarity to define isequal, then
    hash them with the norm.
    - Alternative is to use local sensitive hashing.

2. Implement crash methods to initialize an appropriate starting solution.

3. (Performance Issue) Profiling needed. Probably the subproblem solving is slow because
the problem is passed to the solver again and again. Should reuse instance.