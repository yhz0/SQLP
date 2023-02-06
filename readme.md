# SD Algorithm

A multi-epigraph version of two-stage regularized SD algorithm.

## Difference from the original SD algorithm

1. This version supports weighted scenario, and multiple weighted epigraph variables. Each epigraph variable maintains its own cut pool. This can be used in conjunction with variance reduction methods.
2. The incumbent cut is separately maintained. The incumbent cut is re-generated each time, and the old one is discarded immediately. In original SD, when the incumbent changes, the cut is still in the pool, but scaled down and then (probably) dropped at a later time when it is slack.
3. The master program sets LOWER BOUND for epigraph variables. In other words, the convex piecewise linear function approximation on the epigraph always includes an implied "horizontal" lower bound piece. When evaluating the cuts, all LOWER BOUNDS on epigraphs are also enforced. In original SD, this piece is not recorded.
