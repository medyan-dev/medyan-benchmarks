In the previous versions of MEDYAN, the curvature computation in Mem3DG mode may lead to minimization problems, because the edge length may shrink to almost zero.

An attempt to fix it is to add a minimum length to edges (using something like a one-sided FENE potential).

The benchmark tests whether the new implementation is faster than the old one. To run the benchmarks, first go to this directory and prepare the environment using:

```shell
julia --project prepare.jl <old/new>
```

Then run
```shell
julia --project filopodia.jl
```
And check for printed messages.
