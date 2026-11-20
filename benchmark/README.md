# Benchmarks

The suite is defined in `benchmarks.jl`. Evaluation counts are read from
`params.json`. Regenerate that file after changing the benchmarks:

```
julia --project=benchmark benchmark/tune.jl
```

To compare a change against the previous state, save a baseline first:

```
julia --project=benchmark benchmark/save_baseline.jl
```

then make the change and run

```
julia --project=benchmark benchmark/compare.jl
```

The baseline goes to `_research/benchmark_baseline.json`, which is not tracked.

`sio_emission` downloads the ExoMol SiO EBJT line list on first use. Its
benchmarks are slow, so their evaluation and sample counts are fixed.

`forwarddiff_alloc` should stay within about 5× of `allocating`. Benchmarks of
the line shapes and synthesis methods are in SpectralLines.jl.
