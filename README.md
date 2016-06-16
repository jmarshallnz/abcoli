## Multi-group epidemic model fit to Ecoli prevalence data via ABC-SMC

The `abc_multigroup.R` code is setup to do ABC fit of a multigroup SIR type
model to *E. coli* prevalence data.

The code isn't particularly well documented, but the bit down the bottom of
`abc_multigroup.R` is what does the running.

I seem to recall that `summary_stats.txt` is generated from an initial model
run without it being present. Am not sure what changes need to be made to
`multigroup()` function to allow for no summary statistics though.

The basic idea is that an initial ABC run is done to get an idea of areas
of the parameter space that are likely, and we then run a linear model
(see `find_summary` or similar) to get coefficients for summary statistics
which are basically powers of the difference in prevalence or some such.

My guess is that it needs to be done using the distance to the data, but
it's been quite a while since this was written.

Ideally the code would be packaged up so that the dll bit of it was just done
via Rcpp instead. I've included binaries of the .dll for running on windows and
the (64 bit?) .so for linux runs.

By the looks there is some nasty looking "parallel" code there where you can have
multiple runs of the code and each one waits for the others to complete before
moving on to the next particle refinement run.

Suggest starting by seeing if `multigroup` can be run by providing some parameter
values and getting a simulation out of it. The `summary_stats` bit of that will need removing
so I suggest breaking up `multigroup` into two functions: one that just runs the model
and the other that looks to compute the distance from the simulation to the data.

(I may well be getting some or all of the above wrong - it's been a long time!)

