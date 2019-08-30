# CUDA Computational Electrodynamics Sample & Gallery

## Non Technical Summary

This is a code sample provided to demonstrate the author's ability to provide programming that is both high-performance and academically rigorous enough to be used in a scientific context.    It also specifically demonstrates the author's ability to generate novel implementations synthesized from multiple scientific papers.

It also includes a gallery that demonstrates an eye for visuals, like this totally gratuitous and very well anti-aliased visualization of a tachyon' shockwave:

<img src="misc/tachyon.png" alt="tachyon" width="400"/>
 
As you can imagine, this sort of code is a complicated web that is generally difficult to produce without a litany of tiny bugs that result in wildly bad simulation results.  A gallery of these results is provided to demonstrate the author's ability to solve complex problems with creative methods:

<img src="blowups/blowup17.png" alt="blowup17" width="200"/>
<img src="blowups/blowup5.png" alt="blowup5" width="200"/>

## Technical Summary
Provided is an implementation of the [WCMOD7 algorithm](http://www.sciencedirect.com/science/article/pii/S0021999107002227) in CUDA 7 for use in a high-speed FDTD simulation on the `Yee Lattice`.  Specifically, this code is required to achieve O(n^2) accuracy because [proper sub-pixel materials anti-aliasing is non-trivial](https://dspace.mit.edu/handle/1721.1/49474) and results in anisotropic material boundaries.

Included is a graph demonstrating that the traditionally O(N\^2) complexity algorith results in O(N\^2) accuracy gains using this result against known solutions.  A couple of different kernels were tested, but each of them fell within expected parameters:

<img src="gallery/graphs.png" alt="gooderror" width="300"/>

Finally, here's a visual comparison with output from a similarly rigorous CPU based algorithm, [meep](https://meep.readthedocs.io/en/latest/):

<img src="gallery/comparison2.png" alt="comparison" width="200"/>

## Bonus Content:

#### Nanophotonics experiment design and simulation graphics:

<img src="misc/ez1.png" alt="vis1" width="200"/>

#### Basic electrical engineering knowledge:

<img src="misc/servo-mk-1.png" alt="servo1" width="200"/>
<img src="misc/photo.JPG" alt="servo2" width="200"/>

#### Assorted contextless nanophotonics graphs:

<img src="misc/pdexample.png" alt="graph1" width="200"/>
<img src="misc/fieldmagout.png" alt="graph2" width="200"/>
<img src="misc/spectraoffres.png" alt="vis1" width="200"/>

#### More Fun Blowups:

<img src="blowups/blowup16.png" alt="blowup16" width="200"/>
<img src="blowups/blowup10.png" alt="blowup10" width="200"/>
<img src="blowups/blowup13.png" alt="blowup13" width="200"/>