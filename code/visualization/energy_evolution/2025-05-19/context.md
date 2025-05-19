These data are extracted from a 50 optimization pass on the `112544.off` file.
You can find out the code that I was using by going on the last version of this
day on git. To be more precise I added a priority queue for the first and
second pass, where I prioritize the edges that are on the boundary, and if both
are on or are not I prioritize the one with highest maximal adjacent energy.
