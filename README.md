# transcendentals
Implementations of transcendental functions (log, exp, sin, cos, etc)

Contained is a set of implementations of the transcendental math functions. The primary goal of this implementation is to produce deterministic results across multiple platforms. The approximations have been chosen (with very little rigor) to be a balance of speed vs accuracy suited for general game simulation code. In most cases, these should perform as fast or faster than "percise" CRT functions, but some functions (e.g. acos and asin) are less optimized right now. These implementations also don't handle every edge case required by the CRT standard.

Optimally, this gets expanded to include multiple implementations of each function at different orders of accuracy letting the user choose a speed tradeoff appropriate for their callsite.
