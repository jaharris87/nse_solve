NSE composition solver
============

Fast, thread-safe solver for composition in nuclear statistical equilibrium (NSE).
Solution is obtained by solving the Saha equation using globally convergent Newton method with backtracking line search.
Corrections from Coulomb interactions (screening) and nuclear partition functions are included.
For theoretical review, see Seitenzahl et al., ADNT 95 (2009) 96.

## Getting started

Configure options in `control` file, then run.


## Example(s)

`Data_SN160` is included as an example directory containing the necessary nuclear data for a 160 species reaction network suited for core-collapse supernovae.
Other data directories can be created using the build_net tool.
