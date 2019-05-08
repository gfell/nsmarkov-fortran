#  Markov-Chain Approximations for Life-Cycle Models (Fortran version)

This repository provides **Fortran** subroutines to construct Markov chain approximations of (non-stationary) AR(1) processes as described in the paper "Markov-Chain Approximations for Life-Cycle Models"  by Giulio Fella, Giovanni Gallipoli and Jutong Pan, _Review of Economic Dynamics_ 34, 2019 ([https://doi.org/10.1016/j.red.2019.03.013](https://doi.org/10.1016/j.red.2019.03.013)). 

## Contains

- *mod_lcrouwenhorst.f90* subroutine to discretise a non-stationary AR(1) using our extension of Rouwenhorst [1995. "Asset pricing implications of equilibrium business cycle models," in `Frontiers of business cycle research', T. F. Cooley ed., Princeton University Press, Chapter 10.]
- *mod_lctauchen.f90* subroutine to discretise a non-stationary AR(1) using our extension of Tauchen [1986. "Finite State Markov-Chain Approximations to Univariate and
                  Vector Autoregressions," _Economics Letters_ 20].

