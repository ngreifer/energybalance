
<!-- README.md is generated from README.Rmd. Please edit that file -->

# energybalance

`energybalance` implements energy balance, a method of estimating
weights for distribution-matching that can be used in a variety of
scenarios, including estimating treatment effects in observational
studies, assessing generalizability and transportability, and mediation
analysis. The method was originally described in Huling and Mak (2020).

The main function is `energybalance()`, which estimates a set of weights
that serve different purposes depending on the inputs given. The
wrappers `eb_ate()`, `eb_att()`, `eb_target()`, and `eb_mediation()`
make estimating weights more straightforward for the specific purposes
of estimating the ATE, the ATT, targeting a specific distribution, and
performing mediation analysis, respectively. The functions `edist()` and
`edist_ate()` compute the (weighted) energy distance between two groups.
