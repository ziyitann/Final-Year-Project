# Final Year Project R Code

This repository contains the R code used to generate the simulated triangles, scenario results, and sensitivity check reported in the dissertation.

## Required R packages
- `ChainLadder`
- `DCL`

## Notes
- When running the DCL and BDCL sections, the `DCL` package may return warnings that a suitable delay function cannot be extracted. This occurs in some sparse simulated scenarios and reflects instability in the fitted delay structure rather than a coding error.
- When running the CLM section, the `ChainLadder` package may return a warning that the loglinear model for estimating `sigma_n` is not appropriate, and that Mack's estimation method is used instead. This is handled automatically by the package and does not stop the script from running.
