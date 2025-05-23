# CrystalStructurePrediction

[![Build Status](https://github.com/Br0kenSmi1e/CrystalStructurePrediction.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Br0kenSmi1e/CrystalStructurePrediction.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Br0kenSmi1e/CrystalStructurePrediction.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Br0kenSmi1e/CrystalStructurePrediction.jl)

CrystalStructurePrediction is a Julia package for predicting crystal structures from a given set of ions and their charges.
The method is based on the integer programming formulation that detailed in [^Gusev2023].

## Usage
The following commands will clone the repository, initialize the project, and run the SrTiO3 example.
Just type the following commands in your terminal:
```bash
git clone https://github.com/Br0kenSmi1e/CrystalStructurePrediction.jl.git
cd CrystalStructurePrediction.jl
make init  # initialize the project
make run-example  # run the SrTiO3 example
```

If everything goes well, you should see the following crystal structure visualization:

![SrTiO3 crystal structure](examples/SrTiO3-structure.png)

The center is Ti, surrounded by Sr (corner) and O (face center).

## References

[^Gusev2023]: Gusev, V.V., Adamson, D., Deligkas, A., Antypov, D., Collins, C.M., Krysta, P., Potapov, I., Darling, G.R., Dyer, M.S., Spirakis, P., Rosseinsky, M.J., 2023. Optimality guarantees for crystal structure prediction. Nature 619, 68–72. https://doi.org/10.1038/s41586-023-06071-y
