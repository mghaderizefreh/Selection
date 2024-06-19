# T5.1
This program produces two executables. The first is `analysis.exe` and the second is `selection.exe`.

## analysis.exe
`analysis.exe` is capable of estimating breeding values using either a random regeression model or a single trait.

### Random regression
The first type is a dataset of phenotyped individuals at different environments to estimate variance components and breeding values using an animal model. The model that is fitted is

	y{j,k} = (mu_int + mu_slo * X{k}) + (A_int{j} + A_slo{j} * X{k}) + (E_int{j,k} + E_slo{j,k} * X{k}).
	
In the above linear random-regression model, `y_{j,k}` is the phenotype of animal `j` and environmental level `k`. The challenge level for environment `k` is denoted by the covariate `X{k}`. The three terms on the right-hand side of the above equation correspond to the mean (`mu`), genetic deviation (`A`) and the residual (`E`). Each term is decomposed into a part for intercept (`int`) and slope (`slo`). Hence, the two fixed effects are the mean intercept and the mean slope. The random effects are genetic and environmental deviation bor both intercept and slope.

In the above model, it is assumed that the residual is heterogeneous. However, no correlation is assumed between the terms in the residual (`E_int{j,k}` and `E_slo_{j,k}`). The genetic terms, `A_int{j}` and `A_slo{j}`, may be correlated. The covariance structure for the genetic terms has to be specified in the inputs (See Usage). 

### Single trait
In the second type, the locations are ignored and the aim is to estimate breeding values for individual using a single trait model

  y{j} = mu + A_j + E_j

where `y_j` is the phenotype of animal `j`, `A_j` is its breeding value and `E_j` is its residual.

## selection.exe
This program simulates selection for a population of animals based on a selection index from a reaction norm model or a single trait analysis. In the case of reaction norm model the environmentla challenge levels may be known or unknown. The inpus for this program must be passed as a text file. See `test` folder for an example.

## version number
V2.5.3

## Copyright
* A copyright message and something about SMARTER project/ Roslin Institute/ University of Edinburgh/ us in general

# Installation
This software is written under FORTRAN90. A cmake project is accompanied that makes it possible to build and compile in any operating system with cmake, makefile and a Fortran compiler.
* Cmake project currently has an error. Use `manual.sh` with appropriate changes for compiler and flag for building.

## Dependencies
The main files are `analysis.f90` and `selection.f90` in `src/reml` folder. Other than the accompanied `f90` files with this release BLAS and LAPACK libraries are required for the compilation. If required, these libraries may be downloaded from [here](http://www.netlib.org/blas/blas.tgz) for BLAS and [here](http://www.netlib.org/lapack/lapack.tgz) for LAPACK and compiled based on the instructions from [here](https://gcc.gnu.org/wiki/GfortranBuild) (for Gfortran). For Intel Fortran these libraries are automatically installed.

## CMake
CMake version 2.8.5 and higher is required to build and compile the project. Navigate to the main folder and create an empty direcetory there, say `build`. Then navigate to `build` in terminal and enter

``` shell
cmake ..
make
```
CMake picks up a compiler and builds the executable. If you have multiple compilers, for example the GNU compiler and the Intel one, and wish to compile with Intel, then use
``` shell
FC=ifx cmake ..
```

The installation can be verified by running the command
``` shell
make test
```
inside the folder build. 

The executables `analysis.exe` and `selection.exe` will be placed in the `bin` folder.

Alternatively, run the executable `restart.sh`.

# Usage
Both program can be executed by simply running the name of the executables on the terminal/command line. The program `analysis.exe` is interactive and requires each input to be given separately. It is also possible to store the inputs in a separate file (e.g. `inputs.txt`) and use the input redirection operator `<`.

## analysis.exe
There are 6 different examples in test folder `input1` to `input6` in the `tests` folder. 

## selection.exe
Due to high number of inputs for program `selection.exe`, this program is not interactive and requires a ready input file to be passed at the beginning. A commented example for `selection.exe` is `input11` in the `tests` folder

## Outputs
`analysis.exe` provides estimated random and fix effects and variance components in the provided files.

`selection.exe` gives information on the genetic gain of breeding values. 


