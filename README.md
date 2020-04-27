# RandomRegression
This software analyses a dataset of phenotyped individuals at different environments to estimate variance components and breeding values using an animal model. The model that is fitted is

	y{j,k} = (mu_int + mu_slo * X{k}) + (A_int{j} + A_slo{j} * X{k}) + (E_int{j,k} + E_slo{j,k} * X{k}).
	
In the above linear random-regression model, `y_{j,k}` is the phenotype of animal `j` and environmental level `k`. The challenge level for environment `k` is denoted by the covariate `X{k}`. The three terms on the right-hand side of the above equation correspond to the mean (`mu`), genetic deviation (`A`) and the residual (`E`). Each term is decomposed into a part for intercept (`int`) and slope (`slo`). Hence, the two fixed effects are the mean intercept and the mean slope. The random effects are genetic and environmental deviation bor both intercept and slope.

In the above model, it is assumed that the residual is heterogeneous. However, no correlation is assumed between the terms in the residual (`E_int{j,k}` and `E_slo_{j,k}`). The genetic terms, `A_int{j}` and `A_slo{j}`, may be correlated. The covariance structure for the genetic terms has to be specified in the inputs (See Usage). 

## version number
This is version 1.2 of the software released on 27 April 2020.

## Copyright
* A copyright message and something about SMARTER project/ Roslin Institute/ University of Edinburgh/ us in general

# Installation
This software is written under FORTRAN90. A cmake project is accompanied that makes it possible to build and compile in any operating system with cmake, makefile and a Fortran compiler.

## Dependencies
The main file is `RRReml.f90` in `src/reml` folder. Other than the accompanied `f90` files with this release BLAS and LAPACK libraries are required for the compilation. If required, these libraries may be downloaded from [here](http://www.netlib.org/blas/blas.tgz) for BLAS and [here](http://www.netlib.org/lapack/lapack.tgz) for LAPACK and compiled based on the instructions from [here](https://gcc.gnu.org/wiki/GfortranBuild) (for Gfortran).

## CMake
CMake version 2.8.5 and higher is required to build and compile the project. Navigate to the main folder and create an empty direcetory there, say `build`. Then navigate to `build` in terminal and enter

``` shell
cmake ..
make
```
CMake picks up a compiler and builds the executable. If you have multiple compilers, for example the GNU compiler and the Intel one, and wish to compile with Intel, then use
``` shell
FC=ifort cmake ..
```

The installation can be verified by running the command (CURRENTLY DOES NOT WORK)
``` shell
make test
```
A successful compilation should output two separate messages for correlated and uncorrelated data. Failing to output the two messages indicate the installation is unsuccessful.

The executable `randreg.exe` will be placed in the `bin` folder.

# Usage
The program can be executed by simply running the name of the program `randreg.exe` on the terminal/command line. The program is interactive and requires each input to be given separately. It is also possible to store the inputs in a separate file (e.g. `inputs.txt`) and use the input redirection operator `<`. Hence,

``` shell
./randreg.exe < inputs.txt
```
runs the program with the inputs from `inputs.txt` file. Each input in the `input.txt` must be on a separate line. The input file does not accept any comments.

The first input is the name of file that contains phenotypes. The second input is the name of the file that contains a relationship matrix. The program stops if any of these files are not in the correct format (see next section for a guide on format of the input file). The rest of the inputs are the initial guess for the estimation of variances. If it is assumed that there is no correlation between genetic part of slope and intercept, the value `0.0` must be given to the program. Otherwise, a non-zero value should be given as the correlation. The program stops if the initial guess for correlation is less than -1 or greater than +1.

## Required files for using and format
Two files are required for analysis: a phenotype file and a matrix file. All individuals that are phenotyped must have their relationship to other individuals explained in the matrix file. There is no minimum number of records per individual as long as the dataset is not too sparse. 

### 1. phenotype file
The phenotype contains information about the phenotype of animals. The file must contain 3 columns. The first column contains id of the animals and is integer. The second column is the level of challenge each animal is facing. If there are multiple records for one animal at different challenge levels, multiple lines should be used with each line for one challenge level. Finally, the third column contains the phenotypic values of the animal at the given challenge level. The second and third column are numeric values.

### 2. matrix file
The matrix file is the Numerator Relationship Matrix (NMR) or the Genomic Matrix Relationship (GMR). NMR may be used for analysing the data based on pedigree and GMR can be used when individuals are genotyped. These matrices should be calculated before running the program. The matrices need to be stored in the packed format because these matrices are symmetric. Hence, the matrix file contains three columns. The first and second column are the index of row and column in matrix, respectively and the value of the matrix at that row and column is given in the third column. Since, the matrices are symmetric, for a _n x n_ matrix, there are _n(n+1)/2_ lines in the matrix file.

## Outputs
The logl for each iteration is printed on the screen as the program runs. In addition, three output files are created. These files include:

### 1. `variances`
This file contains estimated variances and the correlation (if the software has been told that there exists a correlation between intercept and slope). This file has headers, therefore, it is self-explanatory.

### 2. `fixedEffects`
This file contains the estimated fixed effect values.

### 3. `randomEffects`
This file contains the estimated random effects for genetic part of slope and the genetic part of the intercept, in that order for all individuals that their relationship is explained in the 'matrix file'. Since we assume the residual has heterogeneous form, it also estimates the breeding values for environmental variance of slope for each phenotyped individual at each environmental level. The latter is appended to the bottom of the file `randomEffect`.

# Future Developments
The software at the moment is not capable of performing ssGblup. Hence, all the individuals need to be genotyped. This is the be developed in the future releases.

Some options, (e.g., maximum number of iterations, convergence criteria) are fixed and not allowed to be changed. This feature may be added in the next version.

