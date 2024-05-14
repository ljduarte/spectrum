# Generate t and F test spectrum
This script is made to demonstrate the use of the spectral representation  of the t and F tests in a $2^2$ factorial design, as discussed in the paper:  
*Determination of individual metabolic abundance changes owing to environmental  impacts: t and F-distribution factorial design spectral representations.* by Gustavo G. Marcheafave, Leonardo J. Duarte, Elis D. Pauli, Ieda S. Scarminio and
Roy E. Bruns. 

## How to use it:
Simply call the script by doing:
```
python script.py
```
The script will ask for the data and the factorial design matrices.

## Structure of the data matrices
The data utilized in the paper is available in this repository (see data_matrix.csv). It consists of a set of 20 NMR spectra ranging from 0.5 ppm to 9.0 ppm and are organized as follow:
The first colum of the csv contains the variables of the spectra. The data are then organized in a $13501 \times 20$, resulting in a csv file with 21 columns. Each column, beside the first one, contain one NMR spectra.

## Structure of the factorial matrices
The 2^2 factorial design matrix is given bellow:

| Factor 1  | Factor 2 | Interaction 12 | 
| --------- |:--------:|:--------------:|
|    -1     | -1       |              1 |
|    +1     | -1       |             -1 |
|    -1     | +1       |             -1 |
|    +1     | +1       |             +1 |

If two replicate of the same experiment is performed, teh matrix becomes:

| Factor 1  | Factor 2 | Interaction 12 | 
| --------- |:--------:|:--------------:|
|    -1     | -1       |              1 |
|    -1     | -1       |              1 |
|    +1     | -1       |             -1 |
|    +1     | -1       |             -1 |
|    -1     | +1       |             -1 |
|    -1     | +1       |             -1 |
|    +1     | +1       |             +1 |
|    +1     | +1       |             +1 |

## Output

## Reference
