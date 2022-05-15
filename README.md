# Repository for singular value approach studies.

## Dependencies
The main scripts needs these python packages installed:
1. [numpy](https://numpy.org/)

Please install it via pip3 before running the main scripts:
```
pip3 install numpy
```

## Main scripts
In this repository three main scripts exists:
1. s3mitigator.py
2. s23mitigator.py
3. s123mitigator.py

All of them accept three positional arguments:
1. name of the output file where the results will be saved,
2. number of matrix with experimental boundaries on T matrix (1 - m > EW ,2 - dm^2 > 100 eV^2, 3 - dm^2 ~ 0.1 - 1 eV^2) saved in the file experimental_boundaries.py,
3. number of random throws of the eigenvalues per one sifting procedure.

Scripts takes T matrix from the experimental_boundaries.py file and constructs the alpha matrices from them. After that for the given set of singular values (i.e. for the s3mitigator.py script which tries to minimalize the third singular value only):
1. `[1.0, 1.0, 0.9999]` for scenario 1,
2. `[1.0, 1.0 , 0.999]` for scenario 2,
3. `[1.0, 1.0, 0.999]` for scenario 3.
The singular value sets are defined as `s` lists in the scripts.


The script starts sifting procedure for given scenario and if the results for fixed singular value set is found then the script lowers the singular value by the `sjump` value defined in the scripts and the whole sifting procedure starts again. The minimalization of the singular value is bounded by the `decraesequantity` variable (i.e if it is set to two the sifting procedure will be run two times).

### Example run
To run the minimalization of third singular value for the mass scenario m > EW with 10^3 random throws of the eigenvalues and save the result to the s3.out run:
```
python3 s3mitigator.py s3.out 1 1000
```
Please remember that accepted arguments are **positional** which means that the correct order of the arguments is expected.


## alpha_utils
alpha_utils.py consists the original script from [http://dx.doi.org/10.1023/A:1021969818438](http://dx.doi.org/10.1023/A:1021969818438) converted to python. One must import this file and use `create_alpha` function which has three parameters:
1. singular values list
2. eigenvalues list
3. verbosity

If verbosity is greater than 0 output from original script will be printed. 

create_alpha function returns lists of 1's if majorization is not satisfied. Otherwise matrix with prescribed eigenvalues and singular values is returned.
