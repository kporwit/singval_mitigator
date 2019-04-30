Repository for singular value approach studies.

alpha_utils.py consists the original script from [1] converted to python. One must import this file and create_alpha function which has three parameters:
1. singular values list
2. eigenvalues list
3. verbosity

If verbosity is greater than 0 output from original script will be printed. 

create_alpha function returns lists of 1's if majorization is not satisfied.Otherwise matrix with prescribed eigenvalues and singular values is returned.
