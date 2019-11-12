# Utility functions for umats
`smsutility.f90`          Needed for use with AceGen
`SolveMatrixEquation.f90` Solve nonlinear system of equations on the form R(x)=0 using Newtons method (Hence the matrix name)
`Tensors_module.f90`      Tensor utility functions

## Rules for how different utilities should be put in this folder
* Only one file should be required to include per utility
* If more files are needed these should be included explicitly
* These should be put in a folder with the same name as the main file.