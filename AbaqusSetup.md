# Setting up Abaqus to work with user subroutines
There are various challenges associated with setting up Abaqus to work with user subroutines. To make matters worse, the procedures vary between Abaqus releases. Instructions in this document should therefore not be considered final. They have not been throuroghly verified. But if you find an error in the instructions, please add that as an issue and I will try to update the instructions accordingly. 

The setup will consist of two parts. The first is related to the compilation and linking on a system level and the second is related to settings in Abaqus' environment files (`abaqus_v6.env`) files

## System setup
### Windows

To be completed

### Linux
To be completed

## Abaqus environment files
To check that the correct options are set in Abaqus do the following:
- Create a folder with a file named `abaqus_v6.env` containing only `print compile_fortran`
- Run `abaqus -information version` (you can run other abaqus commands also, this just gives reasonably little output) from the terminal/command prompt in that folder

The first line being returned should look something like `['ifort', '/fpp', '/heap-arrays:1', '/Qmkl:sequential', '/include:%I']` on Windows and similar on Linux (but with slash `/` replaced by hyphen `-`). Note that in this list there will likely be many more options as well, that is ok (and probably beneficial for performance etc.)

If the output does not contain the item `'/Qmkl:sequential'`, add (or append to) a file `abaqus_v6.env` in your home directory (`%HOME%` on Windows and `~/` on Linux). Add the line `compile_fortran.append('/Qmkl:sequential')` (replace `/` by hyphen `-` on Linux)

If your ifort version > 16, you will also need to add `compile_fortran.append('/nostandard-realloc-lhs')` to the `abaqus_v6.env` in your home directory.

If the command returns the error `NameError: name 'compile_fortran' is not defined` then you should create the file `abaqus_v6.env` in your home directory (`%HOME%` on Windows and `~/` on Linux). This file should contain `compile_fortran = ['ifort', '/fpp', '/heap-arrays:1', '/Qmkl:sequential', '/include:%I']`. If the file already exists, just append this line. However, if this is the case, you will probably need more variables to be defined as well and I would guess something is wrong with your installation. Hints for possible variables to be defined are `link_sl`, `link_exe` but the contents of these lists are system dependent and it is therefore not, to my knowledge, possible to give more specific advice.

