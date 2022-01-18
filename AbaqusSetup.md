# Setting up Abaqus to work with user subroutines
There are various challenges associated with setting up Abaqus to work with user subroutines. To make matters worse, the procedures vary between Abaqus releases. Instructions in this document should therefore not be considered final. They have not been throuroghly verified. But if you find an error in the instructions, please add that as an issue and I will try to update the instructions accordingly. 

The setup will consist of two parts. The first is related to the compilation and linking on a system level and the second is related to settings in Abaqus' environment files (`abaqus_v6.env`) files

## System setup
### Windows 10, Abaqus 2019
Instructions probably work for other versions as well

1)	Install [Visual Studio Community Edition 2022](https://visualstudio.microsoft.com/vs/community/)  
    Include “Desktop development with C++” under “Desktop&Mobile” under “Workloads”
    
2)	Install [Intel® oneAPI Base Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit.html)  
    Choose to connect with installed “visual studio 2022” and “Build Tools 2017” when asked

3)	Install [Intel® oneAPI HPC Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/hpc-toolkit.html)  
    Choose to connect with installed “visual studio 2022” and “Build Tools 2017” when asked
4)	Install Abaqus 2019 (including support for subroutines, i.e. CAA)
5)	Edit `C:\SIMULIA\Commands\abq2019.bat` (Change “2019” to your release. If you modify only abaqus.bat, then it might not work from inside CAE. At least for 2019, `abq2019.bat` is called when opening CAE, but not abaqus.bat) by adding the following on the line after `@echo off`  
    ```bat
    if not defined oneapi_setvars_has_been_called (
        @call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
        set oneapi_setvars_has_been_called=1
    )
    ```  
    The path to `oneAPI\setvars.bat` should be right, but be necessary to check this
6)	Open `cmd`, and run `abaqus verify –user_std` to check that everything is working.  
    If `abaqus` is not found, add `C:\SIMULIA\Commands` to your `path` system variable (“Start – Edit Environment Variables for your account”)


### Linux
To be completed, contributions are welcome!

## Abaqus environment files
To check that the correct options are set in Abaqus do the following:
- Create a folder with a file named `abaqus_v6.env` containing only `print compile_fortran`
- Run `abaqus -information version` (you can run other abaqus commands also, this just gives reasonably little output) from the terminal/command prompt in that folder

The first line being returned should look something like `['ifort', '/fpp', '/heap-arrays:1', '/Qmkl:sequential', '/include:%I']` on Windows and similar on Linux (but with slash `/` replaced by hyphen `-`). Note that in this list there will likely be many more options as well, that is ok (and probably beneficial for performance etc.)

If the output does not contain the item `'/Qmkl:sequential'`, add (or append to) a file `abaqus_v6.env` in your home directory (`%HOME%` on Windows and `~/` on Linux). Add the line `compile_fortran.append('/Qmkl:sequential')` (replace `/` by hyphen `-` on Linux)

If your ifort version > 16, you will also need to add `compile_fortran.append('/nostandard-realloc-lhs')` to the `abaqus_v6.env` in your home directory.

If the command returns the error `NameError: name 'compile_fortran' is not defined` then you should create the file `abaqus_v6.env` in your home directory (`%HOME%` on Windows and `~/` on Linux). This file should contain `compile_fortran = ['ifort', '/fpp', '/heap-arrays:1', '/Qmkl:sequential', '/include:%I']`. If the file already exists, just append this line. However, if this is the case, you will probably need more variables to be defined as well and I would guess something is wrong with your installation. Hints for possible variables to be defined are `link_sl`, `link_exe` but the contents of these lists are system dependent and it is therefore not, to my knowledge, possible to give more specific advice.

