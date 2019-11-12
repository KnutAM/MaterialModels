# MaterialModels
This repository contains user material models for Abaqus (UMAT) written by Knut Andreas Meyer. 
It with different model categories, located in the models folder. Each category contain (at least) three items:
- src: A folder containing the source code necessary to build the model
- doc: A folder containing documentation for the models
- ref.bib: Latex bibliography file with reference to the appropriate paper(s)

If you publish any results based on any of these models, please cite the relevant paper. 

Additionally, the models folder contain a umat_utils folder, with general purpose code used by several models. 

The scripts folder contain two scripts, one for building models for Abaqus (UMAT) and one for building a general shared library with the umat interface. These scripts take the following input:
Input 1: Which model category to build (see folder names in the models folder)
Input 2-: Which model number(s) to build within the given category (see cmake_compile_models.txt or abaqus_compile_files.txt) in the corresponding src folder

The build models for Abaqus and for the general case are put in folders compiled_abaqus and compiled respectively. 

## AceGen
Some of the included code is generated using AceGen, and this should be referenced according to the instructions on their homepage (http://symech.fgg.uni-lj.si/). See also Korelc J., (2002), Multi-language and Multi-environment Generation of Nonlinear Finite Element Codes,  Engineering with Computers, 2002, vol. 18, n. 4,str. 312-327
