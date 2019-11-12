# MaterialModels
This repository contains user material models for Abaqus (UMAT) written by Knut Andreas Meyer. 
It with different model categories, located in the models folder. Each category contain (at least) three items:
- src: A folder containing the source code necessary to build the model
- doc: A folder containing documentation for the models
- ref.bib: Latex bibliography file with reference to the appropriate paper(s)

If you publish any results based on any of these models, please cite the relevant paper. 

Additionally, the models folder contain a umat_utils folder, with general purpose code used by several models. 

## Compiling the models
The scripts folder contain two scripts, one for building models for Abaqus (UMAT) and one for building a general shared library with the umat interface. These scripts take the following input:
Input 1: Which model category to build (see folder names in the models folder)
Input 2-: Which model number(s) to build within the given category (see cmake_compile_models.txt or abaqus_compile_files.txt) in the corresponding src folder

The build models for Abaqus and for the general case are put in folders compiled_abaqus and compiled respectively. 

## Overview of models
### GenSmallStrain
This category includes the "standard" plasticity models for small strain with non-linear isotropic and kinematic hardening (i.e. the Chaboche model). 
Different multiaxial extensions of the Armstrong-Frederick kinematic hardening rule is included in the form of the Ohno-Wang model, as well as the multiaxial modification based on Burlet and Cailletaud (1986)
The models are available in rate-independent as well as rate-dependent forms, the latter with three different overstress functions. 

### GenFiniteStrain
This category includes extensions of the rate-independent small strain models to finite strains, using a hyperelasto-plastic framework. The implementation is from Meyer et. al (2018)

## AceGen
Some of the included code is generated using AceGen, and this should be referenced according to the instructions on their homepage (http://symech.fgg.uni-lj.si/). See also Korelc (2002)

## References
* H. Burlet and G. Cailletaud, “Numerical techniques for cyclic plasticity at variable temperature,” Eng. Comput., vol. 3, no. 2, pp. 143–153, Feb. 1986.
* K. A. Meyer, M. Ekh, and J. Ahlström, “Modeling of kinematic hardening at large biaxial deformations in pearlitic rail steel,” Int. J. Solids Struct., vol. 130–131, pp. 122–132, 2018.
* Korelc J., (2002), Multi-language and Multi-environment Generation of Nonlinear Finite Element Codes,  Engineering with Computers, 2002, vol. 18, n. 4,str. 312-327