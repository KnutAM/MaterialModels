The crystal plasticity model is made up by several modules that can be combined to form the complete model. 
Please see ../ref.bib for which reference if you either directly use this model, or base your work upon it.

# Single crystal model
AceGen is used to generate the different variations of the single crystal models. In the folder ../AceGen all files to generate the variations of the local model are located. The settings are found in MainFile.nb and are described in the subsections below. By running this file, a fortran subroutine is generated. This can be put in the ../src folder in order to create a model for that particular type. Several pre-generated codes are already available, and these are described at the end of this document. 

## Model Name
This parameter only controls the output name of the .f90 file that should be copied to the ../src folder

## Crystal type
Controls the number and orientations of the slip systems. Currently, BCC, FCC, and BCC12 are available. The latter indicate the BCC crystal with only the 12 slip systems {110}<111>, while the regular BCC crystal contain all 48 slip systems. The FCC crystal only has 12 slip systems. This is defined in "CrystalPlasticity.nb"
## Elastic law
Either St Venant Isotropic law using the Green Lagrange strain on the intermediate configuration, or the same law but using a crystal elasticity tensor with cubic symmetry (3 parameters). This is defined in "ElasticLaw.nb".

## Isotropic hardening
Isotropic hardening law (can also include cross-hardening) for the slip systems. Several different hardening laws are available, and more can be added if desired. These are implemented in "IsotropicHardening.nb"

## Kinematic hardening
Kinematic hardening on each slip system. Two different numerical implementations of the Armstrong-Frederick type is included, using either a backward Euler time integration or an exact solution. They are defined in "KinematicHardening.nb"

## Overstress function
Again several different viscoplastic overstress functions are implemented. These are defined in "OverstressFunction.nb". 

# Homogenized response
This particular implementation homogenize multiple crystal orientations using the Taylor (or Voigt) assumption of a homogeneous strain field. When compiling the model, one can choose between different modules for orientations, which are located in the src folder. One example is the ../src/uniform_mod.f90 which contains 2, 4, 8, 16, 32, 64, 128 and 256 uniformly oriented crystals. For all cases, the number of orientations is set as a material parameter. To avoid excessive large binary files, it was avoided to put all orientations in one module, and different models should be compiled depending on use. 

# Compiling a model

## Abaqus
The file ../src/abaqus_compile_files.txt contain the compilation information used to compile for Abaqus. Each line describes one model to be compiled (please see the overall readme for this repository for how to compile for Abaqus). To add a model, the following must be specified: 
- *CompiledModelName*: The name of the resulting library
- *SingleCrystalModel*.f90: The file generated from AceGen for the single crystal model
- *crystal_orientation_mod*.f90: The module file for the particular desired orientations

Hence, replace these in the following line:
*CompiledModelName*,"../../umat_utils/SolveMatrixEquation.f90","../../umat_utils/smsutility.f90","../../umat_utils/Tensors_module.f90",*SingleCrystalModel*.f90,"umat_module.f90","umat.f90",*crystal_orientation_mod*.f90,"taylor_model.f90"

## Stand-alone shared/dynamic library (for e.g. matmodfit)
The file ../src/cmake_compile_models.txt contain the compilation information used to compile as a standalone .dll or .so for use by other softwares, such as matmodfit. Each line should be of the following format:
*CompiledModelName*:-D model=*SingleCrystalModel* -D orientation=*crystal_orientation_mod*
where the parts indicated by the *parameter* notation is described in the above subsection. 