!DEC$ FREEFORM
!Available models to be created by changing the basemod and ratemod variables in cmake
!Parameters ({*} are for rate dependent models
!Chaboche      ! props: E, v, Y0, Hiso, invYiso, {tstar, nexp,} Hk1, invYk1, Hk2, invYk2, ...                     (1-4 backstresses)
!Delobelle     ! props: E, v, Y0, Hiso, invYiso, {tstar, nexp,} delta, Hk1, invYk1, Hk2, invYk2, ...              (1-4 backstresses)
!OhnoWang      ! props: E, v, Y0, Hiso, invYiso, {tstar, nexp,} delta, Hk1, invYk1, mk1, Hk2, invYk2, mk2, ...    (1-4 backstresses)
!
!   Definitions of state variables nstatv = 2 + 6*nback {-1}
!   statev(1)   kappa           Isotropic hardening stress
!   statev(2)   mu              Plastic multiplier (only for rate independent models)
!   statev(3:)  beta            Back-stresses in groups of 6 components (statev(2:) for rate dependent models)
!      
!   Variables that always should be updated
!   stress      [1 x ntens]     stress tensor at beginning of increment, to be updated
!   statev      [1 x nstatv]    Array of solution dependent state variables, to be updated (and rotated using drot)
!   ddsdde      [ntens x ntens] The jacobian dsigma/epsilon. 
!   sse         [1]             Specific elastic strain energy
!   spd         [1]             Specific plastic dissipation
!   scd         [1]             Specific creep dissipation
!
!   Variables that should only be updated in fully coupled thermal-stress analysis
!   rpl         [1]             Volumetric heat generation per unit time at the end of increment (see guide if geostatic pores)
!   ddsddt      [1 x ntens]     Variation of the stress increments with respect to temperature
!   drplde      [1 x ntens]     Variation of rpl with respect to the strain increments
!   drpldt      [1]             Variation of rpl with respect to temperature
!     
!   Variables passed for information (DO NOT CHANGE)
!   stran       [1 x ntens]     Strains at the beginning of the increment. Strain components rotated to account for RBM before UMAT call
!   dstran      [1 x ntens]     Strain increments
!   time        [1 x 2]         Value of [step time, total time] at the beginning of the current increment
!   dtime       [1]             Time increment
!   temp        [1]             Temperature at the start of the increment
!   dtemp       [1]             Increment of temperature
!   predef      [1 x ?]         Array of interpolated values of predefined field variables at this point at the start of the increment
!   dpred       [1 x ?]         Array of increments of predefined field variables
!   cmname      [string]        User-defined material name, left justified. 
!   ndi         [int]           Number of direct stress components at this point
!   nshr        [int]           Number of engineering shear stress components at this point
!   ntens       [int]           Size of the stress or strain component array (ndi + nshr)
!   nstatv      [int]           Number of solution-dependent state variables that are associated with this material type
!   props       [1 x nprops]    User-specified array of material constants associated with this user material
!   nprops      [int]           User-defined number of material constants associated with this user material
!   coords      [1 x 3]         An array containing the coordinates of this point. 
!   drot        [3 x 3]         Rotation increment matrix
!   pnewdt      [1]             CAN BE CHANGED: Ratio of suggested new time increment to the time increment being used
!   celent      [1]             Characteristic element length
!   dfgrd0      [3,3]           Deformation gradient at the beginning of the increment
!   dfgrd1      [3,3]           Deformation gradient at the end of the increment
!   noel        [int]           Element number
!   npt         [int]           Integration point number
!   layer       [int]           Layer number (for composite shells and layered solids)
!   kspt        [int]           Section point number within the current layer
!   kstep       [int]           Step number
!   kinc        [int]           Increment number
! 
!
! Includes needed if compiled for abaqus 
! (Should always be commented away during commits)
! Note that order must be perserved for compilation to work
! For all models (ace_gen needs):
!include '../../umat_utils/smsutility.f90'
! 
! For specific models:
!include '<base_model>[_<rate_model>]_acegen_mod.f90'
!include '<base_model>[_rdep].f90'
! Where <base_model> is chaboche, delobelle or ohnowang
! And [] parts are only included for rate dependent models
! <rate_model> is norton, cowsym or delobelle
! 
! For all models:
!include '../../umat_utils/SolveMatrixEquation.f90'
!include 'gss_module.f90'

SUBROUTINE UMAT(stress,statev,ddsdde,sse,spd,scd, &
rpl,ddsddt,drplde,drpldt,&
stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,&
ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,&
celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
!DEC$ ATTRIBUTES DLLEXPORT :: UMAT
use SolveMatrixEquation
use gss_module
use model_module

implicit none
!
character*80 cmname
integer:: ndi,nshr,ntens,nstatv,nprops,noel,npt,&
layer, kspt, kstep, kinc
double precision:: dtime,temp,dtemp,pnewdt,celent,sse,spd,scd,rpl,drpldt
double precision:: stress(ntens),statev(nstatv),&
ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),&
stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),&
props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
!
!
double precision                :: stress_old(6)
double precision, allocatable   :: jac_inv(:,:), r0(:), mpar(:), x0(:), xold(:), param(:)
integer                         :: k1, info, max_iter
logical                         :: yielding, lconv
double precision                :: newton_tolerance
double precision                :: stressi(6), ddsddei(6,6) ! Internal stress and stiffness (always full tensor size)

if (ndi<3) then
    write(*,*) 'plane stress is not supported'
    stop
endif

stressi(1:ntens) = stress

newton_tolerance = 1.d-8*props(3)  !Scaling with initial yield limit because residual is a stress measure
max_iter = 20

! Check input
call checkinput(nprops, nstatv, props, mpar)
allocate(jac_inv(nvar, nvar))


stress_old = stressi
call elastic(mpar, statev, dstran, stress_old, stressi, ddsddei, yielding)

if (yielding) then
    !Solve for plasticity:
    
    !Set parameters for "Solve Matrix Equation"
    call sme_setparam(itertol=newton_tolerance, maxiter=max_iter)
    
    !Initiate starting guess and parameters
    call setup(x0, xold, param, stress_old, stressi, statev, dtime)

    !Solve local problem using newton raphson iterations and analytical jacobian
    call newton_raphson_ana(residual, anajac, x0, jac_inv, lconv, info, mpar, xold, param, printres=0)
    !call newton_raphson_num(residual, x0, jac_inv, lconv, info, mpar, statev, param, printres=1)
    
    !Calculate ddsdde, statev and stress from solution of local problem
    call plastic_output(ddsddei, statev, stressi, x0, jac_inv, mpar, stress_old, stran, dstran)
    
    !Check that solution converged
    call check_analysis(lconv, pnewdt)
    
endif

ddsdde = ddsddei(1:ntens, 1:ntens)
stress = stressi(1:ntens)

END SUBROUTINE
