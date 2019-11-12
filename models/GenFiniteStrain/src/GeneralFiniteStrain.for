!DEC$ FREEFORM
!   Definitions of material parameters nprops = 6 + 3*nback
!   props(1)    G               Shear modulus
!   props(2)    K               Bulk modulus
!   props(3)    Y0              Initial yield limit
!   props(4)    Hiso            Isotropic hardening modulus
!   props(5)    invYiso         (1/Yiso), Yiso is isotropic saturation limit
!   props(6)    delta           Cailletaud fraction
!   props(7)    Hk,1            Kinematic hardening modulii nr. 1
!   props(8)    invYk,1         (1/Yk,2), Yk,1 is kinematic saturation limit 1
!   props(9)    m,1             Exponent in OhnoWang 
!   props(10)   Hk,2            Kinematic hardening modulii nr. 2
!   props(11)   invYk,2         (1/Yk,2), Yk,2 is kinematic saturation limit 2
!   props(12)   m,2             Exponent in OhnoWang 
!
!   Definitions of state variables nstatv = 10 + 9*nback
!   statev(1)   lambda          Total plastic multiplier
!   statev(2:10)Fp-I2           Plastic deformation gradient
!   statev(11:) Fk-I2           Kinematic ("elastic") deformation gradient
!      
!   ======  Start of constant umat header (Knut Andreas)  ========
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

SUBROUTINE UMAT(stress,statev,ddsdde,sse,spd,scd, &
rpl,ddsddt,drplde,drpldt,&
stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,&
ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,&
celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
!DEC$ ATTRIBUTES DLLEXPORT :: UMAT
use tensors_module
use GFS_module
use SolveMatrixEquation

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
!     ======  End of constant umat header excl. includes (Knut Andreas)  ========
double precision                :: Fnew(9)
double precision                :: param(9), x0(nstatv), ddsdde_elastic(6,6)
double precision, allocatable   :: jac_inv(:,:), r0(:)
integer                         :: k1, info, max_iter
logical                         :: yielding, lconv, use_el_stiff
double precision                :: newton_tolerance
double precision                :: stressi(6), ddsddei(6,6) ! Internal stress and stiffness (always full tensor size)

if (ndi<3) then
    write(*,*) 'plane stress is not supported'
    stop
endif

stressi(1:ntens) = stress

newton_tolerance = 1.d-8  !
max_iter = 20
! Check input
call checkinput(nprops, nstatv)
allocate(jac_inv(nstatv, nstatv))

! Update statevariables to "physical" quantities
call import_statevar(statev)

! Check elastic response
Fnew = m_2_v9(dfgrd1)
! If response is elastic, stress and ddsdde has been correctly calculated. Otherwise they are zero.
call elastic(props, statev, Fnew, yielding, stressi, ddsdde_elastic, use_el_stiff)

if (yielding) then
    !Solve for plasticity:
    
    !Set parameters for "Solve Matrix Equation"
    call sme_setparam(itertol=newton_tolerance, maxiter=max_iter)
    
    !Initiate starting guess and parameters
    call setup(x0, param, Fnew, statev)

    !Solve local problem using newton raphson iterations and analytical jacobian
    call newton_raphson_ana(residual, anajac, x0, jac_inv, lconv, info, props, statev, param, printres=0)
    !call newton_raphson_num(residual, x0, jac_inv, lconv, info, props, statev, param, printres=1)
    
    !Calculate ddsdde, statev and stress from solution of local problem
    call plastic_output(ddsddei, statev, stressi, x0, jac_inv, Fnew, props)
    
    !Check that solution converged
    call check_analysis(lconv,pnewdt)
    
endif

! Update statevars so they go back to being initially zero and scaled with dtime
call export_statevar(statev)

!If desired, return elastic tangent stiffness
if (use_el_stiff) then
    ddsddei = ddsdde_elastic
endif

ddsdde = ddsddei(1:ntens, 1:ntens)
stress = stressi(1:ntens)

END SUBROUTINE