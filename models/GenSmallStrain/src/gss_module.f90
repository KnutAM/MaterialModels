module gss_module
use acegen_mod
use model_module
use SolveMatrixEquation
implicit none

   !double precision, parameter, dimension(6):: I2 = (/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/)
    double precision, parameter :: p23 = 2.d0/3.d0
    double precision, parameter :: m13 = -1.d0/3.d0
    double precision, parameter :: p12 = 1.d0/2.d0
    double precision, parameter, dimension(6):: Z6 = [0.d0,0.d0,0.d0,0.d0,0.d0,0.d0]
    double precision, parameter, dimension(6):: I2 = [1.d0,1.d0,1.d0,0.d0,0.d0,0.d0]
    double precision, parameter, dimension(6,6):: IoI = reshape([I2,I2,I2,Z6,Z6,Z6], [6,6])
    double precision, parameter, dimension(6,6):: I4dev = reshape([ p23, m13, m13, 0.0d0,0.0d0,0.0d0 ,&
                                                                    m13, p23, m13, 0.0d0,0.0d0,0.0d0 ,&
                                                                    m13, m13, p23, 0.0d0,0.0d0,0.0d0 ,&
                                                                    0.d0,0.d0,0.d0,0.5d0,0.0d0,0.0d0 ,&
                                                                    0.d0,0.d0,0.d0,0.0d0,0.5d0,0.0d0 ,&
                                                                    0.d0,0.d0,0.d0,0.0d0,0.0d0,0.5d0 ], &
                                                                  [6,6])
    
    contains
    
!Calculate the elastic response, yielding logical, and if latter the stress and ddsdde variables
subroutine elastic(mpar, statev, deps, stress_old, stress, ddsdde, yielding)
    implicit none
    
    double precision, intent(in)    :: mpar(:), statev(:), deps(:), stress_old(:) ! Input variables
    double precision, intent(inout) :: stress(6)    ! Output trial stress
    double precision, intent(out)   :: ddsdde(6,6)  ! Output stiffness
    logical, intent(out)            :: yielding     ! Output if yielding
    ! Internal variables
    double precision                :: phi          ! Yield function (<0 if elastic)
    double precision                :: beta(6)      ! Total back-stress
    double precision                :: sred(6)      ! Reduced stress
    double precision                :: kappa        ! Isotropic hardening stress
    integer                         :: k1           ! Looping iterator
    double precision                :: G, K, la, Y0 ! Material parameters
    !double precision                :: IoI(6,6), I4dev(6,6)
    
    ! Get material parameters
    G = mpar(1)
    K = mpar(2)
    la = K - (2.d0/3.d0)*G
    Y0 = mpar(3)
    
    ! Calculate the total back-stress, beta
    beta = 0.d0
    do k1=1,nback
        beta = beta + statev((3+nvar_add+(k1-1)*6):(1+nvar_add+k1*6))
    enddo
    kappa = statev(1)
    
    
!    stress = stress_old + K*sum(deps(1:3))*I2 + 2*G*matmul(I4dev, deps)
    stress(1:3) = stress_old(1:3) + la*sum(deps(1:3)) + 2*G*deps(1:3)
    stress(4:6) = stress_old(4:6) + G*deps(4:6)
    
    sred = stress-beta
    phi = yield_function(sred, Y0+kappa)
    
    yielding = (phi>0)
    
    ! Calculate elastic stiffness if not above yield limit
    if (.not.yielding) then
        ddsdde = K*IoI + 2*G*I4dev
    endif
    
end subroutine elastic
 
!Setup for solving the local problem in the case of plasticity
subroutine setup(x0, xold, param, stress_old, stress_trial, statev, dtime)
    implicit none
    double precision, allocatable   ::  param(:), x0(:), xold(:)
    double precision, intent(in)    ::  stress_old(6), stress_trial(6), statev(:), dtime
    
    allocate(x0(nvar), xold(nvar), param(7))
    
    xold(1:6)     = get_deviatoric_part(stress_old)
    xold(7:nvar)  = statev
    
    x0 = xold
    param(1:6) = get_deviatoric_part(stress_trial)
    param(7) = dtime
    
end subroutine setup

!Residual function (wrapper for AceGen code)
subroutine residual(r, x, props, xold, param)
    implicit none
    double precision, intent(in)                :: x(:) 
    double precision, intent(in), optional      :: props(:), xold(:), param(:)      !Residual function parameters
    double precision, allocatable, intent(out)  :: r(:)

    allocate(r(size(x)))
    
    if (nback==1) then
        call RF1(x,props,xold,param(1:6),param(7),r)
    elseif (nback==2) then
        call RF2(x,props,xold,param(1:6),param(7),r)
    elseif (nback==3) then
        call RF3(x,props,xold,param(1:6),param(7),r)
    elseif (nback==4) then
        call RF4(x,props,xold,param(1:6),param(7),r)
    else
        write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (RF)'
        stop
    endif
    
end subroutine residual

!Jacobian function (wrapper for AceGen code)
subroutine anajac(ajacobian, x, props, xold, param)
    implicit none
    double precision, intent(in)                :: x(:)
    double precision, intent(in), optional      :: props(:), xold(:), param(:)      !parameters
    double precision, allocatable, intent(out)  :: ajacobian(:,:)
    
    
    allocate(ajacobian(size(x),size(x)))
    
    if (nback==1) then
        call dRdX1(x,props,xold,param(1:6),param(7),ajacobian)
    elseif (nback==2) then
        call dRdX2(x,props,xold,param(1:6),param(7),ajacobian)
    elseif (nback==3) then
        call dRdX3(x,props,xold,param(1:6),param(7),ajacobian)
    elseif (nback==4) then
        call dRdX4(x,props,xold,param(1:6),param(7),ajacobian)
    else
        write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (RFJ)'
        stop
    endif
    
end subroutine anajac

!Given the solution to the local problem, calculate updated state variables and tangent stiffness
subroutine plastic_output(ddsdde, statev, stress, xsol, invjac, mpar, stress_old, stran, dstran)
    implicit none
    !Interface parameters
    double precision, intent(inout) ::  ddsdde(6,6), statev(:), stress(6)
    double precision, intent(in)    ::  xsol(:), invjac(:,:), mpar(:), stran(6), stress_old(6), dstran(6)
    !Internal parameters
    double precision                ::  K, G
    !double precision                ::  I4dev(6,6), IoI(6,6)
    
    !call identities(IoI, I4dev)
    G = mpar(1)
    K = mpar(2)
    
    stress = get_deviatoric_part(xsol(1:6)) + K*sum(stran(1:3)+dstran(1:3))*I2
    
    ddsdde = -2*G*I4dev  ! dR/deps at constant x
    ddsdde = -matmul(invjac(1:6, 1:6), ddsdde)  ! dsigma_dev/deps
    ddsdde = K*IoI + ddsdde ! dsigma/deps
    
    statev = xsol(7:nvar)
    
end subroutine plastic_output

! Check if problem converged
subroutine check_analysis(lconv, pnewdt)
    implicit none
    logical, intent(in) :: lconv
    double precision    :: pnewdt
    if (.not.lconv) then
        write(*,*) 'Local problem did not converge'
        pnewdt = 5.d-1
        return
    endif
end subroutine check_analysis

function yield_function(sred, yrad) result(phi)
implicit none
    double precision    :: sred(6), yrad
    double precision    :: phi
    double precision    :: sdev(6)
    
    sdev = get_deviatoric_part(sred)
    
    phi = sqrt(3.d0/2.d0)*sqrt(sum(sdev*sdev) + sum(sdev(4:6)*sdev(4:6))) - yrad
    
end function

function get_deviatoric_part(voigt6) result(dev)
implicit none
    double precision    :: voigt6(6)
    double precision    :: dev(6)
    
    dev = voigt6 - I2*sum(voigt6(1:3))/3.d0
    
end function
    
end MODULE
