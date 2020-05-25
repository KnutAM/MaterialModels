module umat_module
use acegen_crystal_plasticity_mod
implicit none

integer :: nvar         ! Number of variables in local problem
integer :: nslip        ! Number of slip systems
integer :: nstatv_model ! Number of state variables in model

    contains
    
!Check that the input to the UMAT is correct
subroutine checkinput(nprops, nstatv)
    implicit none
    integer, intent(in) :: nprops, nstatv
    integer             :: nprops_model
    logical             :: error
    
    call model_size(nprops_model,nstatv_model,nvar,nslip)
    
    error = .false.
    if (nprops.ne.nprops_model) then
        error = .true.
    elseif (nstatv.ne.nstatv_model) then
        error = .true.
    endif
        
    ! Quit analysis if invalid input
    if (error) then
        write(*,*) 'Invalid input to umat subroutine'
        write(*,*) 'You specified ', nprops, ' material parameters: ', nprops_model, ' should be given'
        write(*,*) 'You specified ', nstatv, ' state variables: ', nstatv_model, ' should be given'
        call quit_analysis()
    endif
    
    
end subroutine checkinput

!Setup for solving the local problem in the case of plasticity
subroutine setup(x0, param, Fnew, TauRedTr, dtime, statev)
    implicit none
    double precision, allocatable   ::  param(:)
    double precision, allocatable   ::  x0(:)
    double precision, intent(in)    ::  TauRedTr(:), statev(:), Fnew(9), dtime
    allocate(param(9+nslip+1))
    allocate(x0(nvar))
    param(1:9) = Fnew
    param(10:(9+nslip)) = TauRedTr
    param(9+nslip+1) = dtime
    x0 = statev((nstatv_model-nvar+1):nstatv_model)
    
end subroutine setup

!Residual function (wrapper for AceGen code)
subroutine residualfun(r, x, props, statev, param)
    implicit none
    double precision, intent(in)                :: x(:) 
    double precision, intent(in), optional      :: props(:), statev(:), param(:)      !Residual function parameters
    double precision, allocatable, intent(out)  :: r(:)
    double precision                            :: Fnew(9), TauRedTr(nslip), dtime
    
    allocate(r(nvar))
    
    Fnew = param(1:9)
    TauRedTr = param(10:(9+nslip))
    dtime = param(9+nslip+1)
    
    call residual(x,props,statev,Fnew,TauRedTr,dtime,r)
    
end subroutine residualfun

!Jacobian function (wrapper for AceGen code)
subroutine anajac(ajacobian, x, props, statev, param)
    implicit none
    double precision, intent(in)                :: x(:)
    double precision, intent(in), optional      :: props(:), statev(:), param(:)      !parameters
    double precision, allocatable, intent(out)  :: ajacobian(:,:)
    double precision                            :: Fnew(9), TauRedTr(nslip), dtime
    
    allocate(ajacobian(nvar,nvar))
    
    Fnew = param(1:9)
    TauRedTr = param(10:(9+nslip))
    dtime = param(9+nslip+1)
    
    call jacobian(x,props,statev,Fnew,TauRedTr,dtime,ajacobian)
        
end subroutine anajac


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

function m3x3_2_v9(mat) result(v9)
    implicit none
    double precision    :: mat(3,3), v9(9)
    
    v9 = (/mat(1,1), mat(2,2), mat(3,3), mat(1,2), mat(2,3), mat(3,1), mat(1,3), mat(2,1), mat(3,2)/)
    
end function

!Safely quit analysis in a correct way
subroutine quit_analysis()  
    implicit none
    !call xit()      !Quit abaqus analysis in a correct way (Comment out if running outside abaqus)
    stop            !Use if running outside abaqus
end subroutine quit_analysis

end MODULE