module umat_module
use acegen_mod
implicit none

integer :: nvar         ! Number of variables in local problem

    contains
    
!Check that the input to the UMAT is correct
subroutine checkinput(nprops, nstatv)
    implicit none
    integer, intent(in) :: nprops, nstatv
    logical             :: error
    integer             :: nstatv_model ! Number of state variables in model
    integer             :: nprops_model ! Number of material params in model
    
    call model_size(nprops_model,nstatv_model,nvar)
    
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

!Residual function (wrapper for AceGen code)
subroutine residualfun(r, x, props, statev, param)
    implicit none
    double precision, intent(in)                :: x(:) 
    double precision, intent(in), optional      :: props(:), statev(:), param(:)      !Residual function parameters
    double precision, allocatable, intent(out)  :: r(:)
    
    allocate(r(nvar))
    
    call residual(x,props,statev,param,r)
    
end subroutine residualfun

!Jacobian function (wrapper for AceGen code)
subroutine anajac(ajacobian, x, props, statev, param)
    implicit none
    double precision, intent(in)                :: x(:)
    double precision, intent(in), optional      :: props(:), statev(:), param(:)      !parameters
    double precision, allocatable, intent(out)  :: ajacobian(:,:)
    
    allocate(ajacobian(nvar,nvar))
    
    call jacobian(x,props,statev,param,ajacobian)
        
end subroutine anajac


! Check if problem converged
subroutine check_analysis(lconv, pnewdt)
    implicit none
    logical, intent(in) :: lconv
    double precision    :: pnewdt
    if (.not.lconv) then
        !write(*,*) 'Local problem did not converge'
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

end module