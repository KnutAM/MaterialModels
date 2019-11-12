module GFS_module
use GFS_acegen
use tensors_module
implicit none
integer, parameter, dimension(4):: nback_allowed = (/1,2,3,4/)
double precision, parameter     :: yield_tol = 1.d-8        !Tolerance for detecting yielding, fraction of yield limit
! Parameters to be given a value
integer :: nback  !Number of backstresses
    contains
    
!Check that the input to the UMAT is correct (must be called first to assign nback)
subroutine checkinput(nprops, nstatv)
    implicit none
    integer, intent(in) ::  nprops, nstatv
    logical             ::  error
    integer             ::  k1
    
    nback = (nprops-6)/3            ! Calculate number of backstresses
    
    ! Check allowed number of backstresses 
    ! Must be first check, as error assumed to occur if a correct backstress is not found
    error = .true.
    do k1=1,size(nback_allowed)
        if (nback==nback_allowed(k1)) then
            error = .false.
        endif
    enddo
    
    ! Check correct number of material parameters
    if (nprops.ne.(6 + 3*nback)) then
        error = .true.
    endif
    
    ! Check correct number of state variables
    if (nstatv.ne.(10 + 9*nback)) then 
        error = .true.
    endif
        
    ! Quit analysis if invalid input
    if (error) then
        write(*,*) 'Invalid input to subroutine'
        write(*,"(A)", advance="no")  'nstatv = '
        write(*,*) nstatv
        write(*,"(A)", advance="no")  'nprops = '
        write(*,*) nprops
        write(*,"(A)", advance="no")  'nback = '
        write(*,*) nback
        
        call quit_analysis()
    endif
    
    
end subroutine checkinput

subroutine import_statevar(statev)
implicit none
double precision :: statev(:)
integer          :: k1

    ! Add identity to state variables
    statev(2:4) = statev(2:4) + 1.d0
    do k1=1,nback
        statev((11+9*(k1-1)):(13+9*(k1-1))) = statev((11+9*(k1-1)):(13+9*(k1-1))) + 1.d0
    enddo

end subroutine

subroutine export_statevar(statev)
implicit none
double precision :: statev(:)
integer          :: k1

    ! Add identity to state variables
    statev(2:4) = statev(2:4) - 1.d0
    do k1=1,nback
        statev((11+9*(k1-1)):(13+9*(k1-1))) = statev((11+9*(k1-1)):(13+9*(k1-1))) - 1.d0
    enddo

end subroutine   

!Calculate the elastic response, yielding logical, and if latter the stress and ddsdde variables
subroutine elastic(props, statev, Fnew, yielding, stress, ddsdde, use_el_stiff)
    implicit none
    
    double precision, intent(in)    :: props(:), statev(:), Fnew(:) ! Input variables
    double precision, intent(out)   :: stress(6), ddsdde(6,6)       ! Output stress and stiffness
    logical, intent(out)            :: yielding, use_el_stiff       ! Output logical decision variables
    double precision                :: sigma(9), dtaudF(9,9), phi       !Internal variables
    
    if (nback==1) then
        call EF1(statev, props, Fnew, sigma, dtaudF, phi)
    elseif (nback==2) then
        call EF2(statev, props, Fnew, sigma, dtaudF, phi)
        !write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (elastic)'
        !call quit_analysis()
    elseif (nback==3) then
        call EF3(statev, props, Fnew, sigma, dtaudF, phi)
        !write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (elastic)'
        !call quit_analysis()
    elseif (nback==4) then
        call EF4(statev, props, Fnew, sigma, dtaudF, phi)
        !write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (elastic)'
        !call quit_analysis()
    else
        write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (elastic)'
        call quit_analysis()
    endif
    
    
    ! Check for yielding
    if (phi>(props(3)*yield_tol)) then    !Yielding, plastic tangent stiffness used. (Tolerance scaled with initial yield limit)
        yielding = .true.
        use_el_stiff = .false.
        stress = 0.d0   !Not used
        ddsdde = 0.d0   !Not used
        
    elseif (phi>0.d0) then          !Yielding, elastic tangent stiffness used
        yielding = .true.
        use_el_stiff = .true.
        stress = 0.d0   !Not used
        ddsdde = dtaudF_2_ddsdde(dtaudF, Fnew)
    else                            !No yielding, elastic tangent stiffness used
        yielding = .false.
        use_el_stiff = .true.
        stress = sigma_v9_2_sigma_abaqus(sigma)
        ddsdde = dtaudF_2_ddsdde(dtaudF, Fnew)
    endif
    
end subroutine elastic
 
!Setup for solving the local problem in the case of plasticity
subroutine setup(x0, param, Fnew, statev)
    implicit none
    double precision, intent(out)   ::  param(:), x0(:)
    double precision, intent(in)    ::  statev(:), Fnew(9)
    
    param(1:9) = Fnew
    x0 = statev
    
end subroutine setup

!Residual function (wrapper for AceGen code)
subroutine residual(r, x, props, statev, param)
    implicit none
    double precision, intent(in)                :: x(:) 
    double precision, intent(in), optional      :: props(:), statev(:), param(:)      !Residual function parameters
    double precision, allocatable, intent(out)  :: r(:)
    double precision                            :: F(9)
    integer                                     :: nback
    nback = (size(props)-6)/3
    allocate(r(size(x)))
    
    F = param(1:9)
    
    if (nback==1) then
        call RF1(x,props,statev,F,r)
    elseif (nback==2) then
        call RF2(x,props,statev,F,r)
        !write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (RF)'
        !call quit_analysis()
    elseif (nback==3) then
        call RF3(x,props,statev,F,r)
        !write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (RF)'
        !call quit_analysis()
    elseif (nback==4) then
        call RF4(x,props,statev,F,r)
        !write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (RF)'
        !call quit_analysis()
    else
        write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (RF)'
        call quit_analysis()
    endif
    
end subroutine residual

!Jacobian function (wrapper for AceGen code)
subroutine anajac(ajacobian, x, props, statev, param)
    implicit none
    double precision, intent(in)                :: x(:)
    double precision, intent(in), optional      :: props(:), statev(:), param(:)      !parameters
    double precision, allocatable, intent(out)  :: ajacobian(:,:)
    double precision                            :: F(9)
    integer                                     :: nback
    nback = (size(props)-6)/3
    allocate(ajacobian(size(x),size(x)))
    F = param(1:9)
    
    if (nback==1) then
        call dRdX1(x,props,statev,F,ajacobian)
    elseif (nback==2) then
        call dRdX2(x,props,statev,F,ajacobian)
        !write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (RFJ)'
        !call quit_analysis()
    elseif (nback==3) then
        call dRdX3(x,props,statev,F,ajacobian)
        !write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (RFJ)'
        !call quit_analysis()
    elseif (nback==4) then
        call dRdX4(x,props,statev,F,ajacobian)
        !write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (RFJ)'
        !call quit_analysis()
    else
        write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (RFJ)'
        call quit_analysis()
    endif
    
end subroutine anajac

!Given the solution to the local problem, calculate updated state variables and tangent stiffness
subroutine plastic_output(ddsdde, statev, stress, xsol, invjac, F, props)
    implicit none
    !Interface parameters
    double precision, intent(inout) ::  ddsdde(6,6), statev(:), stress(6)
    double precision, intent(in)    ::  xsol(:), invjac(:,:), F(9), props(:)
    !Internal parameters
    double precision, dimension(9)  ::  sigma, Fp
    double precision, dimension(9,9)::  dtau_dF, pdtau_pdF, pdtau_pdFp, dFp_dF
    double precision                ::  dRdF_X(size(xsol),9), eprops(2)
    integer                         ::  k1
    
    Fp = xsol(2:10)     !Plastic def gradient
    eprops = props(1:2) !Elastic properties
    call EoutF(Fp, eprops, F, sigma, pdtau_pdF, pdtau_pdFp)
    
    if (nback==1) then
        call dRdF1(xsol,props,statev,F,dRdF_X)
    elseif (nback==2) then
        call dRdF2(xsol,props,statev,F,dRdF_X)
        !write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (plastic_output)'
        !call quit_analysis()
    elseif (nback==3) then
        call dRdF3(xsol,props,statev,F,dRdF_X)
        !write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (plastic_output)'
        !call quit_analysis()
    elseif (nback==4) then
        call dRdF4(xsol,props,statev,F,dRdF_X)
        !write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (plastic_output)'
        !call quit_analysis()
    else
        write(*,"(A, I1, A)") 'nback = ', nback, ' is not supported (plastic_output)'
        call quit_analysis()
    endif
    
    dFp_dF  = -matmul(invjac(2:10,:), dRdF_X)
    dtau_dF = pdtau_pdF + matmul(pdtau_pdFp, dFp_dF)
    
    ddsdde = dtaudF_2_ddsdde(dtau_dF, F)
    stress = sigma_v9_2_sigma_abaqus(sigma)
    statev = xsol
    
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

!Safely quit analysis in a correct way
subroutine quit_analysis()  
    implicit none
    !call xit()      !Quit abaqus analysis in a correct way (Comment out if running outside abaqus)
    stop            !Use if running outside abaqus
end subroutine quit_analysis

end MODULE