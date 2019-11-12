! Material parameter definition for the Delobelle model, nprops = 6+2*nback
! Input:    props = E, v, Y0, Hiso, invYiso, delta, Hk1, invYk1, Hk2, ...
! Internal: mpar  = G, K, Y0, Hiso, invYiso, delta, Hk1, invYk1, {mk1}, Hk2, ...
    
module model_module
use acegen_mod
implicit none
    
    integer, parameter, dimension(4):: nback_allowed = (/1,2,3,4/)
    integer, parameter              :: nvar_add = 0 !Addition of variables relative rate independent model (hence zero here)
    integer                         :: nback    !Number of backstresses
    integer                         :: nvar     !Number of variables in local problem
    
    contains

!Check that the input to the UMAT is correct (must be called first to assign nback)
subroutine checkinput(nprops, nstatv, props, mpar)
    implicit none
    integer, intent(in)             :: nprops, nstatv
    double precision, intent(in)    :: props(:) ! Material properties input to umat
    double precision, allocatable   :: mpar(:)  ! Internally used material properties, to be consistent for different models
    logical                         :: error
    integer                         :: k1
    
    nback = (nstatv-2)/6            ! Calculate number of backstresses
    nvar = nstatv + 6               ! Calculate number of variables in local problem
    
    error = .false.
    ! Check allowed number of backstresses 
    if (.not.any(nback==nback_allowed)) then
        error = .true.
    endif
    
    ! Check correct number of material parameters
    if (nprops.ne.(6 + 2*nback)) then
        error = .true.
    endif
    
    ! Check correct number of state variables
    if (nstatv.ne.(2 + 6*nback)) then 
        error = .true.
    endif
        
    ! Quit analysis if invalid input
    if (error) then
        write(*,*) 'Invalid input to the Delobelle model (umat), it should be'
        write(*,*) 'nprops = 6+2*nback and'
        write(*,*) 'nstatv = 2+6*nback, but'
        write(*,"(A)", advance="no")  'nstatv = '
        write(*,*) nstatv
        write(*,"(A)", advance="no")  'nprops = '
        write(*,*) nprops
        write(*,"(A)", advance="no")  'nback = '
        write(*,*) nback
        stop
    endif
    
    ! Put material parameters at consistent location to allow using the same AceGen code for multiple models:
    
    allocate(mpar(6+3*nback))
    mpar = 0.d0
    mpar(1) = props(1)/(2*(1+props(2)))     ! Convert from E and v to G = E/(2*(1+v))
    mpar(2) = props(1)/(3*(1-2*props(2)))   ! Convert from E and v to K = E/(3*(1-2v))
    mpar(3:6) = props(3:6)
    do k1 = 1,nback
        mpar(7+3*(k1-1)) = props(7+2*(k1-1))
        mpar(8+3*(k1-1)) = props(8+2*(k1-1))
    enddo

end subroutine checkinput

end module
