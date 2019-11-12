! SolveMatrixEquation   Routines for solution of R(X)=0
!
! Main subroutines
!   sme_setparam        Set parameters for solving of matrix equations. Required to call, but can be called without parameters
!   numjac              Calculate the numerical jacobian
!   solve_step          Solve a step, i.e. solve for dx so that A(x0+dx)=0
!   invert_jacobian     Calculate the inverted jacobian from the lu-factorized coming from solve_step
!   newton_raphson_num  Solve the equation system R(X)=0 using numerical jacobian
!   newton_raphson_ana  Solve the equation system R(X)=0 using analytical jacobian
! 
! --------------------------------------------------------------------------------------------------------------
! Test subroutines (Use for testing and example)
!   test_residual       An example residual function, 3 unknowns
!   test_anajac         Analytical jacobian for test_residual
! End of SolveMatrixEquation description

module SolveMatrixEquation
implicit none
! Default parameters
integer, parameter          :: DEF_maxiter = 10         !Maximum number of iterations for newton-raphson iterations
double precision, parameter :: DEF_itertol = 1.d-6      !Tolerance on norm(residual) for newton-raphson iterations
double precision, parameter :: DEF_numdeps = 1.d-8      !Pertubation for numerical jacobian
logical, parameter          :: DEF_symnumjac = .false.  !Set true to use symmetrical calculation of jacobian (slower)
integer, parameter          :: DEF_line_search_start = 0!Default to start imediately
double precision, parameter :: DEF_alpha_min = 1.d-3     ! Default is to take minimum 1/100 of step size
double precision, parameter :: DEF_alpha_mult = 5.d-1     ! Default is to take half step at every new attempt

double precision, parameter :: infty = huge(1.d0)       !Set value to compare to check for infinity

! Define actual parameters
integer          :: SME_maxiter     !Maximum number of iterations for newton-raphson iterations
double precision :: SME_itertol     !Tolerance on norm(residual) for newton-raphson iterations
double precision :: SME_numdeps     !Pertubation for numerical jacobian
logical          :: SME_symnumjac   !Set true to use symmetrical calculation of jacobian (slower)
integer          :: SME_line_search_start   ! From what iteration to initiate line search if residual increases
double precision :: SME_alpha_min           ! Minimum value of alpha before failing
double precision :: SME_alpha_mult          ! Multiplication factor for reducing alpha

    contains
!Setup parameters for SolveMatrixEquation
subroutine sme_setparam(maxiter, itertol, pertubation, symnumjac, line_search_start, alpha_min, alpha_mult)
    !Set parameters for SolveMatrixEquation module (see specification statements of module)
    !Note: All options are optional, if not specified default values will be used
    !To call routine, use call sme_setparam(maxiter=NN, itertol=nn.d-n, pertubation=nn.d-n, symnumjac=.true.), 
    !noting that you don't need to specify all
    implicit none
    integer, optional           :: maxiter
    double precision, optional  :: itertol
    double precision, optional  :: pertubation
    logical, optional           :: symnumjac
    integer, optional           :: line_search_start
    double precision, optional  :: alpha_min
    double precision, optional  :: alpha_mult
    
    SME_maxiter             = DEF_maxiter
    SME_itertol             = DEF_itertol
    SME_numdeps             = DEF_numdeps
    SME_symnumjac           = DEF_symnumjac
    SME_line_search_start   = DEF_line_search_start
    SME_alpha_min           = DEF_alpha_min
    SME_alpha_mult          = DEF_alpha_mult
    
    !Set maxiter
    if (present(maxiter)) then
        SME_maxiter = maxiter
    endif
    
    !Set itertol
    if (present(itertol)) then
        SME_itertol = itertol
    endif
    
    !Set pertubation (numdeps)
    if (present(pertubation)) then
        SME_numdeps = pertubation
    endif
    
    !Set condition for symmetric jacobian
    if (present(symnumjac)) then
        SME_symnumjac = symnumjac
    endif
    
    ! Line search start
    if (present(line_search_start)) then
        SME_line_search_start = line_search_start
    endif
    
    ! Line search minimum alpha value
    if (present(alpha_min)) then
        if (alpha_min<=0.d0) then  ! Check needed to prevent infinite loop
            write(*,*) 'alpha_min must be between 0 and 1 (strictly larger than 0!'
            stop
        endif
        
        SME_alpha_min = alpha_min
    endif
    
    ! Line search alpha reduction factor
    if (present(alpha_mult)) then
        if (alpha_mult>=1.d0) then  ! Check needed to prevent infinite loop
            write(*,*) 'alpha_mult must be strictly less than 1.0'
            stop
        endif
        SME_alpha_mult = alpha_mult
    endif
    
end subroutine sme_setparam

!Calculate numerical jacobian
subroutine numjac(residual_function, x0, jacobian, props, statev, param)
    !residual: Residual subroutine
    !x0      : Input variable to residual subroutine where the jacobian should be calculated
    !param   : Set of parameters used in the residual subroutine
    !jacobian: Output from subroutine, giving the numerical jacobian of the residual function at x0
    
    implicit none
    double precision, intent(in)                :: x0(:)
    double precision, optional, intent(in)      :: props(:), statev(:), param(:)!Input parameters
    double precision, allocatable, intent(out)  :: jacobian(:,:)                !Output arguments
    double precision, allocatable               :: dx(:), r0(:), rc(:)          !Internal allocatable float variables
    integer                                     :: nvar, k1, ninput             !Internal integers
    
    interface   !Interface for residual function
        subroutine residual_function(r, x, props, statev, param)
            double precision, intent(in)                :: x(:)
            double precision, optional, intent(in)      :: props(:), statev(:), param(:)      !Input arguments
            double precision, allocatable, intent(out)  :: r(:)
        end subroutine residual_function
    end interface
    
    nvar  = size(x0)  ! Size of problem
    
    !Allocate internal variables
    allocate(dx(nvar), r0(nvar), rc(nvar), jacobian(nvar,nvar))
    
    if (SME_symnumjac) then 
        do k1=1,nvar
            dx = 0.d0
            dx(k1) = SME_numdeps/2.d0
            call residual_function(rc, x0+dx, props, statev, param)    !Get residual rc at x0+dx
            call residual_function(r0, x0-dx, props, statev, param)    !Get residual r0 at x0-dx
            jacobian(:,k1) = (rc-r0)/SME_numdeps   !Calculate k1 column in jacobian
        end do
    else
        call residual_function(r0, x0, props, statev, param)    !Get initial residual r0
        do k1=1,nvar
            dx = 0.d0
            dx(k1) = SME_numdeps
            call residual_function(rc, x0+dx, props, statev, param) !Get current residual rc
            jacobian(:,k1) = (rc-r0)/SME_numdeps   !Calculate k1 column in jacobian
        end do
    endif
    
end subroutine numjac

!Solve one step in newton iteration
subroutine solve_step(jacobian, r0, dx, lconv, ipiv, info)
! Solve jacobian*dx+r0=0 => jacobian*dx = -r0
! Note: jacobian is returned as lu-factorized, as well as pivot indicies ipiv
! usable in dgetri directly instead of calling dgetrf first
! dgesv documentation: https://software.intel.com/en-us/node/468876
    implicit none
    double precision, intent(inout) :: jacobian(:,:)    !Outputted ad LU-factorized
    double precision, intent(in)    :: r0(:)            !Residual at current x0
    double precision, intent(out)   :: dx(:)            !Suggested step
    logical, intent(out)            :: lconv            !Tells if solution succeeded
    integer, intent(out)            :: info             !Information about solution
    integer                         :: nvar             !Number of variables
    integer,intent(out)             :: ipiv(:)          !pivot indices for LU decomposition
    
    dx = -r0    !Input b in dgesv is overwritten to answer x
    nvar = size(r0)

    !    dgesv( n,   nrhs,        a,  lda, ipiv,  b,  ldb, info )
    call dgesv(nvar,    1, jacobian, nvar, ipiv, dx, nvar, info)
    lconv = (info==0)
end subroutine solve_step

!Invert matrix (typically jacobian)
subroutine invert_jacobian(lu_jacobian, ipiv)
    !lu_jacobian is the jacobian lu-factorized using either solve_step, or dgesv/dgetrf directly.
    !ipiv (indicies of pivot elements) are also output from solve_step, dgesv and dgetrf required
    !lu_jacobian is modified in this subroutine to become the inverse of the jacobian matrix
    implicit none
    double precision, intent(in) :: lu_jacobian(:,:)
    integer                      :: nvar, lwork, info
    integer, intent(out)         :: ipiv(:)
    double precision, allocatable:: work(:)
    nvar = size(lu_jacobian, dim=1)
    lwork = 2*nvar  !2 to be safe, to optimize see dgetri documentation
    !Probably not critical, as this routine is normally only called once per material execution
    
    allocate(work(lwork))
    !    dgetri(   n,           a,  lda, ipiv, work, lwork, info )
    call dgetri(nvar, lu_jacobian, nvar, ipiv, work, lwork, info)
    
end subroutine invert_jacobian

!Solve equation system using numerical jacobian
subroutine newton_raphson_num(residual_function, x0, jac_inv, lconv, info, props, statev, param, rvalues, xvalues, printres)
!Main routine for use to solve equation system with numerical jacobian:
!Specify the residual_function (by giving the actual name of the subroutine defining the residual function for the local problem)
!The inverse jacobian will be given as output for further use. lconv is true if convergence achieved, otherwise it is false 
implicit none
    double precision,allocatable,intent(out):: jac_inv(:,:)         !Inverted jacobian at stationary point returned (different uses during subroutine)
    double precision, intent(inout)         :: x0(:)                !Start guess
    double precision, intent(in), optional  :: props(:), statev(:), param(:)    !Residual function parameters
    integer,          intent(in), optional  :: printres             !Logical if print residual or not
    logical,          intent(out)           :: lconv                !Logical if convergence or not
    integer                                 :: k1, k2               !Looping iterator
    integer                                 :: ipiv(size(x0))       !Indicies of pivot elements
    double precision                        :: dx(size(x0))         !Increment
    double precision, allocatable           :: r0(:)                !Residual
    integer                                 :: info                 !Information about solution of step (dgesv)
    logical                                 :: diagnostic           !Logical to describe if diagnostic information should be given or not
    integer                                 :: printinfo            !Integer to decide what info to print
    double precision, optional, allocatable :: rvalues(:,:)         !Optional diagnostics information - residual values
    double precision, optional, allocatable :: xvalues(:,:)         !Optional diagnostics information - x-values
    double precision, allocatable           :: rvalues_full(:,:)    !Used to fill calculate rvalues
    double precision, allocatable           :: xvalues_full(:,:)    !Used to fill calculate xvalues
    character(len=20)                       :: jacobian_format      !Used to print jacobian nicer
    
    interface   !interface for residual function
        subroutine residual_function(r, x, props, statev, param)
            double precision, intent(in)                :: x(:)
            double precision, intent(in), optional      :: props(:), statev(:), param(:)      !Residual function parameters
            double precision, allocatable, intent(out)  :: r(:)
        end subroutine residual_function
    end interface
    
    diagnostic = present(rvalues).or.present(xvalues)   !Check if diagnostic info should be saved
    ! If diagnostic, then allocate the full 2d arrays for saving information
    if (diagnostic) then
        allocate(rvalues_full(size(x0), SME_maxiter+1), xvalues_full(size(x0), SME_maxiter+1))
    endif
    
    ! Should residual norm be printed out for each iterations?
    if (present(printres)) then
        printinfo = printres       !Same if present
    else
        printinfo = 0              !Don't print if not specified
    endif
    
    call residual_function(r0, x0, props, statev, param)
    
    ! If diagnostic, then save the current state
    if (diagnostic) then
        rvalues_full(:, 1) = r0
        xvalues_full(:, 1) = x0
    endif
    
    ! If print residual, print
    if (printinfo==1) then
        write(*,*) 'Printing iteration nr and residual values:'
        write(*,"(I2,E15.5)") 0, sqrt(sum(r0*r0))
    elseif (printinfo==2) then
        write(*,*) 'Printing full variable and residual'
        write(*,*) 'x'
        write(*,*) x0
        write(*,"(A,E15.5)") 'r(x): ', sqrt(sum(r0*r0))
        write(*,*) r0
        write(*,*) ' '
    elseif (printinfo==3) then
        write(*,*) 'Printing full variable, residual and jacobian (transpose)'
        write(jacobian_format,"(I0)") size(r0)
        jacobian_format = '('//trim(adjustl(jacobian_format))//'E15.5)'
    endif
    
    do k1 = 1,SME_maxiter
        call numjac(residual_function, x0, jac_inv, props, statev, param)  !Get jacobian
        if (badcontent2(jac_inv)) then
            write(*,*) 'Jacobian NaN or Inf. Printing x'
            write(*,*) x0
            lconv = .false.
            return
        endif
        if (printinfo==3) then
            write(*,"(A,I2)") 'x_', k1
            write(*,*) x0
            write(*,"(A,E15.5)") 'r(x): ', sqrt(sum(r0*r0))
            write(*,*) r0
            write(*,"(A)") 'Jacobian: '
            do k2=1,size(r0)
                write(*,jacobian_format) jac_inv(:,k2)
            enddo
            write(*,*) ' '
        endif
        call solve_step(jac_inv, r0, dx, lconv, ipiv, info) !Get step increment dx
        x0 = x0 + dx                                        !Update step
        call residual_function(r0, x0, props, statev, param)               !Get new residual
        if (badcontent1(r0)) then
            write(*,*) 'Residual NaN or Inf. Printing x and r'
            write(*,*) x0
            write(*,*) r0
            lconv = .false.
            return
        endif
        ! If diagnostic, then save the current state
        if (diagnostic) then
            rvalues_full(:, 1+k1) = r0
            xvalues_full(:, 1+k1) = x0
        endif
        
        ! If print residual, print        
        if (printinfo==1) then
            write(*,"(I2,E15.5)") k1, sqrt(sum(r0*r0))
        elseif (printinfo==2) then
            write(*,"(A,I2)") 'x_', k1
            write(*,*) x0
            write(*,"(A,E15.5)") 'r(x): ', sqrt(sum(r0*r0))
            write(*,*) r0
            write(*,*) ' '
        endif
        
        if (sqrt(sum(r0*r0))<SME_itertol) then              !Check for convergence
            lconv = .true.                                  !Solution converged
            call invert_jacobian(jac_inv, ipiv)             !Calculate inverse jacobian
            
            ! If diagnostic, allocate and save the diagnostic information to rvalues and xvalues as requested
            if (diagnostic) then
                if (present(rvalues)) then
                    allocate(rvalues(size(x0), 1+k1))
                    rvalues = rvalues_full(:, 1:(1+k1))
                endif
                if (present(xvalues)) then
                    allocate(xvalues(size(x0), 1+k1))
                    xvalues = xvalues_full(:, 1:(1+k1))
                endif
            endif
            
            return                                          !End subroutine when convergence achieved
        endif
    enddo
    lconv = .false. !If subroutine didn't end from convergence, r0>SME_itertol

end subroutine newton_raphson_num

!Solve equation system using analytical jacobian
subroutine newton_raphson_ana(residual_function, jacobian_function, x0, jac_inv, lconv, info, props, statev, param, rvalues, xvalues, printres)
!Main routine for use to solve equation system with analytical jacobian:
!Specify the residual_function (by giving the actual name of the subroutine defining the residual function for the local problem)
!Specify the jacobian_function (by giving the actual name of the subroutine defining the analytical jacobian for the local problem)
!The inverse jacobian will be given as output for further use. lconv is true if convergence achieved, otherwise it is false 
implicit none
    double precision,allocatable,intent(out):: jac_inv(:,:)  !Inverted jacobian at stationary point returned (different uses during subroutine)
    double precision, intent(inout)         :: x0(:)         !Start guess
    double precision, intent(in), optional  :: props(:), statev(:), param(:)      !Residual function parameters
    integer,          intent(in), optional  :: printres             !Logical if print residual or not
    logical,          intent(out)           :: lconv         !Logical if convergence or not
    integer                                 :: k1, k2            !Looping iterator
    integer                                 :: ipiv(size(x0))!Indicies of pivot elements
    double precision                        :: dx(size(x0))  !Increment
    double precision, allocatable           :: r0(:)         !Residual
    integer                                 :: info          !Information about solution of step (dgesv)
    logical                                 :: diagnostic    !Logical to describe if diagnostic information should be given or not
    integer                                 :: printinfo    
    double precision, optional, allocatable :: rvalues(:,:)  !Optional diagnostics information - residual values
    double precision, optional, allocatable :: xvalues(:,:)  !Optional diagnostics information - x-values
    double precision, allocatable           :: rvalues_full(:,:)       !Used to fill calculate rvalues
    double precision, allocatable           :: xvalues_full(:,:)       !Used to fill calculate xvalues
    double precision                        :: residual, residual_old   ! Current and previous residual (scalar = sqrt(sum(r0**2)) )
    double precision                        :: alpha
    character(len=20)                       :: jacobian_format      !Used to print jacobian nicer
    
    interface   !Interface for residual function
        subroutine residual_function(r, x, props, statev, param)
            double precision, intent(in)                :: x(:)
            double precision, intent(in), optional      :: props(:), statev(:), param(:)      !Residual function parameters
            double precision, allocatable, intent(out)  :: r(:)
        end subroutine residual_function
    end interface
    
    interface   !Interface for jacobian function
        subroutine jacobian_function(j, x, props, statev, param)
            double precision, intent(in)                :: x(:)
            double precision, intent(in), optional      :: props(:), statev(:), param(:)      !Residual function parameters
            double precision, allocatable, intent(out)  :: j(:,:)
        end subroutine jacobian_function
    end interface
    diagnostic = present(rvalues).or.present(xvalues)   !Check if diagnostic info should be saved
    
    ! Should residual norm be printed out for each iterations?
    if (present(printres)) then
        printinfo = printres       !Same if present
    else
        printinfo = 0              !Don't print if not specified
    endif
    
    ! If diagnostic, then allocate the full 2d arrays for saving information
    if (diagnostic) then
        allocate(rvalues_full(size(x0), SME_maxiter+1), xvalues_full(size(x0), SME_maxiter+1))
    endif
    
    !Get initial residual
    call residual_function(r0, x0, props, statev, param)
    residual_old = sqrt(sum(r0**2))
    
    ! If diagnostic, then save the current state
    if (diagnostic) then
        rvalues_full(:, 1) = r0
        xvalues_full(:, 1) = x0
    endif
    
    ! If print residual, print
    if (printinfo==1) then
        write(*,*) 'Printing iteration nr and residual values:'
        write(*,"(I2,E15.5)") 0, residual_old
    elseif (printinfo==2) then
        write(*,*) 'Printing full variable and residual'
        write(*,*) 'x'
        write(*,*) x0
        write(*,"(A,E15.5)") 'r(x): ', residual_old
        write(*,*) r0
        write(*,*) ' '
    elseif (printinfo==3) then
        write(*,*) 'Printing full variable, residual and jacobian (transpose)'
        write(jacobian_format,"(I0)") size(r0)
        jacobian_format = '('//trim(adjustl(jacobian_format))//'E15.5)'
    endif
        
    do k1 = 1,SME_maxiter
        call jacobian_function(jac_inv, x0, props, statev, param)      !Get jacobian
        if (diagnostic) then
            if (badcontent2(jac_inv)) then  ! Can take some time, so only check when in diagnostic mode
                write(*,*) 'Jacobian NaN or Inf. Printing x'
                write(*,*) x0
                lconv = .false.
                return
            endif
        endif
        if (printinfo==3) then
            write(*,"(A,I2)") 'x_', k1
            write(*,*) x0
            write(*,"(A,E15.5)") 'r(x): ', sqrt(sum(r0*r0))
            write(*,*) r0
            write(*,"(A)") 'Jacobian: '
            do k2=1,size(r0)
                write(*,jacobian_format) jac_inv(:,k2)
            enddo
            write(*,*) ' '
        endif
        call solve_step(jac_inv, r0, dx, lconv, ipiv, info) !Get step increment dx
        x0 = x0 + dx                                        !Update step
        
        call residual_function(r0, x0, props, statev, param)               !Get new residual
        
        residual = sqrt(sum(r0**2))
        if (badcontent0(residual)) then
            write(*,*) 'Residual NaN or Inf. Printing x and r'
            write(*,*) x0
            write(*,*) r0
            lconv = .false.
            return
        endif
        
        alpha = 1.d0
        do while ( (k1>SME_line_search_start).and.(residual>residual_old).and.(alpha>=(SME_alpha_min+1.d-14)) )
            x0 = x0 - alpha*dx                                      ! Reset step
            alpha = alpha*SME_alpha_mult                            ! Calculate new alpha
            if (alpha<SME_alpha_min) then
                alpha = SME_alpha_min
            endif
            x0 = x0 + alpha*dx                                      !Update step
            call residual_function(r0, x0, props, statev, param)    !Get new residual
            
            residual = sqrt(sum(r0**2))
            if (badcontent0(residual)) then
                write(*,*) 'Residual NaN or Inf'
                lconv = .false.
                return
            endif
        enddo
        residual_old = residual
        
        
        ! If diagnostic, then save the current state
        if (diagnostic) then
            rvalues_full(:, 1+k1) = r0
            xvalues_full(:, 1+k1) = x0
        endif
        
        ! If print residual, print
        if (printinfo==1) then
            write(*,"(I2,E15.5)") k1, residual
        elseif (printinfo==2) then
            write(*,"(A,I2)") 'x_', k1
            write(*,*) x0
            write(*,"(A,E15.5)") 'r(x): ', residual
            write(*,*) r0
            write(*,*) ' '
        endif
        
        if (residual<SME_itertol) then                      !Check for convergence
            lconv = .true.                                  !Solution converged
            call invert_jacobian(jac_inv, ipiv)             !Calculate inverse jacobian
            
            ! If diagnostic, allocate and save the diagnostic information to rvalues and xvalues as requested
            if (diagnostic) then
                if (present(rvalues)) then
                    allocate(rvalues(size(x0), 1+k1))
                    rvalues = rvalues_full(:, 1:(1+k1))
                endif
                if (present(xvalues)) then
                    allocate(xvalues(size(x0), 1+k1))
                    xvalues = xvalues_full(:, 1:(1+k1))
                endif
            endif
            
            return  !End subroutine when convergence achieved
        endif !End check of convergence
    enddo   !End loop
    lconv = .false. !If subroutine didn't end from convergence, r0>SME_itertol    

end subroutine newton_raphson_ana

function badcontent0(scalar) result(isbad)
implicit none
    double precision scalar
    logical          isbad
    isbad = (scalar/=scalar).or.(abs(scalar)>=infty)
    !           NaN /= NaN      infty largest representable double precision (def as parameter)
end function

function badcontent1(vector) result(isbad)
implicit none
    double precision vector(:)
    logical          isbad
    isbad = (any(vector/=vector)).or.(any(abs(vector)>=infty))
    !           NaN /= NaN                          infty largest representable double precision (def as parameter)
end function

function badcontent2(matrix) result(isbad)
implicit none
    double precision matrix(:,:)
    logical          isbad
    isbad = (any(matrix/=matrix)).or.(any(abs(matrix)>=infty))  ! See badcontent1
end function


! ================================================================================================
! Test subroutines (use for reference and testing of modifications)
subroutine test_residual(r, x, props, statev, param)
    implicit none
    double precision, intent(in)                :: x(:) 
    double precision, intent(in), optional      :: props(:), statev(:), param(:)      !Residual function parameters
    double precision, allocatable, intent(out)  :: r(:)
    allocate(r(size(x)))
    if (param(1)==1.d0) then
        r(1) = x(1)+2*x(2)
        r(2) = 3*x(1)+4*x(2)
        r(3) = x(3)
    else
        r(1) = x(1)*x(2)+1.0
        r(2) = x(1)**2-x(2)**2
        r(3) = x(3)
    endif
end subroutine test_residual

subroutine test_anajac(ajacobian, x, props, statev, param)
    implicit none
    double precision, intent(in)                :: x(:)
    double precision, intent(in), optional      :: props(:), statev(:), param(:)      !parameters
    double precision, allocatable, intent(out)  :: ajacobian(:,:)
    
    allocate(ajacobian(size(x),size(x)))
    if (param(1)==1.d0) then
        ajacobian(1,:) = (/1.d0, 2.d0, 0.d0/)
        ajacobian(2,:) = (/3.d0, 4.d0, 0.d0/)
        ajacobian(3,:) = (/0.d0, 0.d0, 1.d0/)
    else
        ajacobian(1,:) = (/x(2), x(1), 0.d0/)
        ajacobian(2,:) = (/2.d0*x(1), -2.d0*x(2), 0.d0/)
        ajacobian(3,:) = (/0.d0, 0.d0, 1.d0/)
    endif
end subroutine test_anajac

end module
