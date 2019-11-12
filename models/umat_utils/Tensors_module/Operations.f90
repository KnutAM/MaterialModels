! OPERATIONS        Operations on itself
! trans_V9          Transpose of a 2nd order tensor with 9 components           
! dev_V9            Deviatoric part of a 2nd order tensor with 9 components     
! tr_V9             Trace of a 9 component voigt vector                         
! det_V9            Determinant of 2nd order 9-component tensor                 
! inv_V9            Inverse of 2nd order 9-component tensor                     
!
! Operations on full matrices
! det_M             Determinant of 1x1, 2x2 or 3x3 matrix                       
! inv_M3            Subroutine to calculate inverse of 3x3 matrix               
! inv_M             Subroutine to calculate inverse of general matrix           
! trans_M           Transpose of square matrix                                  
! tr_M              Trace of a general square matrix                            
! norm_V3           Norm of a 1st order tensor with 3 components
! norm_V9           Norm of a general vector

Module Operations
use Conversions
use Contractions
Implicit none
	CONTAINS
!============================================================================== 
FUNCTION trans_V9(b) RESULT (a)
Implicit none

    DOUBLE PRECISION, INTENT (IN) :: b(9)
    DOUBLE PRECISION :: a(9)

    a(1:3)=b(1:3)
    a((/4,5,6,7,8,9/))=b((/8,9,7,6,4,5/))

END FUNCTION
!==============================================================================

!============================================================================== 
FUNCTION dev_V9(A) RESULT (Ad)
Implicit none

    DOUBLE PRECISION, INTENT(IN)  :: A(1:9)
    DOUBLE PRECISION              :: Ad(1:9),Avol


    Avol = SUM(A(1:3))
    Ad=A; Ad(1:3)=A(1:3)-Avol/3.d0

END FUNCTION
!============================================================================== 

!============================================================================== 
FUNCTION tr_V9(a) RESULT(tra)
Implicit none

      DOUBLE PRECISION, INTENT (IN) :: a(9)
      DOUBLE PRECISION :: tra
      tra=sum(a(1:3)); 
      
END FUNCTION 
!==============================================================================

!==============================================================================       
FUNCTION det_V9(A) RESULT(detA)
Implicit none

    DOUBLE PRECISION, INTENT (IN) :: A(9)
    DOUBLE PRECISION :: detA
    detA=-A(7)*A(2)*A(6)+A(4)*A(5)*A(6)+  &
          A(7)*A(8)*A(9)-A(1)*A(5)*A(9)-  &
          A(4)*A(8)*A(3)+A(1)*A(2)*A(3)  
           
END FUNCTION 
!==============================================================================

!==============================================================================
FUNCTION inv_V9(A) RESULT (inv_A) 
Implicit none

    DOUBLE PRECISION      ::    A(9), inv_A(9),det_A
  
    det_A =  -A(7)*A(2)*A(6)+A(4)*A(5)*A(6)+&
              A(7)*A(8)*A(9)-A(1)*A(5)*A(9)-&
              A(4)*A(8)*A(3)+A(1)*A(2)*A(3)

    IF ( ABS(det_A) .LT. 1.0e-25) THEN
       write(*,*) ''
       write(*,*) 'WARNING! Matrix close to being singular. (inv_V9b)'
       write(*,*) ''
       det_A=1.d-15
    ENDIF

    inv_A(1)=(-A(5)*A(9)+A(2)*A(3))/det_A
    inv_A(4)=( A(7)*A(9)-A(4)*A(3))/det_A
    inv_A(7)=( A(5)*A(4)-A(2)*A(7))/det_A

    inv_A(8)=( A(5)*A(6)-A(8)*A(3))/det_A
    inv_A(2)=(-A(7)*A(6)+A(1)*A(3))/det_A
    inv_A(5)=(-A(5)*A(1)+A(8)*A(7))/det_A

    inv_A(6)=(-A(2)*A(6)+A(8)*A(9))/det_A
    inv_A(9)=( A(4)*A(6)-A(1)*A(9))/det_A
    inv_A(3)=( A(2)*A(1)-A(8)*A(4))/det_A    

END FUNCTION
!==============================================================================

!==============================================================================       
FUNCTION invar_V9(a) RESULT(i)
Implicit none
      
      DOUBLE PRECISION, INTENT (IN) :: a(9)
      DOUBLE PRECISION :: i(3),aa(9)
      aa=V9_d_V9(a,a)
      i(1) = sum(a(1:3))
      i(2) = sum(a*a)
      i(3) = sum(aa*a); 
      
END FUNCTION 
!============================================================================== 


! OPERATIONS ON MATRICES

!==============================================================================
FUNCTION det_M(A) RESULT (DET_A)
Implicit none

    INTEGER              :: N
    DOUBLE PRECISION     :: A(:,:), DET_A
      
    DET_A=0.d0
    N = size(A,1)
    
    SELECT CASE (N)
        CASE(1)
           DET_A=A(1,1)
        CASE(2)
           DET_A=A(1,1)*A(2,2)-A(1,2)*A(2,1)
        CASE(3)
           DET_A=A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+&
           A(1,3)*A(2,1)*A(3,2)-A(1,1)*A(2,3)*A(3,2)-&
           A(1,2)*A(2,1)*A(3,3)-A(1,3)*A(2,2)*A(3,1)         
    END SELECT
      
END FUNCTION
!==============================================================================

!==============================================================================
recursive SUBROUTINE inv_M3(A,inv_A) 
Implicit none
      
    DOUBLE PRECISION A(3,3), inv_A(3,3), det_A
      
    det_A =  -A(1,3)*A(2,2)*A(3,1)+A(1,2)*A(2,3)*A(3,1)+&
              A(1,3)*A(2,1)*A(3,2)-A(1,1)*A(2,3)*A(3,2)-&
              A(1,2)*A(2,1)*A(3,3)+A(1,1)*A(2,2)*A(3,3)


    IF ( ABS(det_A) .LT. 1.0e-25) THEN
       write(*,*) ''
       write(*,*) 'WARNING! Matrix close to being singular. (inv_M3)'
       write(*,*) ''
       det_A=1.d-15
    ENDIF

    inv_A(1,1)=(-A(2,3)*A(3,2)+A(2,2)*A(3,3))/det_A
    inv_A(1,2)=( A(1,3)*A(3,2)-A(1,2)*A(3,3))/det_A
    inv_A(1,3)=( A(2,3)*A(1,2)-A(2,2)*A(1,3))/det_A

    inv_A(2,1)=( A(2,3)*A(3,1)-A(2,1)*A(3,3))/det_A
    inv_A(2,2)=(-A(1,3)*A(3,1)+A(1,1)*A(3,3))/det_A
    inv_A(2,3)=(-A(2,3)*A(1,1)+A(2,1)*A(1,3))/det_A

    inv_A(3,1)=(-A(2,2)*A(3,1)+A(2,1)*A(3,2))/det_A
    inv_A(3,2)=( A(1,2)*A(3,1)-A(1,1)*A(3,2))/det_A
    inv_A(3,3)=( A(2,2)*A(1,1)-A(2,1)*A(1,2))/det_A
    

END SUBROUTINE
!==============================================================================

!==============================================================================
recursive SUBROUTINE inv_M(A,inv_A)
!   Purpose: Inverts a matrix (taken from Fortran 90 Cooper Redwine).
Implicit none

    INTEGER ::  n
    DOUBLE PRECISION, INTENT(in)  :: A(:,:)
    DOUBLE PRECISION, INTENT(out) :: inv_A(:,:) 
    DOUBLE PRECISION, ALLOCATABLE :: temp(:)
    DOUBLE PRECISION, ALLOCATABLE :: C(:,:)
    INTEGER ::  i, j, m
    
    n = size(A,1) 
	m = size(A,2) 
     
    IF (n/=m) THEN
      write(*,*) 'ERROR! Matrix is not square. (inv_M)'
    END IF

    ALLOCATE(C(n,2*n), temp(2*n))
	C(:,1:n) = A      
	C(:,n+1:2*n) = 0.0
    DO i=1,n
      C(i,n+i) = 1.0
    END DO

    ! Gauss-Jordan elimination

    DO j=1,n !Loop over first n columns of matrix
      
      ! Find row index m of element with largest magnitude
      m = j  ! Start on the diagonal
      DO i=j+1,n
        IF (ABS(C(i,j)) > ABS(C(m,j))) m = i
      END DO

      ! Exchange row m with row j
      temp = C(j,:); C(j,:) = C(m,:); C(m,:) = temp

      C(j,:) = C(j,:)/C(j,j)  ! Divide row j by a(j,j)

      DO i=1,n  ! Loop over each row of matrix
         ! Subtract a(i,j) times row j from row i
         IF (i/=j) C(i,:) = C(i,:)-C(i,j)*C(j,:)
      END DO
      
    END DO
    
    inv_A = C(:,n+1:2*n)
	
END SUBROUTINE 
!==============================================================================

!==============================================================================     
FUNCTION trans_M(a) RESULT(at)
Implicit none

    DOUBLE PRECISION, INTENT (IN) :: a(:,:)
    DOUBLE PRECISION              :: at(size(a,2),size(a,1))
    INTEGER                       :: I,J

    DO I=1,size(a,2)
       DO J=1,size(a,1)
          at(I,J)=a(J,I)
       ENDDO
    ENDDO
     
END FUNCTION
!============================================================================== 

!==============================================================================       
FUNCTION tr_M(a) RESULT (tra)
Implicit none

      INTEGER i
      DOUBLE PRECISION, DIMENSION(:,:), INTENT (IN) :: a
      DOUBLE PRECISION :: tra      
      if(size(a,2).ne.size(a,1))then
        write(*,*) 'not square matrix in us_trM_f';read(*,*)
      else
        tra=0.d0;DO i=1,size(a,1);tra=tra+a(i,i);ENDDO
      endif 
            
END FUNCTION 
!==============================================================================

!============================================================================== 
FUNCTION norm_V3(A) RESULT (nA)
Implicit none

    DOUBLE PRECISION, INTENT(IN) ::A(3)
    DOUBLE PRECISION :: nA
      
    nA=sqrt((A(1)**2+A(2)**2+A(3)**2))

END FUNCTION
!==============================================================================

FUNCTION norm_V9(A) RESULT (nA)
Implicit none

    DOUBLE PRECISION, INTENT(IN) ::A(1:9)
    DOUBLE PRECISION :: nA
      
    nA=dsqrt( dabs(sum(A*A)) )

END FUNCTION
!============================================================================== 

END MODULE