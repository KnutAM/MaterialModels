! CONTRACTIONS
! V9_d_v3           Contraction between 2nd order tensor and 1st order tensor
! V9_d_V9           Contraction between two 2nd order tensors with 9 components
! V9_dd_V9          Double contraction between two 2th order tensor of size (9)
! V9x9_dd_V9x9      Double contraction between two 4th order tensor of size (9,9)
! V9x9_dd_V9        Double contraction between a 4th order and a 2nd order tensor of size (9,9) and (9) respectively
! V9_dd_V9x9        Double contraction between a 2nd order and a 4th order tensor of size (9) and (9,9) respectively
! V9x9x9_dd_V9      Double contraction between a 6th order and a 2nd order tensor of size (9,9,9) and (9) respectively
!
! Special cases
! V9dV9dV9          Contraction between three 2nd order tensors with 9 components
! V9_d_V9x9         Single contraction between a 2nd and a 4th order tensor
! V9x9_d_V9         Single contraction between a 4th and a 2nd order tensor
!
! Operations on ordinary matrices
! M_x_ve            Multiplication between a matrix and vector
! ve_x_M            Multiplication between a vector and matrix

Module Contractions
use Conversions

Implicit none
    CONTAINS
    
!==============================================================================     
FUNCTION V9_d_v3(a,b) RESULT (c)
Implicit none

    INTEGER i,j
    DOUBLE PRECISION, INTENT (IN) :: a(9),b(3)
    DOUBLE PRECISION :: c(3)
!                              |1 4 7|
!   [1 2 3 4 5 6 7 8 9]   =>   |8 2 5|
!                              |6 9 3|
    c(1)=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
    c(2)=a(8)*b(1)+a(2)*b(2)+a(5)*b(3)
    c(3)=a(6)*b(1)+a(9)*b(2)+a(3)*b(3)
      
END FUNCTION 
!==============================================================================     

!============================================================================== 
FUNCTION V9_d_V9(a,b) RESULT (c)
Implicit none

    INTEGER i,j
    DOUBLE PRECISION, INTENT (IN) :: a(9),b(9)
    DOUBLE PRECISION :: c(9)

    c(1)=a(1)*b(1) + a(7)* b(6) + a(4)* b(8)
    c(2)=a(2)*b(2) + a(8)* b(4) + a(5)* b(9)
    c(3)=a(3)*b(3) + a(9)* b(5) + a(6)* b(7)
    c(4)=a(4)*b(2) + a(1)* b(4) + a(7)* b(9)
    c(5)=a(5)*b(3) + a(2)* b(5) + a(8)* b(7)
    c(6)=a(6)*b(1) + a(3)* b(6) + a(9)* b(8)
    c(7)=a(7)*b(3) + a(4)* b(5) + a(1)* b(7)
    c(8)=a(8)*b(1) + a(5)* b(6) + a(2)* b(8)
    c(9)=a(9)*b(2) + a(6)* b(4) + a(3)* b(9)

END FUNCTION
!============================================================================== 

!============================================================================== 
FUNCTION V9_dd_V9(a,b) RESULT (c)
Implicit none
    DOUBLE PRECISION, INTENT (IN) :: a(9),b(9)
    DOUBLE PRECISION :: c
    
    !c=tr_V9(V9_d_V9(trans_V9(a),b))
    c = sum(a*b)
    
END FUNCTION
!==============================================================================

!==============================================================================
FUNCTION V9x9_dd_V9x9(a,b) RESULT (c)
Implicit none
! Did a comparison between looping manually(1), matmul(2) and dgemm(3) for 10^6 calls
! (1) ~14.6s
! (2) ~ 8.5s 
! (3) ~ 1.4s
! Clearly, dgemm should be used as much as possible for 9x9.
! Knut Andreas Meyer, 2017-02-08
!
      DOUBLE PRECISION, INTENT (IN) :: a(9,9),b(9,9)
      DOUBLE PRECISION :: c(9,9)
      INTEGER I,J,K
!      c = 0.d0
!      do I=1,9
!        do K=1,9
!          do J=1,9
!            c(I,K)=c(I,K)+a(I,J)*b(J,K)
!          enddo
!        enddo
!      enddo
      call dgemm('N','N', 9, 9, 9, 1.d0,a, 9, b, 9, 0.d0, c, 9)
      
END FUNCTION
!============================================================================== 

!============================================================================== 
FUNCTION V9x9_dd_V9(a,b) RESULT (c)
Implicit none

      DOUBLE PRECISION,INTENT(IN)::a(9,9),b(9)
      DOUBLE PRECISION::c(9)
      INTEGER I,J
!      c=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/);
!      do I=1,9
!        do J=1,9
!            c(I)=c(I)+a(I,J)*b(J)
!        enddo
!      enddo
      call dgemv('N', 9, 9, 1.d0, a, 9, b, 1, 0.d0, c, 1)

END FUNCTION
!==============================================================================

!==============================================================================        
FUNCTION V9_dd_V9x9(a,b) RESULT (c)
Implicit none

      DOUBLE PRECISION, INTENT (IN) :: a(9),b(9,9)
      DOUBLE PRECISION              :: c(9)
      INTEGER I,K
      !c=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/);
      !do K=1,9;do I=1,9;c(K)=c(K)+a(I)*b(I,K);enddo;enddo
      
      call dgemv('T', 9, 9, 1.d0, b, 9, a, 1, 0.d0, c, 1)
      
END FUNCTION
!==============================================================================

!==============================================================================       
FUNCTION V9x9x9_dd_V9(a,b) RESULT (c)
!============================================================================== 
! Purpose: 6th order tensor double contracted with a 2nd order tensor
!
! Written by: Nasim Larijani 
!
! Function that makes c_mnkl = a_mnpqkl * b_pq, NL 2011 Nov 15
!==============================================================================  
Implicit none

      DOUBLE PRECISION, INTENT (IN) :: a(9,9,9),b(9)
      DOUBLE PRECISION              :: c(9,9)
      INTEGER M,K,P
      
      do M=1,9
        do K=1,9
              c(M,K)=0.d0
              do P=1,9
                    c(M,K)=c(M,K)+a(M,P,K)*b(P)
              enddo
        enddo
      enddo
      
           
END FUNCTION
!==============================================================================

!! Special contractions
!==============================================================================     
FUNCTION V9dV9dV9(a,b,c) RESULT (d)
! Old name: us_9d9d9
implicit none

      DOUBLE PRECISION, INTENT (IN) :: a(9),b(9),c(9)
      DOUBLE PRECISION :: d(9)
      d=V9_d_V9(a,V9_d_V9(b,c))
      
END FUNCTION
!==============================================================================

!==============================================================================
FUNCTION V9_d_V9x9(a,b) RESULT (c)
Implicit none

      DOUBLE PRECISION,INTENT(IN)   :: a(9), b(9,9)
      DOUBLE PRECISION              :: c(9,9)
      DOUBLE PRECISION              :: a_2(3,3), b_4(3,3,3,3), c_4(3,3,3,3)
      INTEGER i, j, k, l
      
      c_4 = 0.d0
      
      ! Convert to a real fouth order tensor
      b_4 = V9x9_2_T4(b)
      a_2 = V9_2_M(a)

      Do i=1,3
        Do j=1,3
          Do k=1,3
            Do l=1,3
      
                c_4(i,j,k,l) = sum(a_2(i,1:3)*b_4(1:3,j,k,l))
         
            End do
          End do
        End do
      End do
    

    ! Convert back to Voigt format
    c = T4_2_V9x9(c_4)
      
END FUNCTION
!==============================================================================

!==============================================================================    
FUNCTION V9x9_d_V9(b, a) RESULT (c)
Implicit none

      DOUBLE PRECISION,INTENT(IN)   :: a(9), b(9,9)
      DOUBLE PRECISION              :: c(9,9)
      DOUBLE PRECISION              :: a_2(3,3), b_4(3,3,3,3), c_4(3,3,3,3)
      INTEGER i, j, k, l
      
      c_4 = 0.d0
      
      ! Convert to a real fouth order tensor
      b_4 = V9x9_2_T4(b)
      a_2 = V9_2_M(a)

      Do i=1,3
        Do j=1,3
          Do k=1,3
            Do l=1,3
      
                c_4(i,j,k,l) = sum( b_4(i,j,k,1:3)*a_2(1:3,l) )
         
            End do
          End do
        End do
      End do
    

    ! Convert back to Voigt format
    c = T4_2_V9x9(c_4)
      
END FUNCTION
!==============================================================================


!! Ordinary matrix operations (should perhaps update function and issue warning saying to use built in functions)?
!============================================================================== 
FUNCTION M_x_ve(a,b) RESULT (c)
! Old name: us_matrix_m_vector
Implicit none
	
    DOUBLE PRECISION, INTENT (IN) :: a(:,:)
    DOUBLE PRECISION, INTENT (IN) :: b(:)
    DOUBLE PRECISION              :: c(size(b))
    
    INTEGER I,J
    DO I=1,SIZE(b)
       c(I)=0.d0
       DO J=1,SIZE(b)
          c(I)=c(I)+a(I,J)*b(j)
       ENDDO
    ENDDO

END FUNCTION 
!==============================================================================

!============================================================================== 
FUNCTION ve_x_M(b,a) RESULT (c)
Implicit none

    DOUBLE PRECISION, DIMENSION(:), INTENT (IN) :: b
    DOUBLE PRECISION, DIMENSION(:,:), INTENT (IN) :: a
    DOUBLE PRECISION, DIMENSION(SIZE(b)) :: c
    INTEGER I,J
    
    DO I=1,SIZE(b)
       c(I)=0.d0
       DO J=1,SIZE(b)
          c(I)=c(I)+a(j,i)*b(j)
       ENDDO
    ENDDO

END FUNCTION
!============================================================================== 

END MODULE