! CONVERSIONS
! V6_2_V9           Transformation of a 2nd order tensor from 6 to 9 components in Voigt format
! V9_2_V6           Transformation of a 2nd order tensor from 9 to 6 components in Voigt format
! V9_2_M            Transformation of a 2nd order tensor with 9 components to matrix form of size(3,3)
! V6_2_M            Transformation of a 2nd order tensor with 6 components to matrix form of size(3,3)
! M_2_V9            Transformation of a matrix of size(3,3) to a second order tensor form with 9 components
! M_2_V6            Transformation of a matrix of size(3,3) to a second order tensor form with 9 components
! V9x9_2_V6x6       Transformation of a 4th order tensor of size(9,9) to size(6,6)
! V6x6_2_Sym_V9x9   Transformation of a 4th order tensor of size(6,6) to size(9,9)
! V9x9_2_T4         Transformation of a 4th order tensor of size(9,9) to pure tensor form of size(3,3,3,3)
! T4_2_V9x9         Transformation of a pure tensor form of size(3,3,3,3) to 4th order tensor of size(9,9) 
! minsym_V9x9       Minor symmetric part of 9x9 4th order tensor, not implemented. Required for V9x9_2_V6x6, thus in this file
! majsym_V9x9       Major symmetric part of 9x9 4th order tensor, not implemented. Put in this file to correspond to minsym_V9x9

Module Conversions
Implicit none
    CONTAINS

!============================================================================== 
FUNCTION V6_2_V9(a) RESULT (a9)
Implicit none

    DOUBLE PRECISION, INTENT (IN) :: a(6)
    DOUBLE PRECISION :: a9(9)
    
    a9(1:6)=a
    a9(7)=a(6)
    a9(8)=a(4)
    a9(9)=a(5)

END FUNCTION
!============================================================================== 

!============================================================================== 
FUNCTION V9_2_V6(a9) RESULT (a6)
Implicit none

    DOUBLE PRECISION, INTENT (IN) :: a9(9)
    DOUBLE PRECISION :: a6(6)

    a6(1:3)=a9(1:3)
    a6(4)=(a9(4)+a9(8))/2.d0
    a6(5)=(a9(5)+a9(9))/2.d0
    a6(6)=(a9(6)+a9(7))/2.d0

END FUNCTION
!==============================================================================

!============================================================================== 
FUNCTION V9_2_M(a) RESULT (b)
Implicit none

    DOUBLE PRECISION, INTENT (IN) :: a(9)
    DOUBLE PRECISION :: b(3,3)
  
    b(1,1)=a(1)
    b(2,2)=a(2)
    b(3,3)=a(3)
    b(1,2)=a(4)
    b(2,1)=a(8)
    b(2,3)=a(5)
    b(3,2)=a(9)
    b(1,3)=a(7)
    b(3,1)=a(6)
    
END FUNCTION
!============================================================================== 

!============================================================================== 
FUNCTION V6_2_M(a) RESULT (b)
Implicit none

    DOUBLE PRECISION, INTENT (IN) :: a(6)
    DOUBLE PRECISION :: b(3,3)
  
    b(1,1)=a(1)
    b(2,2)=a(2)
    b(3,3)=a(3)
    b(1,2)=a(4)
    b(2,1)=a(4)
    b(2,3)=a(5)
    b(3,2)=a(5)
    b(1,3)=a(6)
    b(3,1)=a(6)

END FUNCTION
!============================================================================== 

!============================================================================== 
FUNCTION M_2_V9(b) RESULT (a)
! old name: US_M_2_V9
Implicit none
    
    DOUBLE PRECISION, INTENT (IN) :: b(3,3)
    DOUBLE PRECISION :: a(9)
      
    a(1)=b(1,1)
    a(2)=b(2,2)
    a(3)=b(3,3)
    a(4)=b(1,2)
    a(9)=b(3,2)
    a(5)=b(2,3)
    a(8)=b(2,1)
    a(7)=b(1,3)
    a(6)=b(3,1)

END FUNCTION 
!============================================================================== 

!============================================================================== 
FUNCTION M_2_V6(b) RESULT (a)
Implicit none

    DOUBLE PRECISION, INTENT (IN) :: b(3,3)
    DOUBLE PRECISION :: a(6)
      
    a(1)=b(1,1)
    a(2)=b(2,2)
    a(3)=b(3,3)
    a(4)=(b(1,2)+b(2,1))/2.d0
    a(5)=(b(3,2)+b(2,3))/2.d0
    a(6)=(b(1,3)+b(3,1))/2.d0

END FUNCTION 
!============================================================================== 

!============================================================================== 
FUNCTION V9x9_2_V6x6(b) RESULT (a)
! Old name: us_9x9_2_6x6
Implicit none

    DOUBLE PRECISION, INTENT (IN) :: b(9,9)
    DOUBLE PRECISION :: a(6,6),bsym(9,9)

    bsym=minsym_V9x9(b)
    a=bsym(1:6,1:6)

END FUNCTION
!============================================================================== 

!============================================================================== 
FUNCTION V6x6_2_Sym_V9x9(b) RESULT (a)
! Old name: us_6x6_2_9x9
Implicit none

    DOUBLE PRECISION, INTENT (IN) :: b(6,6)
    DOUBLE PRECISION :: a(9,9)

    a(1:6,1:6) = b
    a(7,1:6)   = b(6,:)
    a(8,1:6)   = b(4,:)
    a(9,1:6)   = b(5,:)
    a(:,7)     = a(:,6)
    a(:,8)     = a(:,4)
    a(:,9)     = a(:,5)
    a(4:9,4:9) = 0.5d0*a(4:9,4:9)
    

END FUNCTION
!============================================================================== 

!============================================================================== 
FUNCTION minsym_V9x9(b) RESULT (a)
Implicit none

    DOUBLE PRECISION, INTENT (IN) :: b(9,9)
    DOUBLE PRECISION :: a(9,9)
    integer:: upper(3), lower(3)
    upper = (/4, 5, 7/) !Indicies on upper diagonal
    lower = (/8, 9, 6/) !Indicies on lower diagonal
    a = b
    
    a(:, upper) = (a(:, upper) + a(:, lower))/2.d0  !Symmetrize columns
    a(:, lower) = a(:, upper)
    
    a(upper, :) = (a(upper, :) + a(lower, :))/2.d0  !Symmetrize rows
    a(lower, :) = a(upper,:)

END FUNCTION
!============================================================================== 

!==============================================================================
FUNCTION majsym_V9x9(b) RESULT (a)
! majsym_V9x9       Major symmetric part of a 4th order tensor of size (9,9)
Implicit none

	DOUBLE PRECISION, INTENT (IN) :: b(9,9)
	DOUBLE PRECISION :: a(9,9)
	
	a = (b+transpose(b))/2.d0

END FUNCTION
!============================================================================== 

!!============================================================================== 
FUNCTION T4_2_V9x9(M) RESULT (V9x9)
! m_2_v9x9          Convert 4th order 3x3x3x3 to 9x9 voigt matrix
    implicit none
    double precision, intent(in) :: M(3,3,3,3)
    double precision :: V9x9(9,9)
    integer:: Mindices(3,3), i, j, k, l
    Mindices(1,:) = (/1, 4, 7/)
    Mindices(2,:) = (/8, 2, 5/)
    Mindices(3,:) = (/6, 9, 3/)
    do i = 1,3
        do j = 1,3
            do k = 1,3
                do l = 1,3
                    V9x9(Mindices(i,j),Mindices(k,l)) = M(i,j,k,l)
                end do
            end do
        end do 
    end do
    
END FUNCTION T4_2_V9x9
!==============================================================================

!============================================================================== 
FUNCTION V9x9_2_T4(V9x9) RESULT (M)
! v9x9_2_m          Convert 4th order 9x9 voigt matrix to 3x3x3x3
    implicit none
    double precision, intent(in) :: V9x9(9,9)
    double precision :: M(3,3,3,3)
    integer:: Mindices(3,3), i, j, k, l
    Mindices(1,:) = (/1, 4, 7/)
    Mindices(2,:) = (/8, 2, 5/)
    Mindices(3,:) = (/6, 9, 3/)
    
    do i = 1,3
        do j = 1,3
            do k = 1,3
                do l = 1,3
                    M(i,j,k,l) = V9x9(Mindices(i,j),Mindices(k,l))
                end do
            end do
        end do 
    end do
    
END FUNCTION V9x9_2_T4
!============================================================================== 


END MODULE