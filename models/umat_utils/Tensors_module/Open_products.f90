! OPEN_PRODUCTS
! op_ve             Open product between two 1st order tensors (3 components) [i x j] giving matrix 3x3 out
! op_ve_9           Open product between two 1st order tensors (3 components) [i x j] giving 9 voigt out
! op_V9             Open product between two 2nd order tensors (9 components) [ij x kl]
! op_a_V9           Non-standard open product "above" between two 2nd order tensors (9 components) [ik x jl]
! op_b_V9           Non-standard open product "below" between two 2nd order tensors (9 components) [il x jk]
! Pop_V9            c_ijkl = 0.25*( a_ik b_jl+a_il b_jk+b_ik a_jl+b_il a_jk )
!
! Open products between 4th and - 2nd order tensors
! v9_op_v9x9        a_ij B_klmn - 2nd order in first position
! v9x9_op_v9        A_ijkl b_mn - 2nd order in last position
! v9_opc_v9x9       c_kl D_ijmn - 2nd order in center position
! ilmn_jk           c_ijklmn = a_ilmn * b_jk
! ikmn_jl           c_ijklmn = a_ikmn * b_jl
!
! 6-component functions
! op_V6             Open product between two 2nd order tensors (6 components)


Module Open_products
Implicit none
CONTAINS

!==============================================================================  
FUNCTION op_ve(a,b) RESULT (c)
implicit none

    DOUBLE PRECISION, INTENT (IN) :: a(3),b(3)
    DOUBLE PRECISION :: c(3,3)
    INTEGER I,J
    DO I=1,3
       DO J=1,3
             c(I,J)=a(I)*b(J)
       ENDDO
    ENDDO  
      
END FUNCTION 
!==============================================================================   

!==============================================================================   
function op_ve_9(a,b) result (c)
    implicit none
    double precision, intent(in) :: a(3), b(3)
    double precision :: c(9)
    c(1) = a(1)*b(1)
    c(2) = a(2)*b(2)
    c(3) = a(3)*b(3)
    c(4) = a(1)*b(2)
    c(5) = a(2)*b(3)
    c(6) = a(3)*b(1)
    c(7) = a(1)*b(3)
    c(8) = a(2)*b(1)
    c(9) = a(3)*b(2)
end function op_ve_9
!==============================================================================  

!==============================================================================  
FUNCTION op_V9(a,b) RESULT (c)
implicit none

    DOUBLE PRECISION, INTENT (IN) :: a(9),b(9)
    DOUBLE PRECISION :: c(9,9)
    INTEGER I,J
    DO I=1,9
       DO J=1,9
             c(I,J)=a(I)*b(J)
       ENDDO
    ENDDO

END FUNCTION 
!==============================================================================    

!==============================================================================   
FUNCTION op_a_V9(b,c) RESULT (a) 
implicit none

    DOUBLE PRECISION, INTENT (IN) :: b(9),c(9)
    DOUBLE PRECISION :: a(9,9)
    
    ! Derived from matlab program derive_expr
    a(1,1)=b(1)*c(1); a(1,2)=b(4)*c(4); a(1,3)=b(7)*c(7);
    a(1,4)=b(1)*c(4); a(1,5)=b(4)*c(7); a(1,6)=b(7)*c(1);
    a(1,7)=b(1)*c(7); a(1,8)=b(4)*c(1); a(1,9)=b(7)*c(4);
    a(2,1)=b(8)*c(8); a(2,2)=b(2)*c(2); a(2,3)=b(5)*c(5);
    a(2,4)=b(8)*c(2); a(2,5)=b(2)*c(5); a(2,6)=b(5)*c(8);
    a(2,7)=b(8)*c(5); a(2,8)=b(2)*c(8); a(2,9)=b(5)*c(2);
    a(3,1)=b(6)*c(6); a(3,2)=b(9)*c(9); a(3,3)=b(3)*c(3);
    a(3,4)=b(6)*c(9); a(3,5)=b(9)*c(3); a(3,6)=b(3)*c(6);
    a(3,7)=b(6)*c(3); a(3,8)=b(9)*c(6); a(3,9)=b(3)*c(9);
    a(4,1)=b(1)*c(8); a(4,2)=b(4)*c(2); a(4,3)=b(7)*c(5);
    a(4,4)=b(1)*c(2); a(4,5)=b(4)*c(5); a(4,6)=b(7)*c(8);
    a(4,7)=b(1)*c(5); a(4,8)=b(4)*c(8); a(4,9)=b(7)*c(2);
    a(5,1)=b(8)*c(6); a(5,2)=b(2)*c(9); a(5,3)=b(5)*c(3);
    a(5,4)=b(8)*c(9); a(5,5)=b(2)*c(3); a(5,6)=b(5)*c(6);
    a(5,7)=b(8)*c(3); a(5,8)=b(2)*c(6); a(5,9)=b(5)*c(9);
    a(6,1)=b(6)*c(1); a(6,2)=b(9)*c(4); a(6,3)=b(3)*c(7);
    a(6,4)=b(6)*c(4); a(6,5)=b(9)*c(7); a(6,6)=b(3)*c(1);
    a(6,7)=b(6)*c(7); a(6,8)=b(9)*c(1); a(6,9)=b(3)*c(4);
    a(7,1)=b(1)*c(6); a(7,2)=b(4)*c(9); a(7,3)=b(7)*c(3);
    a(7,4)=b(1)*c(9); a(7,5)=b(4)*c(3); a(7,6)=b(7)*c(6);
    a(7,7)=b(1)*c(3); a(7,8)=b(4)*c(6); a(7,9)=b(7)*c(9);
    a(8,1)=b(8)*c(1); a(8,2)=b(2)*c(4); a(8,3)=b(5)*c(7);
    a(8,4)=b(8)*c(4); a(8,5)=b(2)*c(7); a(8,6)=b(5)*c(1);
    a(8,7)=b(8)*c(7); a(8,8)=b(2)*c(1); a(8,9)=b(5)*c(4);
    a(9,1)=b(6)*c(8); a(9,2)=b(9)*c(2); a(9,3)=b(3)*c(5);
    a(9,4)=b(6)*c(2); a(9,5)=b(9)*c(5); a(9,6)=b(3)*c(8);
    a(9,7)=b(6)*c(5); a(9,8)=b(9)*c(8); a(9,9)=b(3)*c(2);

END FUNCTION
!============================================================================== 

!============================================================================== 
FUNCTION op_b_V9(b,c) RESULT (a)
implicit none

    DOUBLE PRECISION, INTENT (IN) :: b(9),c(9)
    DOUBLE PRECISION :: a(9,9)
    
    ! Derived from matlab program derive_expr
    a(1,1)=b(1)*c(1);a(1,2)=b(4)*c(4);a(1,3)=b(7)*c(7);
    a(1,4)=b(4)*c(1);a(1,5)=b(7)*c(4);a(1,6)=b(1)*c(7);
    a(1,7)=b(7)*c(1);a(1,8)=b(1)*c(4);a(1,9)=b(4)*c(7);
    a(2,1)=b(8)*c(8);a(2,2)=b(2)*c(2);a(2,3)=b(5)*c(5);
    a(2,4)=b(2)*c(8);a(2,5)=b(5)*c(2);a(2,6)=b(8)*c(5);
    a(2,7)=b(5)*c(8);a(2,8)=b(8)*c(2);a(2,9)=b(2)*c(5);
    a(3,1)=b(6)*c(6);a(3,2)=b(9)*c(9);a(3,3)=b(3)*c(3);
    a(3,4)=b(9)*c(6);a(3,5)=b(3)*c(9);a(3,6)=b(6)*c(3);
    a(3,7)=b(3)*c(6);a(3,8)=b(6)*c(9);a(3,9)=b(9)*c(3);
    a(4,1)=b(1)*c(8);a(4,2)=b(4)*c(2);a(4,3)=b(7)*c(5);
    a(4,4)=b(4)*c(8);a(4,5)=b(7)*c(2);a(4,6)=b(1)*c(5);
    a(4,7)=b(7)*c(8);a(4,8)=b(1)*c(2);a(4,9)=b(4)*c(5);
    a(5,1)=b(8)*c(6);a(5,2)=b(2)*c(9);a(5,3)=b(5)*c(3);
    a(5,4)=b(2)*c(6);a(5,5)=b(5)*c(9);a(5,6)=b(8)*c(3);
    a(5,7)=b(5)*c(6);a(5,8)=b(8)*c(9);a(5,9)=b(2)*c(3);
    a(6,1)=b(6)*c(1);a(6,2)=b(9)*c(4);a(6,3)=b(3)*c(7);
    a(6,4)=b(9)*c(1);a(6,5)=b(3)*c(4);a(6,6)=b(6)*c(7);
    a(6,7)=b(3)*c(1);a(6,8)=b(6)*c(4);a(6,9)=b(9)*c(7);
    a(7,1)=b(1)*c(6);a(7,2)=b(4)*c(9);a(7,3)=b(7)*c(3);
    a(7,4)=b(4)*c(6);a(7,5)=b(7)*c(9);a(7,6)=b(1)*c(3);
    a(7,7)=b(7)*c(6);a(7,8)=b(1)*c(9);a(7,9)=b(4)*c(3);
    a(8,1)=b(8)*c(1);a(8,2)=b(2)*c(4);a(8,3)=b(5)*c(7);
    a(8,4)=b(2)*c(1);a(8,5)=b(5)*c(4);a(8,6)=b(8)*c(7);
    a(8,7)=b(5)*c(1);a(8,8)=b(8)*c(4);a(8,9)=b(2)*c(7);
    a(9,1)=b(6)*c(8);a(9,2)=b(9)*c(2);a(9,3)=b(3)*c(5);
    a(9,4)=b(9)*c(8);a(9,5)=b(3)*c(2);a(9,6)=b(6)*c(5);
    a(9,7)=b(3)*c(8);a(9,8)=b(6)*c(2);a(9,9)=b(9)*c(5);

END FUNCTION
!==============================================================================

!============================================================================== 
FUNCTION Pop_V9(b,c) RESULT (a)
! P operator: c_ijkl = 0.25*( a_ik b_jl+a_il b_jk+b_ik a_jl+b_il a_jk )
implicit none

    DOUBLE PRECISION, INTENT (IN) :: b(9),c(9)
    DOUBLE PRECISION :: a(9,9)

    a=0.25d0*( op_a_V9(b,c) + op_a_V9(c,b) + op_b_V9(b,c) + op_b_V9(c,b) )
    
END FUNCTION 
!==============================================================================

! Open products between 4th and - 2nd order tensors
!==============================================================================  
function v9_op_v9x9(a,b) result (c)
    implicit none
    double precision, intent(in) :: a(9),b(9,9)
    double precision             :: c(9,9,9)
    integer k1, k2, k3
    
    do k1=1,9
        do k2=1,9
            do k3=1,9
                c(k1,k2,k3)=a(k1)*b(k2,k3)
            enddo
        enddo
    enddo
    
end function v9_op_v9x9
!==============================================================================

!==============================================================================  
function v9x9_op_v9(a,b) result (c)
    implicit none
    double precision, intent(in) :: a(9,9),b(9)
    double precision             :: c(9,9,9)
    integer k1, k2, k3
    
    do k1=1,9
        do k2=1,9
            do k3=1,9
                c(k1,k2,k3)=a(k1,k2)*b(k3)
            enddo
        enddo
    enddo
    
end function v9x9_op_v9
!==============================================================================

!==============================================================================  
function v9_opc_v9x9(a,b) result (c)
    implicit none
    double precision, intent(in) :: a(9),b(9,9)
    double precision             :: c(9,9,9)
    integer k1, k2, k3
    
    do k1=1,9
        do k2=1,9
            do k3=1,9
                c(k1,k2,k3)=a(k2)*b(k1,k3)
            enddo
        enddo
    enddo
    
end function v9_opc_v9x9
!==============================================================================

!==============================================================================
function ilmn_jk(a,b) result (c)
    !Calculates c_ijklmn = a_ilmn * b_jk
    implicit none
    double precision, intent(in) :: a(9,9),b(9)
    double precision             :: c(9,9,9)
    integer k1, k2, k3, v_voigt(9,2), m_voigt(3,3), ind1, ind2
    
    ! Define voigt numbering
    v_voigt(:,1) = (/1, 2, 3, 1, 2, 3, 1, 2, 3/)
    v_voigt(:,2) = (/1, 2, 3, 2, 3, 1, 3, 1, 2/)
    m_voigt(1,:) = (/1, 4, 7/)
    m_voigt(2,:) = (/8, 2, 5/)
    m_voigt(3,:) = (/6, 9, 3/)
    
    do k1 = 1,9 !ij
        do k2 = 1,9 !kl
            do k3 = 1,9 !mn
                ind1 = m_voigt(v_voigt(k1,1), v_voigt(k2,2)) !il (i is 1st index in k1, l is 2nd index in k2)
                ind2 = m_voigt(v_voigt(k1,2), v_voigt(k2,1)) !jk (j is 2nd index in k1, k is 1st index in k2)
                c(k1,k2,k3) = a(ind1,k3)*b(ind2)
            enddo
        enddo
    enddo
    
end function ilmn_jk
!==============================================================================

!==============================================================================
function ikmn_jl(a,b) result (c)
    !Calculates c_ijklmn = a_ikmn * b_jl
    implicit none
    double precision, intent(in) :: a(9,9),b(9)
    double precision             :: c(9,9,9)
    integer k1, k2, k3, v_voigt(9,2), m_voigt(3,3), ind1, ind2
    
    ! Define voigt numbering
    v_voigt(:,1) = (/1, 2, 3, 1, 2, 3, 1, 2, 3/)
    v_voigt(:,2) = (/1, 2, 3, 2, 3, 1, 3, 1, 2/)
    m_voigt(1,:) = (/1, 4, 7/)
    m_voigt(2,:) = (/8, 2, 5/)
    m_voigt(3,:) = (/6, 9, 3/)
    
    do k1 = 1,9 !ij
        do k2 = 1,9 !kl
            do k3 = 1,9 !mn
                ind1 = m_voigt(v_voigt(k1,1), v_voigt(k2,1)) !ik (i is 1st index in k1, k is 1st index in k2)
                ind2 = m_voigt(v_voigt(k1,2), v_voigt(k2,2)) !jl (j is 2nd index in k1, l is 2nd index in k2)
                c(k1,k2,k3) = a(ind1,k3)*b(ind2)
            enddo
        enddo
    enddo
    
end function ikmn_jl
!==============================================================================

! 6-component voigt functions
!==============================================================================  
FUNCTION op_V6(a,b) RESULT (c)
! Old name: us_open_stand6
implicit none

    DOUBLE PRECISION, INTENT (IN) :: a(6),b(6)
    DOUBLE PRECISION :: c(6,6)
    INTEGER I,J
    DO I=1,6
       DO J=1,6
          IF ((I.GE.4).AND.(J.GE.4)) THEN
             c(I,J)=a(I)*b(J)*2.d0
          ELSE
             c(I,J)=a(I)*b(J)
          ENDIF
       ENDDO
    ENDDO

END FUNCTION 
!==============================================================================  

END MODULE