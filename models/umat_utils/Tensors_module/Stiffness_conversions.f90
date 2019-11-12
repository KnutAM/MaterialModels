! STIFFNESS_CONVERSIONS
! dSdE_2_dPdF       Convert dSdE to dPdF
! dSdE_2_dtaudF     Convert dSdE to d(tau)dF
! dPdF_2_dtaudF     Convert dPdF to d(tau)dF

Module Stiffness_Conversions
use Conversions
use Open_products
use Contractions
use Standard_tensors
use Operations
Implicit none
    CONTAINS

!==============================================================================     
FUNCTION dSdE_2_dPdF(dSdE,F,S) result(dPdF)
   implicit none
   double precision, intent(in) :: dSdE(9,9), F(9), S(9)
   double precision             :: dPdF(9,9), I2(9), Ft(9)
    I2 = eye_V9()
    Ft = trans_V9(F)
    dPdF = op_a_V9(I2, S) + (1/2)*V9x9_dd_V9x9(V9x9_dd_V9x9(op_a_V9(F, I2), dSdE), op_b_V9(I2, Ft)+op_a_V9(Ft,I2))

END FUNCTION dSdE_2_dPdF
!==============================================================================    

!==============================================================================     
FUNCTION dSdE_2_dtaudF(dSdE,F,S) result(dtaudF)
    !S = 2nd Piola Kirchoff
    !E = 0.5(C-I2), C=F^T*F
    !F = Deformation gradient
   implicit none
   double precision, intent(in) :: dSdE(9,9), F(9), S(9)
   double precision             :: dtaudF(9,9), I2(9), Ft(9)
    I2 = eye_V9()
    Ft = trans_V9(F)
    dtaudF = op_a_V9(I2, V9_d_V9(F, S)) + op_b_V9(V9_d_V9(F, S), I2) &
        + 0.5*V9x9_dd_V9x9(V9x9_dd_V9x9(op_a_V9(F, F), dSdE), op_b_V9(I2, Ft)+op_a_V9(Ft, I2))

END FUNCTION dSdE_2_dtaudF
!==============================================================================  

!==============================================================================     
FUNCTION dPdF_2_dtaudF(dPdF,F,P) result(dtaudF)
   implicit none
   double precision, intent(in) :: dPdF(9,9), F(9), P(9)
   double precision             :: dtaudF(9,9), I2(9), Ft(9)
    I2 = eye_V9()
    Ft = trans_V9(F)
    dtaudF = V9x9_dd_V9x9(op_a_V9(I2, F), dPdF) + op_b_V9(P, I2)

END FUNCTION dPdF_2_dtaudF
!============================================================================== 

end module