! ABAQUS CONVERSIONS
!   Functions               Input           Description
! sigma_v9_2_sigma_abaqus   sigma_v9        Convert sigma in 9-voigt to abaqus sigma 6-voigt and abaqus index order
! sigma_abaqus_2_sigma_v9   sigma_abaqus    Convert sigma in abaqus sigma 6-voigt and abaqus index order to 9-voigt
! eps_abaqus_2_eps_v9       eps_aba         Convert epsilon in abaqus output format to tensor module format
! eps_v9_2_eps_abaqus       eps_v9          Convert epsilon in tensor module format to abaqus output format
!
! dtaudF_2_DDSDDE           dtaudF, F       Convert dtaudF in 9x9voigt to Abaqus DDSDDE in 6x6 voigt and abaqus index order

module ABAQUS_Conversions
use Conversions
use Open_products
use Contractions
use Standard_tensors
use Stiffness_Conversions
Implicit none
    contains

!==============================================================================     
function sigma_v9_2_sigma_abaqus(sigma_v9) result(sigma_abaqus)
    !Convert sigma in tensor module format to abaqus output format
   implicit none
   double precision, intent(in) :: sigma_v9(9)
   double precision             :: sigma_abaqus(6)
   double precision             :: tmp(9)
    tmp = 0.5*(sigma_v9 + trans_V9(sigma_v9))   !Symmetrize
    sigma_abaqus = tmp((/1, 2, 3, 4, 6, 5/))      !Convert indicies

end function sigma_v9_2_sigma_abaqus
!==============================================================================    

!==============================================================================    
function sigma_abaqus_2_sigma_v9(sigma_abaqus) result(sigma_v9)
    !Convert sigma in abaqus output format to tensor module format
    implicit none
    double precision, intent(in) :: sigma_abaqus(6)
    double precision             :: sigma_v9(9)
    sigma_v9  = sigma_abaqus((/1, 2, 3, 4, 6, 5, 5, 4, 6/))      !Convert indicies

end function sigma_abaqus_2_sigma_v9
!==============================================================================    

!==============================================================================    
function eps_abaqus_2_eps_v9(eps_aba) result(eps_v9)
    !Convert epsilon in abaqus output format to tensor module format
    implicit none
    double precision, intent(in) :: eps_aba(6)
    double precision             :: eps_v9(9)
    eps_v9  = eps_aba((/1, 2, 3, 4, 6, 5, 5, 4, 6/))      !Convert indicies
    eps_v9(4:9) = eps_v9(4:9)/2

end function eps_abaqus_2_eps_v9
!==============================================================================    

!==============================================================================    
function eps_v9_2_eps_abaqus(eps_v9) result(eps_aba)
    !Convert epsilon in tensor module format to abaqus output format
    implicit none
    double precision, intent(in) :: eps_v9(9)
    double precision             :: eps_aba(6)
    eps_aba(1:3) = eps_v9((/1, 2, 3/))
    eps_aba(4:6) = eps_v9((/4, 7, 5/)) + eps_v9((/8, 6, 9/))

end function eps_v9_2_eps_abaqus
!==============================================================================    

!==============================================================================     
function dtaudF_2_ddsdde(dtaudF, F) result(ddsdde)
    implicit none
    double precision, intent(in) :: dtaudF(9,9), F(9)
    double precision             :: tmp(9,9), ddsdde(6,6)
    integer                      :: posvec(6)
    
    tmp = (1.d0/det_V9(F)) * V9x9_dd_V9x9(dtaudF, op_a_V9(eye_V9(), trans_V9(F)))
    tmp = minsym_V9x9(tmp)
    
    posvec = (/1, 2, 3, 4, 6, 5/)   ! Index conversion
    ddsdde = tmp(posvec, posvec)
    
end function dtaudF_2_ddsdde
!==============================================================================  

end module