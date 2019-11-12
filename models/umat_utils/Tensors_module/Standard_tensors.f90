! STANDARD_TENSORS (PROJECTION TENSORS)
! eye_V9            2nd order identity tensor (9 components)
! idev_v9x9         4th order deviatoric identity tensor
module standard_tensors
use open_products
implicit none
    contains

!==============================================================================     
function eye_v9() result(eye)
   implicit none
   double precision             :: eye(9)

   eye = 0.0d0
   
   eye(1:3) = (/ 1.0d0, 1.0d0, 1.0d0 /)

end function eye_V9
!==============================================================================     

!==============================================================================     
function idev_v9x9() result(idev)
    implicit none
    double precision             :: idev(9,9), eye(9)
    
    eye = eye_v9()
    idev = op_a_v9(eye,eye) - op_v9(eye, eye)/3.d0
   
end function idev_v9x9
!==============================================================================     

end module