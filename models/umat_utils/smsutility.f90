!==================================================================================================================================
! MODUL: SMS,  Jože Korelc, modifikacija M Kegl,  Verzija 2004.07.18
! SMS utility routines
!==================================================================================================================================

module SMSUtility

!_Uporabljeni moduli

!_Default typing
  implicit none
!_Default exporting
  private
!_Default memory
  save

!_Procedure
  public :: SMSDeltaPart
  public :: SMSDot
  public :: SMSSum
  public :: SMSKDelta
  public :: SMSMove
  !public :: SMSZero

!==================================================================================================================================
contains
!==================================================================================================================================


pure function SMSDeltaPart(a,i,j,k)
  doubleprecision :: SMSDeltaPart
  doubleprecision, intent(in) :: a(:)
  integer, intent(in) :: i,j,k
  
  integer :: l

  l=i/j
  IF(mod(i,j).ne.0 .or. l.gt.k) then
      SMSDeltaPart = 0d0
  else
      SMSDeltaPart = a(l)
  endif

end function


!==================================================================================================================================
pure function SMSDot(a,b,n)
  doubleprecision :: SMSDot
  doubleprecision, intent(in) :: a(:),b(:)
  integer, intent(in) :: n

  SMSDot = sum(a(1:n)*b(1:n))

end function


!==================================================================================================================================
pure function SMSSum(a,n)
  doubleprecision :: SMSSum
  doubleprecision, intent(in) :: a(:)
  integer, intent(in) :: n

  SMSSum = sum(a(1:n))

end function


!==================================================================================================================================
pure function SMSKDelta(i,j)
  doubleprecision :: SMSKDelta
  integer, intent(in) :: i,j

  if( i .eq. j )then
   SMSKDelta = 1.0d0
  else
   SMSKDelta = 0.0d0
  endif

end function


!==================================================================================================================================
pure subroutine SMSMove(a,b,n)
  doubleprecision, intent(in) :: a(:)
  doubleprecision, intent(out) :: b(:)
  integer, intent(in) :: n

  b(1:n) = a(1:n)

end subroutine


!====================================================================3
end module