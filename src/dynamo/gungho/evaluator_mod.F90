!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at
! https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
! Abstract base evaluator type.
!-------------------------------------------------------------------------------
!> @brief Abstract base type for for evaluator
module evaluator_mod
use constants_mod, only: i_def
implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, abstract :: evaluator_type
  private

end type

end module evaluator_mod
