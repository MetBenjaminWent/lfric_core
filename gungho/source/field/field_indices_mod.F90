!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you should have
! received as part of this distribution.
!-----------------------------------------------------------------------------
!
!> @brief Define indices for the prognostic field vectors
!>
!> @details Define and set indices for the prognostic field vectors
module field_indices_mod

  implicit none

  ! For gungho model
  integer, parameter :: igh_u = 1  ! wind
  integer, parameter :: igh_t = 2  ! potential temperature
  integer, parameter :: igh_d = 3  ! density
  integer, parameter :: igh_p = 4  ! exner pressure
  integer, parameter :: igh_uv = igh_u ! uv wind if split - must be the same as igh_u
  integer, parameter :: igh_w= 5   ! w wind if split

  ! For gravity wave mini app
  integer, parameter :: igw_u = 1  ! wind
  integer, parameter :: igw_p = 2  ! pressure
  integer, parameter :: igw_b = 3  ! buoyancy

end module field_indices_mod
