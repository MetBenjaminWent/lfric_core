!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the finite-volume divergence in either the x or
!>        y direction.
!> @details The Cosmic scheme updates density in the x, y and z directions
!>          separately.
!>          This code calculates the divergence for either the x or y direction.
!>          The z direction update has not been developed yet.
!>          The scheme is a simple finite difference of the fluxes at opposite
!>          cell edges and is designed to work only with lowest order W2 and W3
!>          spaces.

module fv_divergence_kernel_mod

use argument_mod,  only : arg_type, func_type,                  &
                          GH_FIELD, GH_WRITE, GH_READ,          &
                          W0, W2, W3, GH_BASIS, CELLS
use constants_mod, only : r_def
use kernel_mod,    only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: fv_divergence_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE, W3),                             &
       arg_type(GH_FIELD,   GH_READ,  W2),                             &
       arg_type(GH_FIELD,   GH_READ,  W3)                              &
       /)
  type(func_type) :: meta_funcs(3) = (/                                &
       func_type(W3, GH_BASIS),                                        &
       func_type(W2, GH_BASIS),                                        &
       func_type(W3, GH_BASIS)                                         &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass ::fv_divergence_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface fv_divergence_kernel_type
   module procedure fv_divergence_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public fv_divergence_code
contains

type(fv_divergence_kernel_type) function fv_divergence_kernel_constructor() result(self)
  return
end function fv_divergence_kernel_constructor

!> @brief Computes the finite-volume divergence in either the x or y direction.
!! @param[in]  nlayers             The number of layers
!! @param[out] mass_divergence     The mass_divergence values in W3 space
!! @param[in]  undf_w3             The number of unique degrees of freedom
!! @param[in]  ndf_w3              The number of degrees of freedom per cell
!! @param[in]  map_w3              The dofmap for the cell at the base of the column
!! @param[in]  cell_orientation    The orientation of the cells
!! @param[in]  undf_w2             The number of unique degrees of freedom
!! @param[in]  ndf_w2              The number of degrees of freedom per cell
!! @param[in]  map_w2              The dofmap for the cell at the base of the column
!! @param[in]  mass_flux           The flux values which are calculated
!! @param[in]  direction           The direction in which to calculate the fluxes
subroutine fv_divergence_code( nlayers,              &
                               mass_divergence,      &
                               undf_w3,              &
                               ndf_w3,               &
                               map_w3,               &
                               cell_orientation,     &
                               undf_w2,              &
                               ndf_w2,               &
                               map_w2,               &
                               mass_flux,            &
                               direction )

  use cosmic_flux_mod, only : dof_to_update

  ! Arguments
  integer, intent(in)                                   :: nlayers
  integer, intent(in)                                   :: ndf_w3
  integer, intent(in)                                   :: undf_w3
  integer, dimension(ndf_w3), intent(in)                :: map_w3
  real(kind=r_def), dimension(undf_w3), intent(out)     :: mass_divergence
  real(kind=r_def), dimension(undf_w3), intent(in)      :: cell_orientation
  integer, intent(in)                                   :: ndf_w2
  integer, intent(in)                                   :: undf_w2
  integer, dimension(ndf_w2), intent(in)                :: map_w2
  real(kind=r_def), dimension(undf_w2), intent(in)      :: mass_flux
  integer, intent(in)                                   :: direction

  integer :: k
  integer :: local_dofs(1:2)

  local_dofs = dof_to_update(int(cell_orientation(map_w3(1))),direction)

  ! This kernel has been designed to work with lowest order W2 and W3 spaces.
  ! As is the case for all of the code associated with the Cosmic transport
  ! scheme.

  do k=0,nlayers-1

    mass_divergence( map_w3(1)+k ) = mass_flux(map_w2(local_dofs(2))+k) -     &
                                            mass_flux(map_w2(local_dofs(1))+k)

  end do

end subroutine fv_divergence_code

end module fv_divergence_kernel_mod
