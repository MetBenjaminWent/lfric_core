!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing the JEDI Geometry emulator class.
!>
!> @details This module holds a JEDI Geometry emulator class that includes only
!>           the functionality that is required to support interface emulation.
!>
module jedi_geometry_mod

  use, intrinsic :: iso_fortran_env, only : real64

  use constants_mod,                 only : i_def
  use jedi_lfric_mesh_interface_mod, only : set_target_mesh,           &
                                            is_mesh_cubesphere,        &
                                            get_domain_top,            &
                                            get_cubesphere_resolution, &
                                            get_nlayers,               &
                                            get_layer_ncells,          &
                                            get_lonlat,                &
                                            get_sigma_w3_levels,       &
                                            get_sigma_wtheta_levels,   &
                                            get_stretching_height
  use log_mod,                       only : log_event, LOG_LEVEL_ERROR

  implicit none

  private

type, public :: jedi_geometry_type
  private

  !> The data map between external field data and LFRic fields
  integer( kind = i_def ), allocatable :: horizontal_map(:)
  !> the LFRic field dimensions
  integer( kind = i_def )              :: n_layers
  integer( kind = i_def )              :: n_horizontal

contains

  !> Field initialiser.
  procedure, public :: initialise

  !> getters
  procedure, public :: get_n_horizontal
  procedure, public :: get_n_layers
  procedure, public :: get_horizontal_map

  !> Finalizer
  final             :: jedi_geometry_destructor

end type jedi_geometry_type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @brief    Initialiser for jedi_geometry_type
!>
subroutine initialise( self )
  ! Access config directly until modeldb ready
  use base_mesh_config_mod, only: prime_mesh_name

  implicit none

  class( jedi_geometry_type ), intent(inout) :: self

  ! Local
  integer :: i_horizontal
  real(real64) :: domain_top, stretching_height
  real(real64), allocatable :: lonlat(:,:),          &
                               sigma_W3_levels(:),   &
                               sigma_Wtheta_levels(:)

  ! Set target mesh for all functions in the interface
  call set_target_mesh( prime_mesh_name )

  if ( .not. is_mesh_cubesphere() ) then
    call log_event( "Working mesh is not a cubesphere", LOG_LEVEL_ERROR )
  end if

  ! Get grid size and layers
  self%n_layers = get_nlayers()
  self%n_horizontal = get_layer_ncells()

  ! Create horizontal_map
  lonlat = get_lonlat()
  allocate( self%horizontal_map( self%n_horizontal ) )

  ! For mock purposes return sequential map
  do i_horizontal=1,self%n_horizontal
    self%horizontal_map( i_horizontal ) = i_horizontal
  end do

  ! Here JEDI deals with physical coordinates
  domain_top = get_domain_top()
  sigma_W3_levels = get_sigma_w3_levels()
  sigma_Wtheta_levels = get_sigma_wtheta_levels()
  stretching_height = get_stretching_height()

end subroutine initialise

!> @brief    Get the number of horizontal points
!>
!> @return n_horizontal The number of horizontal points
function get_n_horizontal(self) result(n_horizontal)

  implicit none

  class( jedi_geometry_type ), intent(in) :: self
  integer(kind = i_def) :: n_horizontal

  n_horizontal = self%n_horizontal

end function get_n_horizontal

!> @brief    Get the number of model layers
!>
!> @return n_layers The number of model layers
function get_n_layers(self) result(n_layers)

  implicit none

  class( jedi_geometry_type ), intent(in) :: self
  integer(kind = i_def)              :: n_layers

  n_layers = self%n_layers

end function get_n_layers

!> @brief    Get a pointer to the horizontal map
!>
!> @return horizontal_map A pointer to the map providing the horizontal index
!>                        of the Atlas field
subroutine get_horizontal_map(self, horizontal_map)

  implicit none

  class( jedi_geometry_type ), target, intent(in) :: self
  integer(i_def), pointer, intent(inout)     :: horizontal_map(:)

  horizontal_map => self % horizontal_map

end subroutine get_horizontal_map

!> @brief    Finalizer for jedi_geometry_type
!>
subroutine jedi_geometry_destructor(self)

  implicit none

  type(jedi_geometry_type), intent(inout)    :: self

  if ( allocated(self % horizontal_map) ) deallocate(self % horizontal_map)

end subroutine jedi_geometry_destructor

end module jedi_geometry_mod
