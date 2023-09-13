!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Custom extrusion selector for jedi_lfric
!> @details Allows um_L70_50t_20s_80km extrusion in the miniapp without
!>          importing GungHo. um_L70_50t_20s_80km in this module is a copy of
!>          the extrusion in gungho/source/mesh/gungho_extrusion_mod.f90 of
!>          the same name.
module jedi_lfric_fake_nl_extrusion_mod
  use constants_mod,        only: r_def, i_def
  use extrusion_mod,        only: extrusion_type,               &
                                  uniform_extrusion_type,       &
                                  geometric_extrusion_type,     &
                                  quadratic_extrusion_type,     &
                                  shifted_extrusion_type,       &
                                  double_level_extrusion_type,  &
                                  PRIME_EXTRUSION
  use base_mesh_config_mod, only: geometry,                     &
                                  geometry_spherical,           &
                                  geometry_planar
  use extrusion_config_mod, only: method,                       &
                                  method_uniform,               &
                                  method_geometric,             &
                                  method_quadratic,             &
                                  method_um_L70_50t_20s_80km,   &
                                  domain_top, number_of_layers
  use planet_config_mod,    only: scaled_radius
  use log_mod,              only: log_event,                    &
                                  LOG_LEVEL_ERROR

  implicit none

  !> @brief Extrudes with specific UM configuration L70_50t_20s_80km
  !>
  type, public, extends(extrusion_type) :: um_L70_50t_20s_80km_extrusion_type
    private
  contains
    private
    procedure, public :: extrude => um_L70_50t_20s_80km_extrude
  end type um_L70_50t_20s_80km_extrusion_type

  interface um_L70_50t_20s_80km_extrusion_type
    module procedure um_L70_50t_20s_80km_extrusion_constructor
  end interface um_L70_50t_20s_80km_extrusion_type

contains

  !> @brief Creates a um_L70_50t_20s_80km_extrusion_type object.
  !>
  !> @param[in] atmosphere_bottom Bottom of the atmosphere in meters.
  !> @param[in] atmosphere_top Top of the atmosphere in meters.
  !> @param[in] number_of_layers Number of layers in the atmosphere.
  !> @param[in] extrusion_id Identifier of extrusion type.
  !>
  !> @return New uniform_extrusion_type object.
  function um_L70_50t_20s_80km_extrusion_constructor( atmosphere_bottom, &
                                                      atmosphere_top,    &
                                                      number_of_layers,  &
                                                      extrusion_id ) result(new)

    implicit none

    real(r_def),    intent(in) :: atmosphere_bottom
    real(r_def),    intent(in) :: atmosphere_top
    integer(i_def), intent(in) :: number_of_layers
    integer(i_def), intent(in) :: extrusion_id

    type(um_L70_50t_20s_80km_extrusion_type) :: new

    call new%extrusion_constructor( atmosphere_bottom, atmosphere_top, &
                                    number_of_layers, extrusion_id )

  end function um_L70_50t_20s_80km_extrusion_constructor

  !> @brief Extrudes the mesh with specific UM configuration L70_50t_20s_80km
  !>
  !> @param[out] eta Nondimensional vertical coordinate.
  subroutine um_L70_50t_20s_80km_extrude( this, eta )

    implicit none

    class(um_L70_50t_20s_80km_extrusion_type), intent(in)  :: this
    real(r_def),                               intent(out) :: eta(0:)


    if (this%get_number_of_layers() /= 70) then
      call log_event( "Extrusion L70_50t_20s_80km reqires 70 levels", log_level_error )
    end if

    eta(0:this%get_number_of_layers()) = (/  &
      0.0000000_r_def,  0.0002500_r_def,  0.0006667_r_def,  0.0012500_r_def, &
      0.0020000_r_def,  0.0029167_r_def,  0.0040000_r_def,  0.0052500_r_def, &
      0.0066667_r_def,  0.0082500_r_def,  0.0100000_r_def,  0.0119167_r_def, &
      0.0140000_r_def,  0.0162500_r_def,  0.0186667_r_def,  0.0212500_r_def, &
      0.0240000_r_def,  0.0269167_r_def,  0.0300000_r_def,  0.0332500_r_def, &
      0.0366667_r_def,  0.0402500_r_def,  0.0440000_r_def,  0.0479167_r_def, &
      0.0520000_r_def,  0.0562500_r_def,  0.0606667_r_def,  0.0652500_r_def, &
      0.0700000_r_def,  0.0749167_r_def,  0.0800000_r_def,  0.0852500_r_def, &
      0.0906668_r_def,  0.0962505_r_def,  0.1020017_r_def,  0.1079213_r_def, &
      0.1140113_r_def,  0.1202745_r_def,  0.1267154_r_def,  0.1333406_r_def, &
      0.1401592_r_def,  0.1471838_r_def,  0.1544313_r_def,  0.1619238_r_def, &
      0.1696895_r_def,  0.1777643_r_def,  0.1861929_r_def,  0.1950307_r_def, &
      0.2043451_r_def,  0.2142178_r_def,  0.2247466_r_def,  0.2360480_r_def, &
      0.2482597_r_def,  0.2615432_r_def,  0.2760868_r_def,  0.2921094_r_def, &
      0.3098631_r_def,  0.3296378_r_def,  0.3517651_r_def,  0.3766222_r_def, &
      0.4046373_r_def,  0.4362943_r_def,  0.4721379_r_def,  0.5127798_r_def, &
      0.5589045_r_def,  0.6112759_r_def,  0.6707432_r_def,  0.7382500_r_def, &
      0.8148403_r_def,  0.9016668_r_def,  1.0000000_r_def /)

  end subroutine um_L70_50t_20s_80km_extrude

  !> @brief  Creates the prime vertical mesh extrusion.
  !> @return  Resulting extrusion object
  function create_extrusion() result(new)

    implicit none

    class(extrusion_type), allocatable :: new

    real(kind=r_def) :: domain_bottom

    if (allocated(new)) deallocate(new)

    select case (geometry)
      case (geometry_planar)
        domain_bottom = 0.0_r_def
      case (geometry_spherical)
        domain_bottom = scaled_radius
      case default
        call log_event("Invalid geometry for mesh initialisation", LOG_LEVEL_ERROR)
    end select

    select case (method)
      case (method_uniform)
        allocate( new, source=uniform_extrusion_type( domain_bottom,    &
                                                      domain_top,       &
                                                      number_of_layers, &
                                                      PRIME_EXTRUSION ) )
      case (method_quadratic)
        allocate( new, source=quadratic_extrusion_type( domain_bottom,    &
                                                        domain_top,       &
                                                        number_of_layers, &
                                                        PRIME_EXTRUSION ) )
      case (method_geometric)
        allocate( new, source=geometric_extrusion_type( domain_bottom,    &
                                                        domain_top,       &
                                                        number_of_layers, &
                                                        PRIME_EXTRUSION ) )
      case (method_um_L70_50t_20s_80km)
        allocate( new, source=um_L70_50t_20s_80km_extrusion_type( domain_bottom,    &
                                                                  domain_top,       &
                                                                  number_of_layers, &
                                                                  PRIME_EXTRUSION ) )
      case default
        call log_event("Invalid method for simple extrusion", LOG_LEVEL_ERROR)
    end select

  end function create_extrusion

end module jedi_lfric_fake_nl_extrusion_mod
