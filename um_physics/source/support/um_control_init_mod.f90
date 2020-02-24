!----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Controls the setting of UM high level variables which are either
!>         fixed in LFRic or derived from LFRic inputs.

module um_control_init_mod

  ! LFRic namelists which have been read
  use extrusion_config_mod,    only : number_of_layers
  use timestepping_config_mod, only : dt, outer_iterations

  ! Other modules used
  use constants_mod,           only : i_um, r_um, rmdi
  use conversions_mod,         only : pi_over_180

  implicit none

  private
  public :: um_control_init

contains

  !>@brief Initialise UM high levels variables with are either fixed in LFRic
  !>        or derived from LFRic inputs.
  !>@details Nothing in this file is ever likely to be promoted to the LFRic
  !>          namelist. Everything is either set from an LFRic variable
  !>          already in the namelist, or is "fixed" from the perspective
  !>          of LFRic (but cannot be made a parameter because it is required
  !>          to be variable in the UM and therefore declared as such in the
  !>          UM modules which contains it).
  subroutine um_control_init()

    ! UM modules containing things that need setting
    use atm_fields_bounds_mod, only: atm_fields_bounds_init
    use atm_step_local, only: rhc_row_length, rhc_rows
    use cderived_mod, only: delta_lambda, delta_phi
    use dynamics_input_mod, only: numcycles
    use gen_phys_inputs_mod, only: l_mr_physics
    use level_heights_mod, only: r_theta_levels, r_rho_levels, eta_theta_levels
    use model_domain_mod, only: model_type, mt_single_column
    use nlsizes_namelist_mod, only: model_levels, cloud_levels, n_cca_lev, &
         bl_levels, row_length, rows, tr_vars
    use timestep_mod, only: timestep
    use trignometric_mod, only: cos_theta_latitude

    implicit none

    ! ----------------------------------------------------------------
    ! Model dimensions - contained in UM module nlsizes_namelist_mod
    ! ----------------------------------------------------------------
    ! Vertical dimensions set from LFRic number of layers.
    model_levels = number_of_layers
    cloud_levels = number_of_layers
    n_cca_lev    = number_of_layers
    bl_levels    = number_of_layers-1
    ! Horizontal dimensions set to 1 to match the number of cells passed
    ! to physics kernels.
    ! This may need to be variable if we pass multiple cells to kernels.
    row_length = 1
    rows       = 1
    ! Number of tracer variables - not implemented in LFRic so fixed at 0.
    ! This may change if we implement tracers ever.
    tr_vars = 0

    ! ----------------------------------------------------------------
    ! More model dimensions, this time from atm_step_local
    ! ----------------------------------------------------------------
    ! Dimensions of critical relative humidity array. Again, set to 1 to
    ! match number of points passed to kernels, but may change...
    rhc_row_length = 1
    rhc_rows       = 1

    ! ----------------------------------------------------------------
    ! Model type - contained in UM module model_domain_mod
    ! ----------------------------------------------------------------
    ! Currently set to single column, as the physics kernels are called 1
    ! cell at a time. This may change if we pass multiple cells to kernels.
    model_type = mt_single_column

    ! Set the field bounds which are used by the UM code based on the
    ! information above.
    ! Hard-wired zeros are halo-sizes in UM code.
    ! Currently set to zero as we don't pass a stencil into the kernels
    ! but may change if we ever do.
    call atm_fields_bounds_init( 0_i_um, 0_i_um, 0_i_um, &
                                 0_i_um, row_length, rows, rows, &
                                 bl_levels_opt = bl_levels )

    ! Timestep used in UM code - contained in UM timestep_mod.
    ! Set from LFRic input timestep.
    timestep = real(dt, r_um)

    ! Number of outer iterations of the dynamics - contained in UM
    ! dynamics_input_mod.
    ! Set from LFRic input number of outer iterations.
    numcycles = outer_iterations

    ! Mixing ratio flag - contained in UM gen_phys_inputs.
    ! Set to true as LFRic only supports mixing ratios.
    l_mr_physics = .true.

    ! The following 3D arrays are used direct from level_heights_mod
    ! throughout the UM code.
    ! We must initialise them here so that they are always available.
    ! But they must be set to appropriate values for the current column
    ! in any kernel whos external code uses the variables.
    ! Ideally the UM code will be changed so that they are passed in
    ! through the argument list.
    allocate(r_theta_levels(row_length,rows,0:number_of_layers), source=rmdi)
    allocate(r_rho_levels(row_length,rows,number_of_layers), source=rmdi)
    allocate(eta_theta_levels(0:number_of_layers), source=rmdi)

    ! The following are used in the calculation of grid-box size in
    ! UM parametrizations.
    ! They are temporarily fixed here to be representative of an N96/C48 model.
    ! They should ultimately be set from the LFRic mesh directly.
    allocate(cos_theta_latitude(row_length,rows), source=1.0_r_um)
    delta_lambda = 1.875_r_um*pi_over_180
    delta_phi    = 1.875_r_um*pi_over_180

  end subroutine um_control_init

end module um_control_init_mod
