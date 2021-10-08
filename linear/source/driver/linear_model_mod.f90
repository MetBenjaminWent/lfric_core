!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>@brief Initialises and finalises the linear model numerical schemes.
module linear_model_mod

  use constants_mod,              only : i_def
  use clock_mod,                  only : clock_type
  use field_mod,                  only : field_type
  use field_collection_mod,       only : field_collection_type
  use gungho_model_data_mod,      only : model_data_type
  use tl_rk_alg_timestep_mod,     only : tl_rk_alg_init, &
                                         tl_rk_alg_final
  use formulation_config_mod,     only : transport_only
  use timestepping_config_mod,    only : method,               &
                                         method_semi_implicit, &
                                         method_rk
  use log_mod,                    only : log_event, &
                                         log_scratch_space, &
                                         LOG_LEVEL_INFO,    &
                                         LOG_LEVEL_TRACE,   &
                                         LOG_LEVEL_ERROR
  implicit none

  private
  public initialise_linear_model, &
         finalise_linear_model

contains

  !> @brief Completes the initialisation of the tangent linear model
  !> @param[in] clock Model time
  !> @param[in] mesh_id The identifier of the primary mesh
  !> @param[in,out] model_data The working data set for the model run
  subroutine initialise_linear_model( clock,   &
                                      mesh_id, &
                                      model_data )
    implicit none

    class(clock_type),       intent(in), pointer   :: clock
    integer(i_def),          intent(in)            :: mesh_id
    type( model_data_type ), intent(inout), target :: model_data

    type( field_collection_type ), pointer :: prognostic_fields => null()
    type( field_type ),            pointer :: mr(:) => null()
    type( field_collection_type ), pointer :: ls_fields => null()
    type( field_type ),            pointer :: ls_mr(:) => null()

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()
    type( field_type), pointer :: ls_theta => null()
    type( field_type), pointer :: ls_u => null()
    type( field_type), pointer :: ls_rho => null()
    type( field_type), pointer :: ls_exner => null()

    ! Get pointers to field collections for use downstream
    prognostic_fields => model_data%prognostic_fields
    mr => model_data%mr
    ls_fields => model_data%ls_fields
    ls_mr => model_data%ls_mr

    ! Get pointers to fields in the prognostic/diagnostic field collections
    ! for use downstream
    theta => prognostic_fields%get_field('theta')
    u => prognostic_fields%get_field('u')
    rho => prognostic_fields%get_field('rho')
    exner => prognostic_fields%get_field('exner')
    ls_theta => ls_fields%get_field('ls_theta')
    ls_u => ls_fields%get_field('ls_u')
    ls_rho => ls_fields%get_field('ls_rho')
    ls_exner => ls_fields%get_field('ls_exner')

    if ( transport_only ) then
      call log_event("TL Transport only not available ",LOG_LEVEL_ERROR)

    else
      select case( method )
        case( method_semi_implicit )  ! Semi-Implicit
          call log_event("TL Semi implicit only not available ",LOG_LEVEL_ERROR)

        case( method_rk )             ! RK
          ! Initialise and output initial conditions for first timestep

          call tl_rk_alg_init(mesh_id, u, rho, theta, exner, &
                              ls_u, ls_rho, ls_theta, ls_exner)

        case default
          call log_event("TL: Incorrect time stepping option chosen, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
      end select

    end if

  end subroutine initialise_linear_model

  !> @brief Finalises the remaining infrastructure and constants used by the linear model
  subroutine finalise_linear_model

    implicit none

    if ( .not. transport_only ) then
      if ( method == method_rk )            call tl_rk_alg_final()
    end if

  end subroutine finalise_linear_model

end module linear_model_mod
