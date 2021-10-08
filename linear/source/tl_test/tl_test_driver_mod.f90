!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief   Drives the execution of the tangent linear model tests.
!>@details The tests are initialised and finalised using a similar
!!         method to gungho, but with the addition of the linearisation state.
module tl_test_driver_mod

  use clock_mod,                  only : clock_type
  use constants_mod,              only : i_def, i_native, imdi
  use gungho_mod,                 only : program_name
  use gungho_model_mod,           only : initialise_infrastructure, &
                                         initialise_model,          &
                                         finalise_infrastructure,   &
                                         finalise_model
  use gungho_model_data_mod,      only : model_data_type,       &
                                         create_model_data,     &
                                         initialise_model_data, &
                                         finalise_model_data
  use io_context_mod,             only : io_context_type
  use log_mod,                    only : log_event,         &
                                         LOG_LEVEL_ALWAYS
  use linear_model_data_mod,      only : linear_create_ls,  &
                                         linear_init_ls
  use tl_test_kinetic_energy_gradient_mod, only : test_kinetic_energy_gradient
  use tl_test_advection_mod,               only : test_advection
  use tl_test_advect_density_field_mod,    only : test_advect_density_field
  use tl_test_advect_theta_field_mod,      only : test_advect_theta_field
  use tl_test_vorticity_mod,               only : test_vorticity_advection
  use tl_test_project_pressure_mod,        only : test_project_pressure
  use tl_test_hydrostatic_mod,             only : test_hydrostatic
  use tl_test_pressure_grad_bd_mod,        only : test_pressure_gradient_bd
  use tl_test_rk_alg_mod,                  only : test_rk_alg

  implicit none

  private
  public initialise,                  &
         finalise,                    &
         run_kinetic_energy_gradient, &
         run_advection,               &
         run_advect_density_field,    &
         run_advect_theta_field,      &
         run_vorticity_advection,     &
         run_project_pressure,        &
         run_hydrostatic,             &
         run_pressure_gradient_bd,    &
         run_rk_alg

  type (model_data_type) :: model_data

  integer(i_def) :: mesh_id              = imdi
  integer(i_def) :: twod_mesh_id         = imdi
  integer(i_def) :: shifted_mesh_id      = imdi
  integer(i_def) :: double_level_mesh_id = imdi

  class(io_context_type), allocatable :: io_context

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief     Sets up the required state in preparation for run.
  !>@param[in] filename            Name of the file containing the desired
  !!                               configuration
  !>@param[in] model_communicator  MPI communicator the model is to use
  subroutine initialise( filename, model_communicator )

    implicit none

    character(*),      intent(in) :: filename
    integer(i_native), intent(in) :: model_communicator

    class(clock_type), pointer :: clock

    ! Initialise infrastructure and setup constants
    call initialise_infrastructure( model_communicator,   &
                                    filename,             &
                                    program_name,         &
                                    io_context,           &
                                    mesh_id,              &
                                    twod_mesh_id,         &
                                    shifted_mesh_id,      &
                                    double_level_mesh_id, &
                                    model_data  )

    clock => io_context%get_clock()
    ! Instantiate the fields stored in model_data
    call create_model_data( model_data,   &
                            mesh_id,      &
                            twod_mesh_id, &
                            clock )

    ! Instantiate the linearisation state
    call linear_create_ls( model_data,    &
                           mesh_id,       &
                           twod_mesh_id )

    ! Initialise the fields stored in the model_data prognostics. This needs
    ! to be done before initialise_model.
    call initialise_model_data( model_data, clock )

    ! Model configuration initialisation
    call initialise_model( clock,         &
                           mesh_id,       &
                           model_data )

    ! Initialise the linearisation state
    call linear_init_ls( mesh_id,         &
                         twod_mesh_id,    &
                         model_data,      &
                         clock )

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tests the tangent linear kinetic energy gradient kernel.
  subroutine run_kinetic_energy_gradient()

    implicit none

    call test_kinetic_energy_gradient( model_data,  &
                                       mesh_id,     &
                                       twod_mesh_id )

  end subroutine run_kinetic_energy_gradient

  subroutine run_advection()

    implicit none

    call test_advection( model_data,  &
                         mesh_id,     &
                         twod_mesh_id )

  end subroutine run_advection

  subroutine run_advect_density_field()

    implicit none

    call test_advect_density_field( model_data,  &
                                    mesh_id,     &
                                    twod_mesh_id )

  end subroutine run_advect_density_field

  subroutine run_advect_theta_field()

    implicit none

    call test_advect_theta_field( model_data,  &
                                  mesh_id,     &
                                  twod_mesh_id )

  end subroutine run_advect_theta_field

  subroutine run_vorticity_advection()

    implicit none

    call test_vorticity_advection( model_data,  &
                                   mesh_id,     &
                                   twod_mesh_id )

  end subroutine run_vorticity_advection

  subroutine run_project_pressure()

    implicit none

    call test_project_pressure( model_data,  &
                                mesh_id,     &
                                twod_mesh_id )

  end subroutine run_project_pressure

  subroutine run_hydrostatic()

    implicit none

    call test_hydrostatic( model_data,  &
                           mesh_id,     &
                           twod_mesh_id )

  end subroutine run_hydrostatic

  subroutine run_pressure_gradient_bd()

    implicit none

    call test_pressure_gradient_bd( model_data,  &
                                    mesh_id,     &
                                    twod_mesh_id )

  end subroutine run_pressure_gradient_bd

  subroutine run_rk_alg()

    implicit none

    call test_rk_alg( model_data,  &
                      mesh_id,     &
                      twod_mesh_id )

  end subroutine run_rk_alg

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tidies up after a run.
  subroutine finalise()

    implicit none

    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Model configuration finalisation
    call finalise_model( mesh_id,    &
                         model_data, &
                         program_name )

    ! Destroy the fields stored in model_data
    call finalise_model_data( model_data )

    ! Finalise infrastructure and constants
    call finalise_infrastructure( program_name )

  end subroutine finalise

end module tl_test_driver_mod
