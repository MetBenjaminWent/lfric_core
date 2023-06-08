!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>
!> @brief Drives the execution of the transport miniapp.
!>
module transport_driver_mod

  use base_mesh_config_mod,             only: prime_mesh_name
  use checksum_alg_mod,                 only: checksum_alg
  use check_configuration_mod,          only: get_required_stencil_depth
  use configuration_mod,                only: final_configuration
  use constants_mod,                    only: i_def, i_native, r_def, &
                                              r_second, str_def
  use driver_fem_mod,                   only: init_fem
  use driver_io_mod,                    only: init_io, final_io
  use driver_mesh_mod,                  only: init_mesh
  use driver_time_mod,                  only: init_time, get_calendar
  use derived_config_mod,               only: set_derived_config
  use diagnostics_io_mod,               only: write_scalar_diagnostic, &
                                              write_vector_diagnostic
  use divergence_alg_mod,               only: divergence_alg
  use field_mod,                        only: field_type
  use formulation_config_mod,           only: l_multigrid
  use geometric_constants_mod,          only: get_chi_inventory, &
                                              get_panel_id_inventory
  use io_context_mod,                   only: io_context_type
  use io_config_mod,                    only: diagnostic_frequency,            &
                                              nodal_output_on_w3,              &
                                              write_diag,                      &
                                              use_xios_io,                     &
                                              subroutine_timers,               &
                                              timer_output_path
  use inventory_by_mesh_mod,            only: inventory_by_mesh_type
  use local_mesh_mod,                   only: local_mesh_type
  use log_mod,                          only: log_event,         &
                                              log_scratch_space, &
                                              LOG_LEVEL_ALWAYS,  &
                                              LOG_LEVEL_INFO,    &
                                              LOG_LEVEL_TRACE
  use mass_conservation_alg_mod,        only: mass_conservation
  use mesh_mod,                         only: mesh_type
  use mesh_collection_mod,              only: mesh_collection
  use model_clock_mod,                  only: model_clock_type
  use mpi_mod,                          only: mpi_type
  use mr_indices_mod,                   only: nummr
  use runtime_constants_mod,            only: create_runtime_constants
  use step_calendar_mod,                only: step_calendar_type
  use timer_mod,                        only: init_timer, timer, output_timer
  use timestepping_config_mod,          only: dt
  use transport_init_fields_alg_mod,    only: transport_init_fields_alg
  use transport_control_alg_mod,        only: transport_prerun_setup, &
                                              transport_init, &
                                              transport_step, &
                                              transport_final
  use transport_runtime_collection_mod, only: init_transport_runtime_collection, &
                                              transport_runtime_collection_final

  implicit none

  private

  public :: initialise_transport, run_transport, finalise_transport

  type(model_clock_type), allocatable :: model_clock

  ! Prognostic fields
  type(field_type) :: wind
  type(field_type) :: density
  type(field_type) :: theta
  type(field_type) :: tracer_con
  type(field_type) :: tracer_adv
  type(field_type) :: constant
  type(field_type) :: mr(nummr)
  type(field_type) :: divergence

  ! Number of moisutre species to transport
  integer(kind=i_def) :: nummr_to_transport

  type(mesh_type), pointer :: mesh => null()

  integer(i_def) :: num_meshes


contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up required state in preparation for run.
  !> @param [in,out] mpi          The structure that holds comms details
  !> @param [in]     program_name An identifier given to the model being run
  subroutine initialise_transport( mpi, program_name )

    implicit none

    class(mpi_type), intent(inout) :: mpi
    character(*),    intent(in)    :: program_name

    character(len=*), parameter :: xios_ctx  = "transport"


    integer(kind=i_def),      allocatable :: local_mesh_ids(:)
    type(local_mesh_type),        pointer :: local_mesh => null()
    type(mesh_type),              pointer :: mesh => null()
    type(inventory_by_mesh_type), pointer :: chi_inventory => null()
    type(inventory_by_mesh_type), pointer :: panel_id_inventory => null()
    character(str_def),       allocatable :: base_mesh_names(:)
    character(str_def),       allocatable :: shifted_mesh_names(:)
    character(str_def),       allocatable :: double_level_mesh_names(:)

    call log_event( program_name//': Runtime default precision set as:', LOG_LEVEL_ALWAYS )
    write(log_scratch_space, '(I1)') kind(1.0_r_def)
    call log_event( '     r_def kind = '//log_scratch_space , LOG_LEVEL_ALWAYS )
    write(log_scratch_space, '(I1)') kind(1_i_def)
    call log_event( '     i_def kind = '//log_scratch_space , LOG_LEVEL_ALWAYS )

    call set_derived_config( .true. )

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------
    if ( subroutine_timers ) then
      call init_timer(timer_output_path)
      call timer( program_name )
    end if

    call init_time( model_clock )

    !-------------------------------------------------------------------------
    ! Work out which meshes are required
    !-------------------------------------------------------------------------

    ! Always have just one base mesh with shifted version
    allocate(base_mesh_names(1))
    allocate(shifted_mesh_names(1))
    allocate(double_level_mesh_names(1))
    base_mesh_names(1) = prime_mesh_name
    shifted_mesh_names(1) = prime_mesh_name
    double_level_mesh_names(1) = prime_mesh_name

    !-------------------------------------------------------------------------
    ! Initialise mesh and FEM infrastructure
    !-------------------------------------------------------------------------

    ! Create the mesh
    call init_mesh( mpi%get_comm_rank(), mpi%get_comm_size(),    &
                    base_mesh_names,                                        &
                    shifted_mesh_names=shifted_mesh_names,                  &
                    double_level_mesh_names=double_level_mesh_names,        &
                    required_stencil_depth=get_required_stencil_depth() )

    ! FEM initialisation
    chi_inventory => get_chi_inventory()
    panel_id_inventory => get_panel_id_inventory()

    call init_fem( mesh_collection, chi_inventory, panel_id_inventory )

    ! Create runtime_constants object.
    call create_runtime_constants( mesh_collection,    &
                                   chi_inventory,      &
                                   panel_id_inventory, &
                                   model_clock )

    ! Set up transport runtime collection type
    ! Transport on only one horizontal local mesh
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    local_mesh => mesh%get_local_mesh()
    allocate(local_mesh_ids(1))
    local_mesh_ids(1) = local_mesh%get_id()
    call init_transport_runtime_collection(local_mesh_ids)

    ! Set transport metadata for primal mesh
    num_meshes = 1_i_def
    call transport_prerun_setup( num_meshes )

    ! Initialise prognostic variables
    call transport_init_fields_alg( mesh, wind, density, theta, &
                                    tracer_con, tracer_adv,     &
                                    constant, mr, divergence )

    ! Initialise all transport-only control algorithm
    call transport_init( density, theta, tracer_con, tracer_adv, constant, mr )

    nummr_to_transport = 1_i_def

    ! I/O initialisation
    call init_io( xios_ctx,           &
                  mpi%get_comm(), &
                  chi_inventory,      &
                  panel_id_inventory, &
                  model_clock,        &
                  get_calendar() )

    ! Output initial conditions
    if (model_clock%is_initialisation() .and. write_diag) then

      call write_vector_diagnostic( 'u', wind, model_clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'rho', density, model_clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'theta', theta, model_clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'tracer_con', tracer_con, model_clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'tracer_adv', tracer_adv, model_clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'constant', constant, model_clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'm_v', mr(1), model_clock, &
                                    mesh, nodal_output_on_w3 )
      call write_scalar_diagnostic( 'divergence', divergence, model_clock, &
                                    mesh, nodal_output_on_w3 )

    end if

    deallocate(base_mesh_names)
    deallocate(shifted_mesh_names)
    deallocate(double_level_mesh_names)
    deallocate(local_mesh_ids)
    nullify(chi_inventory, panel_id_inventory, mesh, local_mesh)

  end subroutine initialise_transport

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Performs time stepping.
  !>
  subroutine run_transport()

    implicit none

    call log_event( 'Miniapp will run with default precision set as:', LOG_LEVEL_INFO )
    write(log_scratch_space, '(I1)') kind(1.0_r_def)
    call log_event( '        r_def kind = '//log_scratch_space , LOG_LEVEL_INFO )
    write(log_scratch_space, '(I1)') kind(1_i_def)
    call log_event( '        i_def kind = '//log_scratch_space , LOG_LEVEL_INFO )

    call mass_conservation( model_clock%get_step(), density, mr )
    call density%log_minmax( LOG_LEVEL_INFO, 'rho' )
    call theta%log_minmax( LOG_LEVEL_INFO, 'theta' )
    call tracer_con%log_minmax( LOG_LEVEL_INFO, 'tracer_con' )
    call tracer_adv%log_minmax( LOG_LEVEL_INFO, 'tracer_adv' )
    call constant%log_minmax( LOG_LEVEL_INFO, 'constant' )
    call mr(1)%log_minmax( LOG_LEVEL_INFO, 'm_v' )

    mesh => mesh_collection%get_mesh(prime_mesh_name)

    !--------------------------------------------------------------------------
    ! Model step
    !--------------------------------------------------------------------------
    do while (model_clock%tick())

      write(log_scratch_space, '("/", A, "\ ")') repeat('*', 76)
      call log_event( log_scratch_space, LOG_LEVEL_TRACE )
      write( log_scratch_space, '(A,I0)' ) &
        'Start of timestep ', model_clock%get_step()
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      if ( subroutine_timers ) call timer( 'transport step' )

      call transport_step( model_clock,                      &
                           wind, density, theta, tracer_con, &
                           tracer_adv, constant, mr,         &
                           nummr_to_transport )

      if ( subroutine_timers ) call timer( 'transport step' )

      ! Write out conservation diagnostics
      call mass_conservation( model_clock%get_step(), density, mr )
      call density%log_minmax( LOG_LEVEL_INFO, 'rho' )
      call theta%log_minmax( LOG_LEVEL_INFO, 'theta' )
      call tracer_con%log_minmax( LOG_LEVEL_INFO, 'tracer_con' )
      call tracer_adv%log_minmax( LOG_LEVEL_INFO, 'tracer_adv' )
      call constant%log_minmax( LOG_LEVEL_INFO, 'constant' )
      call mr(1)%log_minmax( LOG_LEVEL_INFO, 'm_v' )


      write( log_scratch_space, &
             '(A,I0)' ) 'End of timestep ', model_clock%get_step()
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
      write(log_scratch_space, '("\", A, "/ ")') repeat('*', 76)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      ! Output wind and density values.
      if ( (mod( model_clock%get_step(), diagnostic_frequency ) == 0) &
           .and. write_diag ) then

        ! Compute divergence
        call divergence_alg( divergence, wind )

        call write_vector_diagnostic( 'u', wind,                &
                                      model_clock, mesh, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'rho', density,           &
                                      model_clock, mesh, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'theta', theta,           &
                                      model_clock, mesh, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'tracer_con', tracer_con,         &
                                      model_clock, mesh, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'tracer_adv', tracer_adv,         &
                                      model_clock, mesh, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'constant', constant,         &
                                      model_clock, mesh, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'm_v', mr(1),             &
                                      model_clock, mesh, nodal_output_on_w3 )
        call write_scalar_diagnostic( 'divergence', divergence, &
                                      model_clock, mesh, nodal_output_on_w3 )


      end if

    end do ! while clock%is_running()

    call transport_final( density, theta, tracer_con, tracer_adv, constant, mr )

    nullify(mesh)

  end subroutine run_transport

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Tidies up after a run.
  !>
  subroutine finalise_transport( program_name )

    implicit none

    character(*), intent(in) :: program_name

    !--------------------------------------------------------------------------
    ! Model finalise
    !--------------------------------------------------------------------------

    ! Write checksums to file
    call checksum_alg( program_name, density, 'rho',  wind, 'u',  &
                       theta, 'theta', tracer_adv, 'tracer',      &
                       field_bundle=mr, bundle_name='mr' )

    if ( subroutine_timers ) then
      call timer( program_name )
      call output_timer()
    end if

    call final_io()

    call transport_runtime_collection_final()

  end subroutine finalise_transport

end module transport_driver_mod
