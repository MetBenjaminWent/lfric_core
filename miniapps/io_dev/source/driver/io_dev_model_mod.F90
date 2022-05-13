!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Handles initialisation and finalisation of infrastructure, constants
!>        and the io_dev model.
!>
module io_dev_model_mod

  ! Infrastructure
  use base_mesh_config_mod,       only : prime_mesh_name
  use constants_mod,              only : i_def, i_native, &
                                         PRECISION_REAL
  use convert_to_upper_mod,       only : convert_to_upper
  use field_mod,                  only : field_type
  use halo_comms_mod,             only : initialise_halo_comms, &
                                         finalise_halo_comms
  use local_mesh_collection_mod,  only : local_mesh_collection, &
                                         local_mesh_collection_type
  use log_mod,                    only : log_event,          &
                                         log_set_level,      &
                                         log_scratch_space,  &
                                         initialise_logging, &
                                         finalise_logging,   &
                                         LOG_LEVEL_ALWAYS,   &
                                         LOG_LEVEL_ERROR,    &
                                         LOG_LEVEL_WARNING,  &
                                         LOG_LEVEL_INFO,     &
                                         LOG_LEVEL_DEBUG,    &
                                         LOG_LEVEL_TRACE
  use mesh_collection_mod,        only : mesh_collection, &
                                         mesh_collection_type
  use mesh_mod,                   only : mesh_type
  use mpi_mod,                    only : store_comm,    &
                                         get_comm_size, &
                                         get_comm_rank
  use timer_mod,                  only : timer, output_timer, init_timer
  ! Configuration
  use configuration_mod,          only : final_configuration
  use io_config_mod,              only : subroutine_timers, &
                                         timer_output_path
  ! IO_Dev driver modules
  use io_dev_mod,                 only : load_configuration
  use io_dev_init_files_mod,      only : init_io_dev_files
  ! Driver modules
  use driver_fem_mod,             only : init_fem, final_fem
  use driver_mesh_mod,            only : init_mesh, final_mesh
  use driver_io_mod,              only : init_io, final_io, &
                                         filelist_populator
  ! External libraries
  use xios,                       only : xios_context_finalize

  implicit none

  private
  public initialise_infrastructure, &
         finalise_infrastructure

contains

  !> @brief Initialises the infrastructure components of the model.
  !>
  !> @param[in]     filename     The name of the configuration namelist file
  !> @param[in]     program_name An identifier given to the model run
  !> @param[in]     communicator The MPI communicator for use within the model
  !>                              (not XIOS' communicator)
  !> @param[in,out] mesh         The current 3d mesh
  !> @param[in,out] twod_mesh    The current 2d mesh
  !> @param[in,out] chi          A size 3 array of fields holding the
  !>                             coordinates of the mesh
  !> @param[in,out] panel_id     A 2D field holding the cubed sphere panel id
  !>
  subroutine initialise_infrastructure( filename,     &
                                        program_name, &
                                        communicator, &
                                        mesh,         &
                                        twod_mesh,    &
                                        chi,          &
                                        panel_id )

    use logging_config_mod, only: run_log_level,          &
                                  key_from_run_log_level, &
                                  RUN_LOG_LEVEL_ERROR,    &
                                  RUN_LOG_LEVEL_INFO,     &
                                  RUN_LOG_LEVEL_DEBUG,    &
                                  RUN_LOG_LEVEL_TRACE,    &
                                  RUN_LOG_LEVEL_WARNING

    implicit none

    ! Arguments
    character(*),           intent(in)    :: filename
    character(*),           intent(in)    :: program_name
    integer(i_native),      intent(in)    :: communicator

    type(mesh_type),        intent(inout), pointer :: mesh
    type(mesh_type),        intent(inout), pointer :: twod_mesh

    type(field_type),       intent(inout) :: chi(3)
    type(field_type),       intent(inout) :: panel_id

    ! Local variables
    character(*), parameter :: xios_context_id = 'io_dev'

    procedure(filelist_populator), pointer :: files_init_ptr => null()

    integer(i_def)    :: total_ranks, local_rank, stencil_depth
    integer(i_native) :: log_level

    ! Save the model's part of the split communicator for later use
    call store_comm( communicator )

    ! Initialise halo functionality
    call initialise_halo_comms( communicator )

    total_ranks = get_comm_size()
    local_rank  = get_comm_rank()

    call initialise_logging( local_rank, total_ranks, program_name )

    call load_configuration( filename )

    select case (run_log_level)
    case( RUN_LOG_LEVEL_ERROR )
      log_level = LOG_LEVEL_ERROR
    case( RUN_LOG_LEVEL_WARNING )
      log_level = LOG_LEVEL_WARNING
    case( RUN_LOG_LEVEL_INFO )
      log_level = LOG_LEVEL_INFO
    case( RUN_LOG_LEVEL_DEBUG )
      log_level = LOG_LEVEL_DEBUG
    case( RUN_LOG_LEVEL_TRACE )
      log_level = LOG_LEVEL_TRACE
    end select

    call log_set_level( log_level )

    write(log_scratch_space,'(A)')                              &
        'Runtime message logging severity set to log level: '// &
        convert_to_upper(key_from_run_log_level(run_log_level))
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    ! Initialise timer
    if ( subroutine_timers ) then
      call init_timer(timer_output_path)
      call timer(program_name)
    end if

    call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Create the mesh
    allocate( local_mesh_collection, &
              source = local_mesh_collection_type() )
    allocate( mesh_collection, &
              source=mesh_collection_type() )

    stencil_depth = 1

    call init_mesh( local_rank, total_ranks, stencil_depth, &
                    mesh, twod_mesh )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh, chi, panel_id )

    ! Set up IO
    files_init_ptr => init_io_dev_files
    call init_io( xios_context_id, communicator, &
                  chi, panel_id,                 &
                  populate_filelist=files_init_ptr )

  end subroutine initialise_infrastructure

  !> @brief Finalises infrastructure and constants used by the model
  !> @param[in] program_name The model run identifier
  subroutine finalise_infrastructure(program_name)

    implicit none

    character(*), intent(in) :: program_name

    ! Finalise XIOS context if we used it for diagnostic output or checkpointing
    call final_io()

    ! Finalise timer
    if ( subroutine_timers ) then
      call timer(program_name)
      call output_timer()
    end if

    ! Finalise aspects of the grid
    call final_mesh()
    call final_fem()

    ! Final logging before infrastructure is destroyed
    call log_event( program_name//' completed.', LOG_LEVEL_ALWAYS )

    ! Finalise namelist configurations
    call final_configuration()

    ! Finalise halo functionality
    call finalise_halo_comms()

    ! Finalise the logging system
    call finalise_logging()

  end subroutine finalise_infrastructure

end module io_dev_model_mod
