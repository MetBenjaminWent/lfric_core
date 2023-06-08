!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Drives the execution of the io_dev miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module io_dev_driver_mod

  use base_mesh_config_mod,       only: prime_mesh_name
  use checksum_alg_mod,           only: checksum_alg
  use clock_mod,                  only: clock_type
  use constants_mod,              only: i_def, i_native, str_def, &
                                        PRECISION_REAL, r_def, r_second
  use convert_to_upper_mod,       only: convert_to_upper
  use driver_time_mod,            only: init_time, get_calendar
  use driver_mesh_mod,            only: init_mesh, final_mesh
  use driver_fem_mod,             only: init_fem, final_fem
  use driver_io_mod,              only: init_io, final_io, &
                                        filelist_populator, &
                                        get_io_context
  use extrusion_mod,              only: TWOD
  use field_mod,                  only: field_type
  use inventory_by_mesh_mod,      only: inventory_by_mesh_type
  use io_dev_config_mod,          only: multi_mesh, alt_mesh_name
  use io_config_mod,              only: write_diag, diagnostic_frequency, &
                                        subroutine_timers, timer_output_path
  use local_mesh_collection_mod,  only: local_mesh_collection, &
                                        local_mesh_collection_type
  use log_mod,                    only: log_event,          &
                                        log_scratch_space,  &
                                        LOG_LEVEL_ALWAYS,   &
                                        LOG_LEVEL_INFO
  use mesh_collection_mod,        only: mesh_collection, &
                                        mesh_collection_type
  use mesh_mod,                   only: mesh_type
  use model_clock_mod,            only: model_clock_type
  use mpi_mod,                    only: mpi_type
  use timer_mod,                  only: timer, output_timer, init_timer
  use io_dev_init_files_mod,      only: init_io_dev_files
  use io_dev_data_mod,            only: io_dev_data_type,          &
                                        create_model_data,         &
                                        initialise_model_data,     &
                                        update_model_data,         &
                                        output_model_data,         &
                                        finalise_model_data

  use io_context_mod, only: io_context_type
  use lfric_xios_context_mod, only: lfric_xios_context_type, advance

  implicit none

  private

  public initialise, run, finalise

  type(model_clock_type), allocatable :: model_clock


  contains

  !> @brief Sets up required state in preparation for run.
  subroutine initialise( model_data, mpi, program_name )

    implicit none

    class(io_dev_data_type), intent(inout) :: model_data
    class(mpi_type),         intent(inout) :: mpi
    character(*),            intent(in)    :: program_name


    type(field_type), pointer :: chi(:) => null()
    type(field_type), pointer :: panel_id => null()
    type(mesh_type),  pointer :: mesh => null()
    type(mesh_type),  pointer :: twod_mesh => null()
    type(mesh_type),  pointer :: alt_mesh => null()

    type(inventory_by_mesh_type)    :: chi_inventory
    type(inventory_by_mesh_type)    :: panel_id_inventory
    character(str_def), allocatable :: base_mesh_names(:)
    character(str_def), allocatable :: alt_mesh_names(:)

    class(io_context_type), pointer :: io_context

    procedure(filelist_populator), pointer :: files_init_ptr => null()

    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    if ( subroutine_timers ) then
      call init_timer(timer_output_path)
      call timer(program_name)
    end if

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------
    call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Initialise model clock and calendar
    call init_time( model_clock )

    ! Create the meshes used to test multi-mesh output
    if (multi_mesh) then
      allocate(base_mesh_names(2))
      base_mesh_names(2) = alt_mesh_name
    else
      allocate(base_mesh_names(1))
    end if

    base_mesh_names(1) = prime_mesh_name

    call init_mesh( mpi%get_comm_rank(), &
                    mpi%get_comm_size(), &
                    base_mesh_names )

    ! Create FEM specifics (function spaces and chi fields)
    call init_fem( mesh_collection, chi_inventory, panel_id_inventory )

    ! Create IO and instantiate the fields stored in model_data
    files_init_ptr => init_io_dev_files

    mesh => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)
    call chi_inventory%get_field_array(mesh, chi)
    call panel_id_inventory%get_field(mesh, panel_id)

    if (multi_mesh) then
      alt_mesh => mesh_collection%get_mesh(alt_mesh_name)
      allocate(alt_mesh_names(1))
      alt_mesh_names(1) = alt_mesh_name
      call create_model_data( model_data, chi, panel_id, &
                              mesh, twod_mesh, alt_mesh )
      call init_io( program_name, mpi%get_comm(),              &
                    chi_inventory, panel_id_inventory,         &
                    model_clock, get_calendar(),               &
                    populate_filelist = files_init_ptr,        &
                    model_data = model_data,                   &
                    alt_mesh_names = alt_mesh_names )

    else
      call create_model_data( model_data, chi, panel_id, &
                              mesh, twod_mesh )
      call init_io( program_name, mpi%get_comm(),              &
                    chi_inventory, panel_id_inventory,         &
                    model_clock, get_calendar(),               &
                    populate_filelist = files_init_ptr,        &
                    model_data = model_data )
    end if

    ! Initialise the fields stored in the model_data
    call initialise_model_data( model_data, model_clock, chi, panel_id )

    ! Write initial output
    io_context => get_io_context()
    if (model_clock%is_initialisation()) then
      select type (io_context)
      type is (lfric_xios_context_type)
          call advance(io_context, model_clock)
      end select
    end if

    nullify(mesh, twod_mesh, chi, panel_id, files_init_ptr)
    deallocate(base_mesh_names)

  end subroutine initialise

  !> @brief Timesteps the model, calling the desired timestepping algorithm
  !>        based upon the configuration
  !> @param [in,out] model_data The structure that holds model state
  subroutine run( model_data, program_name )

    implicit none

    class(io_dev_data_type), intent(inout) :: model_data
    character(*),             intent(in) :: program_name

    write(log_scratch_space,'(A)') 'Running '//program_name//' ...'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    ! Model step
    do while( model_clock%tick() )

      ! Update fields
      call update_model_data( model_data, model_clock )

      ! Write out diagnostics
      if (write_diag) then
        if ( (mod( model_clock%get_step(), diagnostic_frequency ) == 0) ) then
          call log_event( program_name//': Writing output', LOG_LEVEL_INFO)
          call output_model_data( model_data, model_clock )
        end if
      end if

    end do

  end subroutine run

  !> @brief Tidies up after a model run.
  !> @param [in,out] model_data The structure that holds model state
  subroutine finalise( model_data, program_name )

    implicit none

    class(io_dev_data_type), intent(inout) :: model_data
    character(*),            intent(in) :: program_name

    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Finalise IO context
    call final_io()

    ! Destroy the fields stored in model_data
    call finalise_model_data( model_data )

    ! Finalise timer
    if ( subroutine_timers ) then
      call timer( program_name )
      call output_timer()
    end if

    ! Finalise aspects of the grid
    call final_mesh()
    call final_fem()

  end subroutine finalise

end module io_dev_driver_mod
