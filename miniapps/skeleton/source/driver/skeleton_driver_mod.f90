!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the skeleton miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module skeleton_driver_mod

  use base_mesh_config_mod,       only : prime_mesh_name
  use checksum_alg_mod,           only : checksum_alg
  use constants_mod,              only : i_def, i_native, str_def, &
                                         PRECISION_REAL, r_def, r_second
  use convert_to_upper_mod,       only : convert_to_upper
  use driver_time_mod,            only : init_time, get_calendar
  use driver_mesh_mod,            only : init_mesh, final_mesh
  use driver_fem_mod,             only : init_fem, final_fem
  use driver_io_mod,              only : init_io, final_io
  use field_mod,                  only : field_type
  use init_skeleton_mod,          only : init_skeleton
  use inventory_by_mesh_mod,      only : inventory_by_mesh_type
  use io_config_mod,              only : write_diag
  use log_mod,                    only : log_event, log_scratch_space, &
                                         LOG_LEVEL_ALWAYS, LOG_LEVEL_INFO
  use mesh_mod,                   only : mesh_type
  use mesh_collection_mod,        only : mesh_collection
  use model_clock_mod,            only : model_clock_type
  use mpi_mod,                    only : mpi_type
  use skeleton_alg_mod,           only : skeleton_alg

  implicit none

  private
  public initialise, run, finalise

  type(model_clock_type), allocatable :: model_clock

  ! Prognostic fields
  type( field_type ) :: field_1

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !>
  subroutine initialise( mpi, program_name )

    implicit none

    class(mpi_type), intent(inout) :: mpi
    character(*),    intent(in)    :: program_name

    real(r_def) :: dt_model

    ! Coordinate field
    type(field_type),             pointer :: chi(:) => null()
    type(field_type),             pointer :: panel_id => null()
    type(mesh_type),              pointer :: mesh => null()
    type(inventory_by_mesh_type)          :: chi_inventory
    type(inventory_by_mesh_type)          :: panel_id_inventory
    character(str_def),       allocatable :: base_mesh_names(:)


    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------

    ! Initialise clock
    call init_time( model_clock )
    dt_model = real(model_clock%get_seconds_per_step(), r_def)

    ! Create the mesh
    allocate(base_mesh_names(1))
    base_mesh_names(1) = prime_mesh_name
    call init_mesh( mpi%get_comm_rank(), mpi%get_comm_size(), base_mesh_names )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh_collection, chi_inventory, panel_id_inventory )

    ! Initialise I/O context
    call init_io( program_name, mpi%get_comm(), chi_inventory, &
                  panel_id_inventory, model_clock, get_calendar() )

    ! Create and initialise prognostic fields
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    call chi_inventory%get_field_array(mesh, chi)
    call panel_id_inventory%get_field(mesh, panel_id)
    call init_skeleton( mesh, chi, panel_id, dt_model, field_1 )

    nullify(mesh, chi, panel_id)
    deallocate(base_mesh_names)

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time steps.
  !>
  subroutine run( program_name )

    implicit none

    character(*), intent(in) :: program_name

    logical :: running

    running = model_clock%tick()

    ! Call an algorithm
    call skeleton_alg(field_1)

    ! Write out output file
    call log_event(program_name//": Writing diagnostic output", LOG_LEVEL_INFO)

    if (write_diag ) then
      ! Calculation and output of diagnostics
      call field_1%write_field('skeleton_field')
    end if

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !>
  subroutine finalise( program_name )

    implicit none

    character(*), intent(in) :: program_name

    !--------------------------------------------------------------------------
    ! Model finalise
    !--------------------------------------------------------------------------

    ! Write checksums to file
    call checksum_alg(program_name, field_1, 'skeleton_field_1')

    call log_event( program_name//': Miniapp completed', LOG_LEVEL_INFO )

    !-------------------------------------------------------------------------
    ! Driver layer finalise
    !-------------------------------------------------------------------------

    ! Finalise IO
    call final_io()

    call final_fem()

    call final_mesh()

  end subroutine finalise

end module skeleton_driver_mod
