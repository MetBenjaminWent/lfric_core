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

  use constants_mod,                  only : i_def
  use cli_mod,                        only : get_initial_filename
  use skeleton_mod,                   only : load_configuration
  use init_mesh_mod,                  only : init_mesh
  use init_fem_mod,                   only : init_fem
  use init_skeleton_mod,              only : init_skeleton
  use ESMF
  use global_mesh_collection_mod,     only : global_mesh_collection, &
                                             global_mesh_collection_type
  use field_mod,                      only : field_type
  use skeleton_alg_mod,               only : skeleton_alg
  use derived_config_mod,             only : set_derived_config
  use log_mod,                        only : log_event,         &
                                             log_set_level,     &
                                             log_scratch_space, &
                                             LOG_LEVEL_ERROR,   &
                                             LOG_LEVEL_INFO
  use output_config_mod,              only : write_nodal_output, &
                                             write_xios_output
  use restart_config_mod,             only : restart_filename => filename
  use restart_control_mod,            only : restart_type
  use io_mod,                         only : output_nodal, &
                                             output_xios_nodal, &
                                             xios_domain_init
  use checksum_alg_mod,               only : checksum_alg

  use xios
  use mpi
  use mod_wait

  implicit none

  private
  public initialise, run, finalise

  type(restart_type) :: restart

  ! Prognostic fields
  type( field_type ) :: field_1

  ! Coordinate field
  type(field_type), target, dimension(3) :: chi

  integer(i_def) :: mesh_id, twod_mesh_id

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !>
  subroutine initialise( filename )

    implicit none

    character(:), intent(in), allocatable :: filename

    character(len=*), parameter   :: xios_id   = "lfric_client"
    character(len=*), parameter   :: xios_ctx  = "skeleton"

    type(ESMF_VM) :: vm

    integer(i_def) :: rc
    integer(i_def) :: total_ranks, local_rank
    integer(i_def) :: petCount, localPET, ierr
    integer(i_def) :: comm = -999
    integer(i_def) :: dtime

    ! Initialise MPI

    call mpi_init(ierr)

    ! Initialise XIOS and get back the split mpi communicator
    call init_wait()
    call xios_initialize(xios_id, return_comm = comm)

    ! Initialise ESMF using mpi communicator initialised by XIOS
    ! and get the rank information from the virtual machine
    call ESMF_Initialize(vm=vm, &
                        defaultlogfilename="skeleton.Log", &
                        logkindflag=ESMF_LOGKIND_MULTI, &
                        mpiCommunicator=comm, &
                        rc=rc)

    if (rc /= ESMF_SUCCESS) call log_event( 'Failed to initialise ESMF.', &
                                            LOG_LEVEL_ERROR )

    call ESMF_VMGet(vm, localPet=localPET, petCount=petCount, rc=rc)
    if (rc /= ESMF_SUCCESS) &
      call log_event( 'Failed to get the ESMF virtual machine.', &
                      LOG_LEVEL_ERROR )

    total_ranks = petCount
    local_rank  = localPET

    ! Currently log_event can only use ESMF so it cannot be used before ESMF
    ! is initialised.
    call log_event( 'skeleton: Running miniapp ...', LOG_LEVEL_INFO )

    call load_configuration( filename )
    call set_derived_config( .true. )

    restart = restart_type( restart_filename, local_rank, total_ranks )

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------

    allocate( global_mesh_collection, &
              source = global_mesh_collection_type() )

    ! Create the mesh
    call init_mesh(local_rank, total_ranks, mesh_id, twod_mesh_id)

    ! Full global meshes no longer required, so reclaim
    ! the memory from global_mesh_collection
    write(log_scratch_space,'(A)') &
        "Purging global mesh collection."
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    deallocate(global_mesh_collection)

    ! Create FEM specifics (function spaces and chi field)
    call init_fem(mesh_id, chi)

    !-------------------------------------------------------------------------
    ! IO init
    !-------------------------------------------------------------------------

    ! If using XIOS for diagnostic output or checkpointing, then set up
    ! XIOS domain and context

    if ( (write_xios_output) .or. (restart%use_xios()) ) then

      dtime = 1

      call xios_domain_init( xios_ctx,   &
                             comm,       &
                             dtime,      &
                             restart,    &
                             mesh_id,    &
                             chi,        &
                             vm,         &
                             local_rank, &
                             total_ranks )

      ! Make sure XIOS calendar is set to timestep 1 as it starts there
      ! not timestep 0.
      call xios_update_calendar(1)

    end if


    ! Create and initialise prognostic fields
    call init_skeleton(mesh_id, chi, field_1)

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time steps.
  !>
  subroutine run()

    implicit none

    ! Call an algorithm
    call skeleton_alg(field_1)
  

    ! Write out output file
    call log_event("skeleton: Writing diagnostic output", LOG_LEVEL_INFO)

    ! Original nodal output
    if ( write_nodal_output)  then
      call output_nodal('skeleton_field', 0, field_1, mesh_id)
    end if

    ! XIOS output
    if (write_xios_output) then
      call output_xios_nodal("skeleton_field", field_1, mesh_id)
    end if

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !>
  subroutine finalise()

    implicit none

    integer(i_def) :: rc
    integer(i_def) :: ierr

    !-----------------------------------------------------------------------------
    ! Model finalise
    !-----------------------------------------------------------------------------

    ! Write checksums to file
    call checksum_alg('skeleton', field_1, 'skeleton_field_1')

    call log_event( 'skeleton: Miniapp completed', LOG_LEVEL_INFO )

    !-------------------------------------------------------------------------
    ! Driver layer finalise
    !-------------------------------------------------------------------------

   ! Finalise XIOS context if we used it for diagnostic output or checkpointing
    if ( (write_xios_output) .or. (restart%use_xios()) ) then
      call xios_context_finalize()
    end if

    ! Finalise XIOS
    call xios_finalize()

    ! Finalise ESMF
    call ESMF_Finalize(endflag=ESMF_END_KEEPMPI,rc=rc)

    ! Finalise mpi
    call mpi_finalize(ierr)

  end subroutine finalise

end module skeleton_driver_mod
