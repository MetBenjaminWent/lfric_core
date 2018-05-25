!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

! A very simply program which just logs an error.
!
program log_mod_error_test

  use ESMF,            only : ESMF_Initialize, ESMF_Finalize
  use iso_fortran_env, only : error_unit
  use mpi_mod,         only : store_comm, get_comm_size
  use log_mod,         only : log_event, log_set_parallel_logging, &
                              LOG_LEVEL_ERROR
  use mpi,             only : MPI_COMM_WORLD

  integer :: condition

  ! ESMF is needed even for a serial run as it is used to determine that it
  ! *is* a serial run.
  !
  ! Since the logging module is being tested it is important not to use it to
  ! handle a failure in ESMF. Instead we do it the old fashioned way.
  call ESMF_Initialize( rc=condition )
  if (condition /=0 )then
    write( error_unit, '("Failed to initialise ESMF: ", I0)') condition
    stop 1
  end if
  call store_comm(MPI_COMM_WORLD)
  if(get_comm_size() > 1) call log_set_parallel_logging(.true.)

  ! Everything else is here purely to support the testing of this line:
  !
  call log_event( 'An error was logged.', LOG_LEVEL_ERROR )

  ! Since the logging module is being tested it is important not to use it to
  ! handle a failure in ESMF. Instead we do it the old fashioned way.
  call ESMF_Finalize( rc=condition )
  if (condition /=0 )then
    write( error_unit, '("Failed to finalise ESMF: ", I0)') condition
    stop 2
  end if

end program log_mod_error_test
