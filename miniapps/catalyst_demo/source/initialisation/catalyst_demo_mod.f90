!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> Catalyst demonstration program support functions.
!>
!> Originally these were "block" constructs within the program but neither
!> GNU or Intel Fortran where properly able to cope with that.
!>
module catalyst_demo_mod

  use configuration_mod, only : read_configuration,   &
                                ensure_configuration, &
                                final_configuration

  use log_mod, only : log_event,         &
                      log_scratch_space, &
                      LOG_LEVEL_ERROR,   &
                      LOG_LEVEL_TRACE,   &
                      LOG_LEVEL_DEBUG

  implicit none

  private
  public :: load_configuration, final_configuration

contains

  !> Loads run-time configuration and ensures everything is ship-shape.
  !>
  subroutine load_configuration( filename )


    implicit none

    character(*), intent(in) :: filename

    character(*), parameter :: &
                            required_configuration(10) = ['base_mesh             ', &
                                                          'planet                ', &
                                                          'extrusion             ', &
                                                          'initial_temperature   ', &
                                                          'initial_wind          ', &
                                                          'visualisation         ', &
                                                          'timestepping          ', &
                                                          'multigrid             ', &
                                                          'gravity_wave_constants', &
                                                          'domain_size           ']

    logical              :: okay
    logical, allocatable :: success_map(:)
    integer              :: i

    allocate( success_map(size(required_configuration)) )

    call read_configuration( filename )

    okay = ensure_configuration( required_configuration, success_map )
    if (.not. okay) then
      write( log_scratch_space, '(A)' ) &
                             'The following required namelists were not loaded:'
      do i = 1,size(required_configuration)
        if (.not. success_map(i)) &
          log_scratch_space = trim(log_scratch_space) // ' ' &
                              // required_configuration(i)
      end do
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    deallocate( success_map )

  end subroutine load_configuration

end module catalyst_demo_mod
