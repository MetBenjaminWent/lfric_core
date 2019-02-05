!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page Miniapp cell_locator program

!> @brief Finds the cell index that contains a target point

!> @details Calls init, run and finalise routines from a driver module

program cell_locator

  use cli_mod,             only : get_initial_filename
  use cell_locator_driver_mod, only : initialise, run, finalise

  implicit none

  character(:), allocatable :: filename

  call get_initial_filename( filename )
  call initialise( filename )
  deallocate( filename )

  call run()

  call finalise()

end program cell_locator
