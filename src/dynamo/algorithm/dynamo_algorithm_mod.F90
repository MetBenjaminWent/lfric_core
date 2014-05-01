!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> A simple algorithm for testing the psy layer.
!>
module dynamo_algorithm_mod

  use lfric
  use log_mod,    only: log_event, log_scratch_space, LOG_LEVEL_INFO
  use psy,        only: invoke_rhs_v3, invoke_v3_solver_kernel

  implicit none

  private
  public :: dynamo_algorithm

contains

  !> A simple algorithm which calls two kernels.
  !>
  subroutine dynamo_algorithm( pressure_density, rhs )

    implicit none

    type( field_type ), intent( in ) :: pressure_density
    type( field_type ), intent( in ) :: rhs

    !Construct PSy layer given a list of kernels. This is the line the code
    !generator may parse and do its stuff.

    call log_event( "Dynamo: calling 1st kernel", LOG_LEVEL_INFO )
    !PSY call invoke ( v3_rhs_kernel_type(rhs) )
    call invoke_rhs_v3( rhs )

    call log_event( "Dynamo:calling 2nd kernel", LOG_LEVEL_INFO )
    !PSY call invoke ( v3_solver_kernel_type(pressure_density,rhs) )
    call invoke_v3_solver_kernel( pressure_density, rhs )

  end subroutine dynamo_algorithm

end module dynamo_algorithm_mod
