!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing the JEDI pseudo model emulator class.
!>
!> @details This module holds a JEDI pseudo model emulator where the
!>          time-stepping uses the read_file state method. An included forecast
!>          application uses the model init, step and final to run a forecast.
!>
module jedi_pseudo_model_mod

  use constants_mod,                 only : i_def, str_def
  use jedi_lfric_datetime_mod,       only : jedi_datetime_type
  use jedi_lfric_duration_mod,       only : jedi_duration_type
  use jedi_state_mod,                only : jedi_state_type
  use log_mod,                       only : log_event,          &
                                            log_scratch_space,  &
                                            LOG_LEVEL_ERROR

  implicit none

  private

type, public :: jedi_pseudo_model_type
  private

  !> A date_time list to read
  type( jedi_datetime_type ), allocatable  :: datetime_states(:)
  !> Index for the current state in date_time_states
  integer( kind=i_def )                    :: current_state
  !> The number of states
  integer( kind=i_def )                    :: n_states
  !> The file prefix for reading
  character( len=str_def )                 :: file_prefix

contains

  !> Model initialiser.
  procedure, public  :: initialise

  !> Methods
  procedure, private :: model_init
  procedure, private :: model_step
  procedure, private :: model_final

  !> Run a forecast
  procedure, public  :: forecast

  !> Finalizer
  final              :: jedi_pseudo_model_destructor

end type jedi_pseudo_model_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> @brief    Initialiser for jedi_pseudo_model_type
!>
!> @param [in] config Configuration used to setup the model class
subroutine initialise( self, config )

  use jedi_pseudo_model_config_mod, only : jedi_pseudo_model_config_type

  implicit none

  class( jedi_pseudo_model_type ),       intent(inout)  :: self
  type( jedi_pseudo_model_config_type ), intent(in)     :: config

  ! Setup the pseudo model
  self%current_state = 1_i_def
  self%n_states = size( config%datetime_states, dim=1 )
  allocate( self%datetime_states(self%n_states) )
  self%datetime_states = config%datetime_states
  self%file_prefix = config%read_file_prefix

end subroutine initialise

!> @brief    Initialise the model
!>
!> @param [inout] state State object to be used in the model initialise
subroutine model_init( self, state )

  implicit none

  class( jedi_pseudo_model_type ), intent(in)    :: self
  type( jedi_state_type ),         intent(inout) :: state

end subroutine model_init

!> @brief    Step the model
!>
!> @param [inout] state State object to propagated
subroutine model_step( self, state )

  implicit none

  class( jedi_pseudo_model_type ), intent(inout) :: self
  type( jedi_state_type ),         intent(inout) :: state

  ! Local
  type( jedi_datetime_type ) :: datetime

  if ( self%current_state > self%n_states ) then
    write ( log_scratch_space, '(A)' ) "self%current_state>self%n_states."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  endif

  datetime = self%datetime_states( self%current_state )
  call state%read_file( datetime, self%file_prefix )

  ! Iterate the current state
  self%current_state = self%current_state + 1_i_def

end subroutine model_step

!> @brief    Finalise the model
!>
!> @param [inout] state State object to be used in the model finalise
subroutine model_final( self, state )

  implicit none

  class( jedi_pseudo_model_type ), intent(in)    :: self
  type( jedi_state_type ),         intent(inout) :: state

end subroutine model_final

!> @brief    Finalize the jedi_pseudo_model_type
!>
subroutine jedi_pseudo_model_destructor( self )

  implicit none

  type( jedi_pseudo_model_type ), intent(inout) :: self

  if ( allocated(self%datetime_states) ) deallocate(self%datetime_states)

end subroutine jedi_pseudo_model_destructor

!------------------------------------------------------------------------------
! OOPS defined forecast method
!------------------------------------------------------------------------------

!> @brief    Run a forecast using the model init, step and final
!>
!> @param [inout] state           The state object to propagate
!> @param [in] datetime_duration  The duration of the forecast
subroutine forecast( self, state, datetime_duration )

  implicit none

  class( jedi_pseudo_model_type ), intent(inout) :: self
  type( jedi_state_type ),         intent(inout) :: state
  type( jedi_duration_type ),      intent(in)    :: datetime_duration

  ! Local
  type( jedi_datetime_type ) :: datetime_end

  ! End time
  datetime_end = state%datetime + datetime_duration

  call self%model_init( state )

  ! initialize the post processor
  ! post.initialize(state, ...);
  ! In a H(X) the following call would be used to do spatial interpolations
  ! on the Atlas/dummy field data
  ! post%process(state)

  ! loop until date_time_end
  do while ( datetime_end%is_ahead( state%datetime ) )
    call self%model_step( state )
    ! post%process(state)
  end do

  ! post.finalize(state);
  call self%model_final( state )

end subroutine forecast

end module jedi_pseudo_model_mod
