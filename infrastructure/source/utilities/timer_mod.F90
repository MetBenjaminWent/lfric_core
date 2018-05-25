!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you should have
! received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief A Simple subroutine timer based upon calls to cpu time
module timer_mod   
   use constants_mod,   only: i_def, r_def, str_def

   implicit none

   private

   integer(i_def),        parameter            :: num_subs       = 100
   integer(i_def)                              :: num_tim_in_use = 0
   character(len=str_def)                      :: routine_name(num_subs)
   real(r_def),           dimension(num_subs)  :: tot_time
   real(r_def),           dimension(num_subs)  :: prev_time
   integer,               dimension(num_subs)  :: num_calls
   logical,               dimension(num_subs)  :: start_stop

   public  :: timer
   public  :: output_timer
   private :: convert_to_lower
   
   ! Routines only needed for unit test
   public :: get_routine_name
   public :: get_routine_total_calls

contains
!=============================================================================!
!> @brief start/stop recording the runtime of a give section
!> @param[in] cname Name of the timing section to start/stop
   subroutine timer(cname)

     use log_mod,    only: log_event,         &
                           LOG_LEVEL_ERROR

     implicit none

     character(len=*),          intent(in)  :: cname
     character(len=str_def)                 :: lowname
     integer                                :: k
     real(r_def)                            :: t

     ! check if this is the first call to timer
     if( num_tim_in_use == 0 ) start_stop(:) = .false.

     ! convert cname to lower case
     lowname = convert_to_lower(cname)

     do k = 1, num_tim_in_use
        if( lowname == routine_name(k) ) exit
      end  do

     if( k > num_tim_in_use ) then
     ! subroutine not in list so initialise
       num_tim_in_use = k
       if( num_tim_in_use > num_subs ) then
         call log_event( "Run out of timers", LOG_LEVEL_ERROR )
       end if
       routine_name(k) = lowname
       call cpu_time(prev_time(k))
       start_stop(k) = .true.
       tot_time(k)   = 0.0_r_def
       num_calls(k)  = 1
     else

       ! Found routine check to see if its the start or end of
       ! a timing section
       if( start_stop(k) ) then
         call cpu_time(t)
         tot_time(k)   = tot_time(k) + t - prev_time(k)
         start_stop(k) = .false.
       else
         start_stop(k) = .true.               
         call cpu_time(prev_time(k))
         num_calls(k)  = num_calls(k) + 1
       endif
     endif

   end subroutine timer

!=============================================================================!
   !> @brief write out timer information to file
   subroutine output_timer()
     use mpi_mod,    only: get_comm_rank
     use log_mod,    only: log_event,         &
                           LOG_LEVEL_ERROR,   &
                           LOG_LEVEL_INFO,    &
                           log_scratch_space
     use scalar_mod, only: scalar_type
     implicit none
     integer(i_def)    :: k
     real(r_def)       :: pc_time
     type(scalar_type) :: time
     integer(i_def)    :: stat

     ! check all timers are closed
     do k = 1, num_tim_in_use
       if( start_stop(k) ) then
         write( log_scratch_space, '(A,A,A)') &
                    'Timer for routine ',trim(routine_name(k)), &
                    ' not closed. Timing information will be incorrect'
         call log_event( log_scratch_space, LOG_LEVEL_INFO )
       end if
     end do

     do k = 1, num_tim_in_use
       time = scalar_type(tot_time(k))
       tot_time(k) = time%get_sum()
     end do
     if ( get_comm_rank() == 0 ) then
       open( 9, file='timer.txt', status="replace", iostat=stat)
       if (stat /= 0) then
         call log_event( "Unable to open timer file", LOG_LEVEL_ERROR )
       end if

       ! Write out timer information in wiki formatted table
       write(9,'(A)')  &
       '||=  Routine =||= total time(s) =||= No. calls =||= %time =||= time per call(s) =||'
       do k = 1, num_tim_in_use
         pc_time = tot_time(k)/tot_time(1)*100.0_r_def
         if ( pc_time > 1.0_r_def ) then
           write(9,('(A,A50,A,f14.2,A,I12,A,f14.2,A,f14.2,A)')) &
                '||', trim(routine_name(k)),'||', &
                      tot_time(k),'||',          &
                      num_calls(k),'||',&
                      pc_time,'||', &
                      tot_time(k)/REAL(num_calls(k)),'||'
         end if
       end do
       close(9)
     end if
   end subroutine output_timer

!=============================================================================!
!> @brief Changes a string to lower case
!> @param[in] str Input string to convert
!> @result string Lower case string
  pure function convert_to_lower (str) Result (string)

    implicit none
    character(*), intent(in) :: str
    character(len(str))      :: string

    integer :: ic, i

    character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    string = str
    do i = 1, len_trim(str)
        ic = index(cap, str(i:i))
        if (ic > 0) string(i:i) = low(ic:ic)
    end do

  end function convert_to_lower
!=============================================================================!

!> @brief return the routine name for a given index
!> @param[in] idx index of routine
!> @result nme name of routine
  pure function get_routine_name(idx) result(nme)
    implicit none
    integer(i_def), intent(in) :: idx
    character(str_def)         :: nme
    nme = routine_name(idx)
  end function get_routine_name
!=============================================================================!

!> @brief return the total number of calls for a given index
!> @param[in] idx index of routine
!> @result n number of calls
  pure function get_routine_total_calls(idx) result(n)
    implicit none
    integer(i_def), intent(in) :: idx
    integer(i_def)             :: n
    n = num_calls(idx)
  end function get_routine_total_calls
!=============================================================================!

end module timer_mod
