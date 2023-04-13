!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @page linear Linear model Program
!> This is a code that uses the LFRic infrastructure to build a model that
!> is the tangent linear/ perturbation forecast for gungho.

!> @brief Main program used to illustrate linear model functionality.

!> @details This top-level code simply calls initialise, run and finalise
!>          routines that are required to run the model.

program linear_model

  use linear_driver_mod, only : initialise, run, finalise

  implicit none

  character(*), parameter :: application_name = "linear_model"

  call initialise( application_name )

  call run( application_name )

  call finalise( application_name )

end program linear_model
