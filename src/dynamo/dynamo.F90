!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @mainpage Dynamo
!> Illustration of the PsyKAl (Parallel-system/Kernel/Algorithm) architecture
!> for Gung Ho. Whilst the computational and optimisation infrastructure is
!> being developed, the science code is being developed using 
!> a hand-rolled Psy layer, Psy-lite. A PsyKAl-lite needs a dynamo!
!> Eventually, PsyKAl-lite will be replaced with the real Psy and Dynamo
!> will become Gung Ho.

!> @brief Main program used to illustrate dynamo functionality.

!> @details Creates the function spaces and calls <code>set_up</code> to
!> populate them (either read in or compute) then individual calls to the
!> psy-layer with kernels as if the code has been pre-processed by Psyclone.

program dynamo

  use dynamo_algorithm_mod, only : dynamo_algorithm
  use field_mod,            only : field_data_type
  use lfric
  use log_mod,              only : log_event, LOG_LEVEL_INFO
  use set_up_mod,           only : set_up

  implicit none

  type( function_space_type )      :: v3_function_space, v2_function_space, &
                                      v1_function_space, v0_function_space
  type( field_data_type )          :: pressure_density, rhs
  type( gaussian_quadrature_type ) :: gq
  integer                          :: num_layers

  call log_event( 'Dynamo running...', LOG_LEVEL_INFO )

  call set_up( v0_function_space, v1_function_space, v2_function_space, &
               v3_function_space, num_layers )

  gq = gaussian_quadrature_type( )

  pressure_density = field_data_type( vector_space = v3_function_space, &
                                      gq = gq,                          &
                                      num_layers = num_layers)

  rhs = field_data_type( vector_space = v3_function_space, &
                         gq = gq,                          &
                         num_layers = num_layers )

  call dynamo_algorithm( pressure_density%new_proxy( ), rhs%new_proxy( ) )

  call print_field( 'RHS field...', rhs )
  call print_field( 'LHS field...', pressure_density )

  call log_event( 'Dynamo completed', LOG_LEVEL_INFO )

contains

  !> Send a field to the log.
  !>
  subroutine print_field( title, field )

    use lfric
    use log_mod, only : log_event, log_scratch_space, LOG_LEVEL_DEBUG

    implicit none

    character( * ),          intent( in ) :: title
    type( field_data_type ), intent( in ) :: field

    integer                   :: cell
    integer                   :: layer
    integer                   :: df
    integer,          pointer :: map( : )

    call log_event( title, LOG_LEVEL_DEBUG )

    do cell=1,field%vspace%get_ncell()
      call field%vspace%get_cell_dofmap(cell,map)
      do df=1,field%vspace%get_ndf()
        do layer=0,field%get_nlayers()-1
          write( log_scratch_space, '( I4, I4, I4, F8.2 )' ) &
              cell, df, layer+1, field%data( map( df ) + layer )
          call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
        end do
      end do
    end do

  end subroutine print_field

end program dynamo
