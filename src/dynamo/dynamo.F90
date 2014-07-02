!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @mainpage Dynamo
!> Illustration of the PSyKAl (Parallel-system/Kernel/Algorithm) architecture
!> for Gung Ho. Whilst the computational and optimisation infrastructure is
!> being developed, the science code is being developed using 
!> a hand-rolled PSy layer, PSy-lite. A PSyKAl-lite needs a dynamo!
!> Eventually, PSyKAl-lite will be replaced with the real PSy and Dynamo
!> will be the implementation of the Gung Ho dynamical core.

!> @brief Main program used to illustrate dynamo functionality.

!> @details Calls Creates function spaces, then fields on those 
!> function spaces, before passing the fields to the algorithm layer

program dynamo

  use dynamo_algorithm_mod,    only : dynamo_algorithm
  use field_mod,               only : field_type
  use function_space_mod,      only : function_space_type, V1, V2, V3
  use log_mod,                 only : log_event, LOG_LEVEL_INFO
  use set_up_mod,              only : set_up
  use gaussian_quadrature_mod, only : gaussian_quadrature_type
  use mesh_mod,                only : num_layers

  implicit none

  type( function_space_type )      :: function_space
  type( field_type )               :: pressure_density, rhs
  type( field_type )               :: flux_velocity, rhs_v2
  type( field_type )               :: rhs_v1, circulation
  type( gaussian_quadrature_type ) :: gq

  call log_event( 'Dynamo running...', LOG_LEVEL_INFO )

  call set_up( )


  pressure_density = field_type( function_space%get_instance(V3),          &
                                 gq%get_instance(),                        &
                                 num_layers = num_layers )

  rhs = field_type( function_space%get_instance(V3),                       &
                    gq%get_instance(),                                     &
                    num_layers = num_layers )

  flux_velocity = field_type( function_space%get_instance(V2),             &
                         gq%get_instance(),                                &
                         num_layers = num_layers )
  
  rhs_v2 = field_type( function_space%get_instance(V2),                    &
                       gq%get_instance(),                                  &
                       num_layers = num_layers )

  circulation = field_type( function_space%get_instance(V1),               &
                            gq%get_instance(),                             &
                            num_layers = num_layers )

  rhs_v1 = field_type( function_space%get_instance(V1),                    &
                       gq%get_instance(),                                  &
                       num_layers = num_layers )

  call dynamo_algorithm( pressure_density, rhs,                            &
                         flux_velocity, rhs_v2,                            &
                         circulation, rhs_v1 )

  call rhs%print_field( "RHS field..." )
  call pressure_density%print_field( "LHS field..." )

  call rhs_v2%print_field( "RHS field...v2" )
  call flux_velocity%print_field( "flux_velocity ..." )

  call rhs_v1%print_field("RHS v1 ...")
  call circulation%print_field("circulation ...")

  call log_event( 'Dynamo completed', LOG_LEVEL_INFO )

end program dynamo
