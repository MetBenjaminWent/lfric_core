!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Provides access to the members of the psy class.

!> @details Accessor functions for the psy class are defined in this module.

!> @param invoke_RHS_V3              Invoke the RHS for a v3 field
!> @param invoke_v3_solver_kernel    Invoke the solver for a v3 field kernel

module psy

  use field_mod, only : field_data_type, field_type, field_data_from_proxy
  use lfric

  implicit none

contains

  subroutine invoke_rhs_v3( right_hand_side_proxy )

    use v3_rhs_kernel_mod, only : rhs_v3_code

    implicit none

    type( field_type ), intent( in ) :: right_hand_side_proxy

    type( field_data_type), pointer :: right_hand_side => null( )
    integer :: cell
    integer, pointer :: map(:)
    integer :: nlayers
    integer :: ndf
    real(kind=dp), pointer  :: v3_basis(:,:,:,:,:)

    right_hand_side => field_data_from_proxy( right_hand_side_proxy )

    ! Unpack data
    nlayers = right_hand_side%get_nlayers( )
    ndf = right_hand_side%vspace%get_ndf( )

    call right_hand_side%vspace%get_basis( v3_basis )
    do cell = 1, right_hand_side%get_ncell( )
       call right_hand_side%vspace%get_cell_dofmap( cell, map )
       call rhs_v3_code( nlayers, &
                         ndf, &
                         map, &
                         v3_basis, &
                         right_hand_side%data, &
                         right_hand_side%gaussian_quadrature )
    end do

  end subroutine invoke_rhs_v3

  subroutine invoke_v3_solver_kernel( pdfield, rhs )

    use v3_solver_kernel_mod, only : solver_v3_code

    type( field_type ), intent( in ) :: pdfield
    type( field_type ), intent( in ) :: rhs

    integer                 :: cell
    integer, pointer        :: map(:)
    integer                 :: nlayers
    integer                 :: ndf
    real(kind=dp), pointer  :: v3_basis(:,:,:,:,:)

    type( field_data_type ), pointer :: pd_data  => null( )
    type( field_data_type ), pointer :: rhs_data => null( )

    pd_data  => field_data_from_proxy( pdfield )
    rhs_data => field_data_from_proxy( rhs )

    nlayers = pd_data%get_nlayers( )
    ndf     = pd_data%vspace%get_ndf( )
    call pd_data%vspace%get_basis( v3_basis )

    do cell = 1, pd_data%get_ncell()
       call pd_data%vspace%get_cell_dofmap( cell, map )
       call solver_v3_code( nlayers, &
                            ndf, &
                            map, &
                            v3_basis, &
                            pd_data%data, &
                            rhs_data%data, &
                            pd_data%gaussian_quadrature )
    end do

  end subroutine invoke_v3_solver_kernel

end module psy
