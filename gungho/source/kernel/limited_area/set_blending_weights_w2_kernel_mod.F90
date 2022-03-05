!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Populates the w2 weights field with the blending weights obtained from the
!>        configuration.
!> @details Uses the onion layer values as array index for blending_weights array
!>          to populate a field of blending weights.
!>          The onion layers are on the discontinuous W3 space.  For the W2 dofs,
!>          average values are taken using the weights from either side of the face.
module set_blending_weights_w2_kernel_mod

  use argument_mod,              only : arg_type,             &
                                        GH_SCALAR, GH_FIELD,  &
                                        GH_READ, GH_REAL,     &
                                        GH_INC,               &
                                        GH_INTEGER, GH_BASIS, &
                                        CELL_COLUMN
  use fs_continuity_mod,         only : W3, W2
  use constants_mod,             only : r_def, i_def
  use kernel_mod,                only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------

  type, public, extends(kernel_type) :: set_blending_weights_w2_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/              &
         arg_type(GH_FIELD,   GH_REAL, GH_INC, W2),  & ! weights_field
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W3), & ! onion_layers
         arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)   & ! depth
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: set_blending_weights_w2_code
  end type

  !-------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------
  public :: set_blending_weights_w2_code

contains

!> @param[in] nlayers      Number of layers
!> @param[in,out] weights_field The output field
!> @param[in] onion_layers The input field
!> @param[in] depth     Depth of the weight array
!> @param[in] ndf_out   Number of degrees of freedom for weights_field
!> @param[in] undf_out  Total number of degrees of freedom for weights_field
!> @param[in] map_out   Dofmap for the cell at the base of the column for weights_field
!> @param[in] ndf_in    Number of degrees of freedom for onion_layers
!> @param[in] undf_in   Total number of degrees of freedom for onion_layers
!> @param[in] map_in    Dofmap for the cell at the base of the column for onion_layers
subroutine set_blending_weights_w2_code( nlayers,   &
                                         weights_field, &
                                         onion_layers,  &
                                         depth,     &
                                         ndf_out,   &
                                         undf_out,  &
                                         map_out,   &
                                         ndf_in,    &
                                         undf_in,   &
                                         map_in)

  use boundaries_config_mod,        only : blending_weights

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in) :: nlayers
  integer(kind=i_def),                     intent(in) :: ndf_out, undf_out
  integer(kind=i_def),                     intent(in) :: ndf_in, undf_in
  real(kind=r_def), dimension(undf_out),   intent(inout) :: weights_field
  real(kind=r_def), dimension(undf_in),    intent(in) :: onion_layers
  integer(kind=i_def),                     intent(in) :: depth
  integer(kind=i_def), dimension(ndf_out), intent(in) :: map_out
  integer(kind=i_def), dimension(ndf_in),  intent(in) :: map_in

  ! Internal variables
  integer(kind=i_def)  :: k, df
  integer(kind=i_def)  :: index

  ! W2 weights should not go right up to the inner region
  if (onion_layers(map_in(1)) > 1.0_r_def)then
    index = depth - INT(onion_layers(map_in(1))) + 1
    do k=0,nlayers-1
      do df=1,ndf_out
        ! W2 weights are averages of neighouring W3 weights
        weights_field(map_out(df)+k) = weights_field(map_out(df)+k) &
                      + 0.5_r_def*blending_weights(index)
      end do
    end do
  end if

end subroutine set_blending_weights_w2_code

end module set_blending_weights_w2_kernel_mod
