!-------------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Convert a higher order field to a multi-data field
!> @details Temporary infrastructure required to convert a higher order
!>          field into a multi-data field for IO purposes. It will be retired
!>          when multi-dimensional fields are properly implemented.
!>          Also used without propoer psyclone support of multi-data fields
!>          see https://github.com/stfc/PSyclone/issues/868

module high_to_multi_kernel_mod

  use argument_mod,  only: arg_type, CELLS,           &
                           GH_FIELD, GH_INTEGER,      &
                           GH_READ, GH_WRITE,         &
                           ANY_DISCONTINUOUS_SPACE_1, &
                           ANY_DISCONTINUOUS_SPACE_2
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: high_to_multi_kernel_type
      private
      type(arg_type) :: meta_args(3) = (/                             &
          arg_type(GH_FIELD,   GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),  & ! multi-data field
          arg_type(GH_FIELD,   GH_READ,  ANY_DISCONTINUOUS_SPACE_2),  & ! high order field
          arg_type(GH_INTEGER, GH_READ               )                & ! ndata
          /)
      integer :: iterates_over = CELLS
  contains
      procedure, nopass :: high_to_multi_code
  end type

  public high_to_multi_code

contains

  !> @param[in]     nlayers       The number of layers
  !> @param[out]    multi_field   Multi-data field to write to
  !> @param[in]     high_field    Higher order field to write from
  !> @param[in]     ndata         Dimension of multi-data field
  !> @param[in]     ndf_multi     Number of DOFs per cell for multi-data field
  !> @param[in]     undf_multi    Number of total DOFs for multi-data field
  !> @param[in]     map_multi     Dofmap for cell for multi-data fields
  !> @param[in]     ndf_high      Number of DOFs per cell for high-order field
  !> @param[in]     undf_high     Number of total DOFs for high-order field
  !> @param[in]     map_high      Dofmap for cell for high-order fields
  subroutine high_to_multi_code(nlayers,                          &
                                multi_field,                      &
                                high_field,                       &
                                ndata,                            &
                                ndf_multi, undf_multi, map_multi, &
                                ndf_high, undf_high, map_high)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, ndata
    integer(kind=i_def), intent(in) :: ndf_multi, undf_multi
    integer(kind=i_def), intent(in) :: map_multi(ndf_multi)
    integer(kind=i_def), intent(in) :: ndf_high, undf_high
    integer(kind=i_def), intent(in) :: map_high(ndf_high)

    real(kind=r_def), intent(out) :: multi_field(undf_multi)
    real(kind=r_def), intent(in)  :: high_field(undf_high)

    integer(kind=i_def) :: i

    ! Convert high-order field to multi-data field
    do i = 1, ndata
      multi_field(map_multi(1)+i-1) = high_field(map_high(i))
    end do

  end subroutine high_to_multi_code

end module high_to_multi_kernel_mod
