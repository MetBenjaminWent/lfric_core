!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Transpose a matrix in local operator representation
!> @details Given a matrix \f$A\f$ calculate the transpose \f$B=A^T\f$

module transpose_matrix_kernel_mod
use argument_mod,            only : arg_type,                                 &
                                    GH_FIELD, GH_OPERATOR, GH_READ, GH_WRITE, &
                                    ANY_SPACE_1, ANY_SPACE_2,                 &
                                    CELLS
use constants_mod,           only : r_def
use kernel_mod,              only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: transpose_matrix_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_OPERATOR, GH_READ,  ANY_SPACE_1, ANY_SPACE_2),      &
       arg_type(GH_OPERATOR, GH_WRITE, ANY_SPACE_2, ANY_SPACE_1)       &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass ::transpose_matrix_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public transpose_matrix_code
contains

!> @brief Computes the transpose of a matrix
!> @param[in]  cell Horizontal cell index
!> @param[in]  nlayers Number of layers
!> @param[in]  ncell_3d Total number of cells
!> @param[in]  mat_in Input matrix
!> @param[in]  ncell_3d_2 Total number of cells (passed in twice)
!> @param[in,out] mat_out Resulting transposed matrix
!> @param[in]  ndf1 Number of degrees of freedom per cell for space 1
!> @param[in]  ndf2 Number of degrees of freedom per cell for space 2
subroutine transpose_matrix_code(cell,        &
                                 nlayers,     &
                                 ncell_3d,    &
                                 mat_in,      &
                                 ncell_3d_2,  &
                                 mat_out,     &
                                 ndf1,        &
                                 ndf2)

  implicit none

  !Arguments
  integer,                   intent(in)    :: cell
  integer,                   intent(in)    :: nlayers
  integer,                   intent(in)    :: ncell_3d
  integer,                   intent(in)    :: ncell_3d_2
  integer,                   intent(in)    :: ndf1
  integer,                   intent(in)    :: ndf2
  real(kind=r_def), dimension(ndf1,ndf2,ncell_3d), intent(in)      :: mat_in
  real(kind=r_def), dimension(ndf2,ndf1,ncell_3d), intent(inout)   :: mat_out

  !Internal variables
  integer                           :: k, ik

  do k = 0, nlayers-1
    ik = (cell-1)*nlayers + k + 1
    mat_out(:,:,ik) = transpose(mat_in(:,:,ik))
  end do

end subroutine transpose_matrix_code

end module transpose_matrix_kernel_mod
