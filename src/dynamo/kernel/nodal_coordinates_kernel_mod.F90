!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

! Kernel to compute the coordinates fields at nodal points of another
! function space. In general this will give a polynomial approximation to the
! mesh, i.e. if the mesh is spherical then the nodal coordinates will not lie on
! spherical shells but will instead represent a polynomial approximation to the
! spherical shell.

!> @detail 
module nodal_coordinates_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,                     &
                                    GH_FIELD, GH_READ, GH_WRITE,             &
                                    W0, ANY_SPACE_1,                         &
                                    GH_BASIS,                                &
                                    CELLS
use constants_mod,           only : r_def

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: nodal_coordinates_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_FIELD*3,   GH_WRITE, ANY_SPACE_1),                  &
       arg_type(GH_FIELD*3,   GH_READ,  W0)                            &
       /)
  type(func_type) :: meta_funcs(1) = (/                                &
       func_type(W0, GH_BASIS)                                         &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass ::nodal_coordinates_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface nodal_coordinates_kernel_type
   module procedure nodal_coordinates_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public nodal_coordinates_code
contains

type(nodal_coordinates_kernel_type) function nodal_coordinates_kernel_constructor() result(self)
  return
end function nodal_coordinates_kernel_constructor

!nlayers Integer the number of layers
!ndf_x The number of degrees of freedom per cell for the output field
!undf_x The number of unique degrees of freedom for the output field
!map_x Integer array holding the dofmap for the cell at the base of the column for the output field
!ndf_chi The number of degrees of freedom per cell for the input field
!undf_chi The number of unique degrees of freedom for the input field
!map_chi Integer array holding the dofmap for the cell at the base of the column for the input field
!nodal_x the array of nodal coordinates in the first direction
!nodal_y the array of nodal coordinates in the second direction
!nodal_z the array of nodal coordinates in the third direction
!chi1 the array of fem coordinates in the first direction
!chi2 the array of fem coordinates in the second direction
!chi3 the array of fem coordinates in the third direction
!basis_chi the basis functions of the chi function space evaluated at the nodal points of the x function space
subroutine nodal_coordinates_code(nlayers,                                    &
                                  nodal_x, nodal_y, nodal_z,                  &
                                  chi1, chi2, chi3,                           &
                                  ndf_x, undf_x, map_x,                       &
                                  ndf_chi, undf_chi, map_chi,                 &
                                  basis_chi                                   &
                                  )
                           
  !Arguments
  integer, intent(in) :: nlayers
  integer, intent(in) :: ndf_x, ndf_chi, undf_x, undf_chi
  integer, dimension(ndf_x),   intent(in) :: map_x
  integer, dimension(ndf_chi), intent(in) :: map_chi
  real(kind=r_def), dimension(undf_x),        intent(out) :: nodal_x, nodal_y, nodal_z
  real(kind=r_def), dimension(undf_chi),      intent(in)  :: chi1, chi2, chi3
  real(kind=r_def), dimension(1,ndf_chi,ndf_x), intent(in)  :: basis_chi

  !Internal variables
  integer          :: df_x, df_chi, k
  real(kind=r_def) :: xyz(3)
  
  do k = 0, nlayers-1
    do df_x = 1,ndf_x
      xyz(:) = 0.0_r_def
      do df_chi = 1, ndf_chi
        xyz(1) = xyz(1) + chi1(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
        xyz(2) = xyz(2) + chi2(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
        xyz(3) = xyz(3) + chi3(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
      end do
      nodal_x(map_x(df_x)+k) = xyz(1)
      nodal_y(map_x(df_x)+k) = xyz(2)
      nodal_z(map_x(df_x)+k) = xyz(3)
    end do
  end do

end subroutine nodal_coordinates_code

end module nodal_coordinates_kernel_mod
