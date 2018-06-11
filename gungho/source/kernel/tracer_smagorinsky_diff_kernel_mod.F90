!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Applies horizontal Smagorinsky diffusion mu * (d2dx2 + d2dy2) to a tracer
!>        variable in the Wtheta space for lowest order elements. Initially using
!>        constant viscosity mu but will eventually use the blended BL-Smagorinsky
!>        diffusion coefficient.
!>
module tracer_smagorinsky_diff_kernel_mod

  use argument_mod,          only : arg_type, func_type,         &
                                    GH_FIELD, GH_READ, GH_WRITE, &
                                    CELLS, STENCIL, CROSS
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : Wtheta
  use kernel_mod,            only : kernel_type
  use mixing_config_mod,     only : viscosity_mu

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: tracer_smagorinsky_diff_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                        &
        arg_type(GH_FIELD,   GH_WRITE,  Wtheta),               &
        arg_type(GH_FIELD,   GH_READ, Wtheta, STENCIL(CROSS)), &
        arg_type(GH_FIELD,   GH_READ, Wtheta)                  &
        /)
    integer :: iterates_over = CELLS
  contains
    procedure, nopass ::tracer_smagorinsky_diff_code
  end type

  !---------------------------------------------------------------------------
  ! Constructors
  !---------------------------------------------------------------------------

  ! Overload the default structure constructor for function space
  interface tracer_smagorinsky_diff_kernel_type
    module procedure tracer_smagorinsky_diff_kernel_constructor
  end interface

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public tracer_smagorinsky_diff_code

contains

type(tracer_smagorinsky_diff_kernel_type) function tracer_smagorinsky_diff_kernel_constructor() result(self)
  implicit none
  return
end function tracer_smagorinsky_diff_kernel_constructor

!> @brief Calculates horizontal Smagorinsky diffusion for a tracer variable
!! @param[in] nlayers Number of layers in the mesh
!! @param[in] theta_inc Diffusion increment for temperature field
!! @param[in] theta_n Input temperature field
!! @param[in] map_wt_size Number of cells in the stencil at the base of the column for wt
!! @param[in] map_wt Array holding the dofmap for the stencil at the base of the column for wt
!! @param[in] delta Edge length on wtheta points
!! @param[in] ndf_wt Number of degrees of freedom per cell for theta space
!! @param[in] undf_wt  Number of unique degrees of freedom  for theta_space
!! @param[in] cell_map_wt Cell dofmap for the wtheta space

subroutine tracer_smagorinsky_diff_code(nlayers,                        &
                                        theta_inc, theta_n,             &
                                        map_wt_size, map_wt,            &
                                        delta,                          &
                                        ndf_wt, undf_wt, cell_map_wt)

  implicit none
  ! Arguments
  integer(kind=i_def), intent(in)                                 :: nlayers
  integer(kind=i_def), intent(in)                                 :: ndf_wt, undf_wt
  integer(kind=i_def), intent(in)                                 :: map_wt_size
  integer(kind=i_def), dimension(ndf_wt,map_wt_size), intent(in)  :: map_wt
  integer(kind=i_def), dimension(ndf_wt),             intent(in)  :: cell_map_wt

  real(kind=r_def), dimension(undf_wt),  intent(inout) :: theta_inc
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: theta_n
  real(kind=r_def), dimension(undf_wt),  intent(in)    :: delta

  ! Internal variables
  integer(kind=i_def) :: k
  real(kind=r_def)    :: d2dx, d2dy
  real(kind=r_def)    :: idx2

  !  ---------- 
  !  |    |   |
  !  |  w | i |
  !  -----x----
  !  |    |   |
  !  | sw | s |
  !  ----------
  !  y
  !  ^
  !  |_> x
  !

  ! Horizontal theta diffusion 

  do k = 0, nlayers
    idx2 = 1.0_r_def/(delta(map_wt(1,1) + k))**2
    d2dx = (theta_n(map_wt(1,2) + k)  - 2.0_r_def*theta_n(map_wt(1,1) + k) + theta_n(map_wt(1,4) + k) )*idx2
    d2dy = (theta_n(map_wt(1,3) + k)  - 2.0_r_def*theta_n(map_wt(1,1) + k) + theta_n(map_wt(1,5) + k) )*idx2
    theta_inc(cell_map_wt(1) + k) = viscosity_mu*(d2dx + d2dy)
  end do

end subroutine tracer_smagorinsky_diff_code

end module tracer_smagorinsky_diff_kernel_mod
