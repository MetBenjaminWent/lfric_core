!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Compute the q3t matrix for analytic elimination of theta. The family
!!        of qXY matrices come from elimination of theta from the mixed solver.
!> @details Operator to map the residual from the thermodynamic equation
!!          to the equation of state: q3t =  inv_m3 * < sigma, detJ/theta * w >
!!          where w is a basis function in the Wtheta space and
!!          sigma is a test function in the W3 space.
!!          For more details, see the solver section of
!!          https://code.metoffice.gov.uk/trac/lfric/wiki/GhaspSupport/Documentation

module eliminated_theta_q3t_kernel_mod

  use argument_mod,            only: arg_type, func_type,             &
                                     GH_OPERATOR, GH_FIELD, GH_REAL,  &
                                     GH_READ, GH_WRITE,               &
                                     GH_BASIS, GH_DIFF_BASIS,         &
                                     CELL_COLUMN, GH_QUADRATURE_XYoZ, &
                                     ANY_DISCONTINUOUS_SPACE_3
  use constants_mod,           only: r_def, i_def
  use coordinate_jacobian_mod, only: pointwise_coordinate_jacobian
  use fs_continuity_mod,       only: W3, Wtheta, Wchi
  use kernel_mod,              only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  type, public, extends(kernel_type) :: eliminated_theta_q3t_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                     &
        arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W3, Wtheta),               &
        arg_type(GH_FIELD,    GH_REAL, GH_READ,  Wtheta),                   &
        arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3, W3),                   &
        arg_type(GH_FIELD*3,  GH_REAL, GH_READ,  Wchi),                     &
        arg_type(GH_FIELD,    GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3) &
        /)
    type(func_type) :: meta_funcs(3) = (/                                   &
        func_type(W3,     GH_BASIS),                                        &
        func_type(Wtheta, GH_BASIS),                                        &
        func_type(Wchi,   GH_BASIS, GH_DIFF_BASIS)                          &
        /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: eliminated_theta_q3t_code
  end type eliminated_theta_q3t_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public eliminated_theta_q3t_code

contains

!> @brief Compute the q3t matrix thaty arises from analytic elimination of theta
!!        in the equation of state:
!!        q3t =  inv_m3 * < sigma, detJ/theta * w >.
!> @param[in]     cell           Horizontal cell index
!> @param[in]     nlayers        Number of layers
!> @param[in]     ncell_3d1      Number of cells in the 3D mesh
!> @param[in,out] q3t_op         Projection matrix
!> @param[in]     theta          Potential temperature
!> @param[in]     ncell_3d2      Number of cells in the 3D mesh
!> @param[in]     inv_m3         Inverse of W3 mass matrix
!> @param[in]     chi1           First component of the coordinate array
!> @param[in]     chi2           Second component of the coordinate array
!> @param[in]     chi3           Third component of the coordinate array
!> @param[in]     panel_id       Field containing the panel ID indicator
!> @param[in]     ndf_w3         Number of degrees of freedom per cell W3 space
!> @param[in]     basis_w3       Basis functions for the W3 space evaluated at quadrature points
!> @param[in]     ndf_wt         Number of degrees of freedom per cell for the theta space
!> @param[in]     undf_wt        Total number of degrees of freedom for the theta space
!> @param[in]     map_wt         Dofmap for the bottom layer in the theta space
!> @param[in]     basis_wt       Basis functions evaluated at quadrature points
!> @param[in]     ndf_chi        Number of degrees of freedom per cell for the coordinate field
!> @param[in]     undf_chi       Number of unique degrees of freedom for chi field
!> @param[in]     map_chi        Dofmap for the cell at the base of the column
!> @param[in]     basis_chi      Wchi basis functions evaluated at quadrature points
!> @param[in]     diff_basis_chi Wchi differential basis functions evaluated at quadrature points
!> @param[in]     ndf_pid        Number of degrees of freedom per cell for panel_id
!> @param[in]     undf_pid       Number of unique degrees of freedom for panel_id
!> @param[in]     map_pid        Dofmap for the cell at the base of the column for panel_id
!> @param[in]     nqp_h          Number of horizontal quadrature points
!> @param[in]     nqp_v          Number of vertical quadrature points
!> @param[in]     wqp_h          Horizontal quadrature weights
!> @param[in]     wqp_v          Vertical quadrature weights
subroutine eliminated_theta_q3t_code(cell, nlayers,                      &
                                     ncell_3d1, q3t_op,                  &
                                     theta,                              &
                                     ncell_3d2, inv_m3,                  &
                                     chi1, chi2, chi3,                   &
                                     panel_id,                           &
                                     ndf_w3, basis_w3,                   &
                                     ndf_wt, undf_wt, map_wt, basis_wt,  &
                                     ndf_chi, undf_chi,                  &
                                     map_chi, basis_chi, diff_basis_chi, &
                                     ndf_pid, undf_pid, map_pid,         &
                                     nqp_h, nqp_v, wqp_h, wqp_v)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: cell, nqp_h, nqp_v
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w3, ndf_chi, ndf_wt, ndf_pid
  integer(kind=i_def), intent(in) :: undf_chi, undf_wt, undf_pid
  integer(kind=i_def), intent(in) :: ncell_3d1, ncell_3d2

  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid
  integer(kind=i_def), dimension(ndf_wt),  intent(in) :: map_wt

  real(kind=r_def), dimension(ndf_w3,ndf_wt,ncell_3d1), intent(inout) :: q3t_op
  real(kind=r_def), dimension(ndf_w3,ndf_w3,ncell_3d2), intent(in)    :: inv_m3

  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v), intent(in) :: basis_chi
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: diff_basis_chi
  real(kind=r_def), dimension(1,ndf_w3, nqp_h,nqp_v), intent(in) :: basis_w3
  real(kind=r_def), dimension(1,ndf_wt, nqp_h,nqp_v), intent(in) :: basis_wt

  real(kind=r_def), dimension(undf_wt),  intent(in) :: theta
  real(kind=r_def), dimension(undf_chi), intent(in) :: chi1
  real(kind=r_def), dimension(undf_chi), intent(in) :: chi2
  real(kind=r_def), dimension(undf_chi), intent(in) :: chi3
  real(kind=r_def), dimension(undf_pid), intent(in) :: panel_id

  real(kind=r_def), dimension(nqp_h), intent(in) :: wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) :: wqp_v

  ! Internal variables
  integer(kind=i_def)                  :: df, df1, df2, k, ik
  integer(kind=i_def)                  :: qp1, qp2
  real(kind=r_def), dimension(ndf_chi) :: chi1_e, chi2_e, chi3_e
  real(kind=r_def)                     :: theta_quad
  real(kind=r_def)                     :: integrand
  real(kind=r_def), dimension(3,3)     :: jac
  real(kind=r_def)                     :: dj

  integer(kind=i_def) :: ipanel

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
    do df = 1, ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do
    ik = 1 + k + (cell-1)*nlayers
    q3t_op(:,:,ik) = 0.0_r_def
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        call pointwise_coordinate_jacobian(ndf_chi, chi1_e, chi2_e, chi3_e, &
                                           ipanel, basis_chi(:,:,qp1,qp2),  &
                                           diff_basis_chi(:,:,qp1,qp2),     &
                                           jac, dj                          )

        theta_quad = 0.0_r_def
        do df = 1,ndf_wt
          theta_quad = theta_quad + theta(map_wt(df)+k)*basis_wt(1,df,qp1,qp2)
        end do
        integrand = wqp_h(qp1)*wqp_v(qp2)/theta_quad*dj
        do df2 = 1, ndf_wt
          do df1 = 1, ndf_w3
            q3t_op(df1,df2,ik) = q3t_op(df1,df2,ik)                &
                               + integrand*basis_w3(1,df1,qp1,qp2) &
                                          *basis_wt(1,df2,qp1,qp2)
          end do
        end do
      end do
    end do
    q3t_op(:,:,ik) = matmul(inv_m3(:,:,ik),q3t_op(:,:,ik))
  end do

end subroutine eliminated_theta_q3t_code

end module eliminated_theta_q3t_kernel_mod
