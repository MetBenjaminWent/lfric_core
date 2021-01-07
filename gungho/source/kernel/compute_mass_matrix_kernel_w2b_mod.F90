!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Provides access to the members of the w2b_kernel class.
module compute_mass_matrix_kernel_w2b_mod

  use argument_mod,            only: arg_type, func_type,       &
                                     GH_OPERATOR, GH_FIELD,     &
                                     GH_READ, GH_WRITE,         &
                                     ANY_SPACE_9,               &
                                     GH_BASIS, GH_DIFF_BASIS,   &
                                     CELLS, GH_QUADRATURE_XYoZ, &
                                     ANY_DISCONTINUOUS_SPACE_3
  use constants_mod,           only: i_def, r_def
  use coordinate_jacobian_mod, only: coordinate_jacobian
  use fs_continuity_mod,       only: W2broken
  use kernel_mod,              only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: compute_mass_matrix_kernel_w2b_type
    private
    type(arg_type) :: meta_args(3) = (/                            &
        arg_type(GH_OPERATOR, GH_WRITE, W2broken, W2broken),       &
        arg_type(GH_FIELD*3,  GH_READ,  ANY_SPACE_9),              &
        arg_type(GH_FIELD,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3) &
        /)
    type(func_type) :: meta_funcs(2) = (/                          &
        func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS),           &
        func_type(W2broken, GH_BASIS)                              &
        /)
    integer :: iterates_over = CELLS
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: compute_mass_matrix_w2b_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public compute_mass_matrix_w2b_code

contains

  !> @brief Computes the mass matrix for the w2b space.
  !!
  !! @param[in] cell      Identifying number of cell.
  !! @param[in] nlayers   Number of layers.
  !! @param[in] ncell_3d  ncell*ndf
  !! @param[in,out] mm    Local stencil or mass matrix.
  !! @param[in] chi1      1st (spherical) coordinate field in Wchi
  !! @param[in] chi2      2nd (spherical) coordinate field in Wchi
  !! @param[in] chi3      3rd (spherical) coordinate field in Wchi
  !! @param[in] panel_id  Field giving the ID for mesh panels.
  !! @param[in] ndf_w2b   Degrees of freedom per cell for W2broken space.
  !! @param[in] basis_w2b Vector basis functions evaluated at quadrature points
  !!                      for W2broken space.
  !! @param[in] ndf_chi   Degrees of freedom per cell for chi field.
  !! @param[in] undf_chi  Unique degrees of freedom for chi field.
  !! @param[in] map_chi   Dofmap for the cell at the base of the column, for the
  !!                      space on which the chi field lives.
  !! @param[in] basis_chi Basis functions evaluated at quadrature points for the
  !!                      space on which the chi field lives.
  !! @param[in] diff_basis_chi Vector differential basis functions evaluated at
  !!                      quadrature points for the space on which the chi
  !!                      field lives.
  !! @param[in] ndf_pid   Number of degrees of freedom per cell for panel_id
  !! @param[in] undf_pid  Number of unique degrees of freedom for panel_id
  !! @param[in] map_pid   Dofmap for the cell at the base of the column for panel_id
  !! @param[in] nqp_h     Number of horizontal quadrature points.
  !! @param[in] nqp_v     Number of vertical quadrature points.
  !! @param[in] wqp_h     Horizontal quadrature weights.
  !! @param[in] wqp_v     Vertical quadrature weights.
  subroutine compute_mass_matrix_w2b_code(cell, nlayers, ncell_3d,     &
                                          mm,                          &
                                          chi1, chi2, chi3, panel_id,  &
                                          ndf_w2b, basis_w2b,          &
                                          ndf_chi, undf_chi, map_chi,  &
                                          basis_chi, diff_basis_chi,   &
                                          ndf_pid, undf_pid, map_pid,  &
                                          nqp_h, nqp_v, wqp_h, wqp_v)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: cell, nqp_h, nqp_v
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ncell_3d
    integer(kind=i_def), intent(in)    :: ndf_w2b
    integer(kind=i_def), intent(in)    :: ndf_chi
    integer(kind=i_def), intent(in)    :: undf_chi
    integer(kind=i_def), intent(in)    :: map_chi(ndf_chi)
    integer(kind=i_def), intent(in)    :: ndf_pid
    integer(kind=i_def), intent(in)    :: undf_pid
    integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)

    real(kind=r_def),    intent(inout) :: mm(ndf_w2b,ndf_w2b,ncell_3d)
    real(kind=r_def),    intent(in)    :: basis_chi(1,ndf_chi,nqp_h,nqp_v)
    real(kind=r_def),    intent(in)    :: diff_basis_chi(3,ndf_chi,nqp_h,nqp_v)
    real(kind=r_def),    intent(in)    :: chi1(undf_chi)
    real(kind=r_def),    intent(in)    :: chi2(undf_chi)
    real(kind=r_def),    intent(in)    :: chi3(undf_chi)
    real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
    real(kind=r_def),    intent(in)    :: wqp_h(nqp_h)
    real(kind=r_def),    intent(in)    :: wqp_v(nqp_v)

    real(kind=r_def), dimension(3, ndf_w2b, nqp_h, nqp_v), intent(in) :: basis_w2b

    ! Internal variables
    integer(kind=i_def)                             :: df, df2, k, ik, ipanel
    integer(kind=i_def)                             :: qp1, qp2

    real(kind=r_def), dimension(ndf_chi)            :: chi1_e, chi2_e, chi3_e
    real(kind=r_def)                                :: integrand
    real(kind=r_def), dimension(nqp_h, nqp_v)       :: dj
    real(kind=r_def), dimension(3, 3, nqp_h, nqp_v) :: jac

    ipanel = int(panel_id(map_pid(1)), i_def)

    ! Loop over layers: Start from 1 as in this loop k is not an offset
    do k = 1, nlayers
      ik = k + (cell - 1) * nlayers

      ! Indirect the chi coord field here
      do df = 1, ndf_chi
        chi1_e(df) = chi1(map_chi(df) + k - 1)
        chi2_e(df) = chi2(map_chi(df) + k - 1)
        chi3_e(df) = chi3(map_chi(df) + k - 1)
      end do

      ! Compute Jacobian
      call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,  &
                               ipanel, basis_chi, diff_basis_chi, jac, dj)

      ! Run over dof extent of W2Broken
      do df2 = 1, ndf_w2b
        do df = df2, ndf_w2b ! Exploiting symmetry of the operator
          ! Initialize
          mm(df, df2, ik) = 0.0_r_def
          do qp2 = 1, nqp_v
            do qp1 = 1, nqp_h
              ! Mass term: dot(w2b_basis_i, w2b_basis_j) * dx (mapped to reference space)
              integrand = wqp_h(qp1) * wqp_v(qp2)                                               &
                        * dot_product(matmul(jac(:, :, qp1, qp2), basis_w2b(:, df, qp1, qp2)),  &
                                      matmul(jac(:, :, qp1, qp2), basis_w2b(:, df2, qp1, qp2))) &
                        / dj(qp1, qp2)
              mm(df, df2, ik) = mm(df, df2, ik) + integrand
            end do
          end do
        end do
        do df = df2, 1, -1
          ! Mass matrix is symmetric
          mm(df, df2, ik) = mm(df2, df, ik)
        end do
      end do

    end do ! End of layer loop

  end subroutine compute_mass_matrix_w2b_code

end module compute_mass_matrix_kernel_w2b_mod
