!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module compute_p0_vert_precon_kernel_mod

  use argument_mod,      only : arg_type,                       &
                                GH_FIELD, GH_OPERATOR, GH_REAL, &
                                GH_READ, GH_WRITE,              &
                                CELLS
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W2, W3, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: compute_p0_vert_precon_kernel_type
    private
    type(arg_type) :: meta_args(10) = (/                  &
        arg_type(GH_FIELD*3,  GH_WRITE, W3),              & ! tri_precon
        arg_type(GH_OPERATOR, GH_READ,  W2,     W3),      & ! div_star
        arg_type(GH_FIELD,    GH_READ,  W2),              & ! hb_lumped_inv
        arg_type(GH_FIELD,    GH_READ,  W2),              & ! w_normalisation
        arg_type(GH_FIELD,    GH_READ,  Wtheta),          & ! t_normalisation
        arg_type(GH_OPERATOR, GH_READ,  W3,     W2),      & ! compound_div
        arg_type(GH_OPERATOR, GH_READ,  W3,     Wtheta),  & ! p3theta
        arg_type(GH_OPERATOR, GH_READ,  Wtheta, W2),      & ! ptheta2v
        arg_type(GH_OPERATOR, GH_READ,  W3,     W3),      & ! m3_exner_star
        arg_type(GH_OPERATOR, GH_READ,  W3,     W3)       & ! m3_inv
        /)
    integer :: iterates_over = CELLS
  contains
    procedure, nopass :: compute_p0_vert_precon_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public compute_p0_vert_precon_code

contains


  !> @brief Compute the terms of the helmholz operator restricted to a single
  !!        column
  !> @details Computes the coefficients of the tri diagonal matrix
  !!          [ tri_m, tri_0, tri_p] that correspond to the helmholtz operator
  !!          for lowest order elements
  !!          restricted to a single column. This is then used as the
  !!          preconditioner. This routine is a condensed version of the
  !!          algorithm level calls to apply the Helmholtz operator for lowest
  !!          order elements. Since it is restricted to a single column it can be
  !!          applied in a single kernel.
  !> @param[in]  cell Horizontal cell index
  !> @param[in]  nlayers Number of layers
  !> @param[in,out] tri_0 Diagonal entry to preconditioner matrix
  !> @param[in,out] tri_p Upper diagonal entry to preconditioner matrix
  !> @param[in,out] tri_m Lower diagonal entry to preconditioner matrix
  !> @param[in]  ncell_3d_1 Total number of cells for divergence matrix
  !> @param[in]  div_star Weighted transpose of the divergence operator
  !> @param[in]  hb_lumped_inv Lumped inverse of the HB (mass matrix + buoyancy)
  !!             term
  !> @param[in]  u_normalisation Field use to normalise the momentum equation
  !> @param[in]  t_normalisation Field use to normalise the thermodynamic equation
  !> @param[in]  ncell_3d_2 Total number of cells for compound_div matrix
  !> @param[in]  compound_div divergence operator weighted by reference density
  !!             and mass matrices
  !> @param[in]  ncell_3d_3 Total number of cells for p3t matrix
  !> @param[in]  p3theta Weighted projection operator from Wtheta to W3
  !> @param[in]  ncell_3d_4 Total number of cells for ptheta2v matrix
  !> @param[in]  ptheta2v Weighted projection operator from the vertical
  !!             components of W2 to Wtheta
  !> @param[in]  ncell_3d_5 Total number of cells for m3_exner_star matrix
  !> @param[in]  m3_exner_star weighted W3 mass matrix
  !> @param[in]  ncell_3d_6 Total number of cells for m3_inv matrix
  !> @param[in]  m3_inv Inverse of W3 mass matrix
  !> @param[in]  ndf_w3 Number of degrees of freedom per cell for the pressure space
  !> @param[in]  undf_w3 Unique number of degrees of freedom  for the pressure space
  !> @param[in]  map_w3 Dofmap for the cell at the base of the column for the pressure space
  !> @param[in]  ndf_w2 Number of degrees of freedom per cell for the velocity space
  !> @param[in]  undf_w2 Unique number of degrees of freedom  for the velocity space
  !> @param[in]  map_w2 Dofmap for the cell at the base of the column for the velocity space
  !> @param[in]  ndf_wt Number of degrees of freedom per cell for the temperature space
  !> @param[in]  undf_wt Unique number of degrees of freedom  for the temperature space
  !> @param[in]  map_wt Dofmap for the cell at the base of the column for the temperature space
  subroutine compute_p0_vert_precon_code(cell,                    &
                                         nlayers,                 &
                                         tri_0, tri_p, tri_m,     &
                                         ncell_3d_1,              &
                                         div_star,                &
                                         hb_lumped_inv,           &
                                         u_normalisation,         &
                                         t_normalisation,         &
                                         ncell_3d_2,              &
                                         compound_div,            &
                                         ncell_3d_3,              &
                                         p3theta,                 &
                                         ncell_3d_4,              &
                                         ptheta2v,                &
                                         ncell_3d_5,              &
                                         m3_exner_star,           &
                                         ncell_3d_6,              &
                                         m3_inv,                  &
                                         ndf_w3, undf_w3, map_w3, &
                                         ndf_w2, undf_w2, map_w2, &
                                         ndf_wt, undf_wt, map_wt)

  implicit none
  ! Arguments
  integer(kind=i_def),                    intent(in) :: cell, nlayers
  integer(kind=i_def),                    intent(in) :: ncell_3d_1, ncell_3d_2
  integer(kind=i_def),                    intent(in) :: ncell_3d_3, ncell_3d_4
  integer(kind=i_def),                    intent(in) :: ncell_3d_5, ncell_3d_6
  integer(kind=i_def),                    intent(in) :: undf_w2, ndf_w2
  integer(kind=i_def),                    intent(in) :: undf_w3, ndf_w3
  integer(kind=i_def),                    intent(in) :: undf_wt, ndf_wt
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt

  real(kind=r_def), dimension(undf_w3), intent(inout) :: tri_0, tri_p, tri_m
  real(kind=r_def), dimension(undf_w2), intent(in)    :: hb_lumped_inv
  real(kind=r_def), dimension(undf_w2), intent(in)    :: u_normalisation
  real(kind=r_def), dimension(undf_wt), intent(in)    :: t_normalisation

  real(kind=r_def), dimension(ndf_w2, ndf_w3, ncell_3d_1), intent(in) :: div_star
  real(kind=r_def), dimension(ndf_w3, ndf_w2, ncell_3d_2), intent(in) :: compound_div
  real(kind=r_def), dimension(ndf_w3, ndf_wt, ncell_3d_3), intent(in) :: p3theta
  real(kind=r_def), dimension(ndf_wt, ndf_w2, ncell_3d_4), intent(in) :: ptheta2v
  real(kind=r_def), dimension(ndf_w3, ndf_w3, ncell_3d_5), intent(in) :: m3_exner_star
  real(kind=r_def), dimension(ndf_w3, ndf_w3, ncell_3d_6), intent(in) :: m3_inv

  ! Internal variables
  integer(kind=i_def) :: k, ik, df

  real(kind=r_def), dimension(6,0:nlayers-1) :: grad_p
  real(kind=r_def), dimension(0:nlayers)     :: t_e
  real(kind=r_def)                           :: div_u, t_at_p

  ! Compute weights of helmholtz operator in a single column and store them in
  ! the tridiagonal matrix
  ! Hard wired optimisation for desired configuration (p=0 elements with pt2
  ! only acting on vertical components of grad_p )

  ! First Stage: Compute hb_lumped_inv*u_normalisation*div_star ~ grad
  do k = 0, nlayers-1
    ik = (cell-1)*nlayers + k + 1
    do df = 1,6
      grad_p(df,k) = div_star(df,1,ik)*hb_lumped_inv(map_w2(df)+k) &
                                      *u_normalisation(map_w2(df)+k)
    end do
  end do
  ! Apply zero flux boundary conditions
  grad_p(5,0)         = 0.0_r_def
  grad_p(6,nlayers-1) = 0.0_r_def

  ! Second stage: Apply pointwise mapping of grad_p to t
  t_e = 0.0_r_def
  do k = 0,nlayers-1
    ik = (cell-1)*nlayers + k + 1
    t_e(k)   = t_e(k)   + ptheta2v(1,5,ik)*grad_p(5,k)
    t_e(k+1) = t_e(k+1) + ptheta2v(2,6,ik)*grad_p(6,k)
  end do

  ! Third stage: Compute D * grad_p + P3t * Mt^-1 * ( Pt2 * grad_p )

  ! Diagonal term
  do k = 0,nlayers-1
    ik = (cell-1)*nlayers + k + 1

    t_at_p = t_normalisation(map_wt(1)+k)*p3theta(1,1,ik)*t_e(k) &
           + t_normalisation(map_wt(2)+k)*p3theta(1,2,ik)*t_e(k+1)

    div_u = compound_div(1,1,ik)*grad_p(1,k) + compound_div(1,2,ik)*grad_p(2,k)&
          + compound_div(1,3,ik)*grad_p(3,k) + compound_div(1,4,ik)*grad_p(4,k)&
          + compound_div(1,5,ik)*grad_p(5,k) + compound_div(1,6,ik)*grad_p(6,k)
    tri_0(map_w3(1)+k) = m3_inv(1,1,ik)*(m3_exner_star(1,1,ik) + div_u + t_at_p)
  end do

  ! Upper diagonal term
  do k = 0,nlayers-2
    ik = (cell-1)*nlayers + k + 1

    t_at_p = t_normalisation(map_wt(2)+k)*p3theta(1,2,ik)&
            *(ptheta2v(2,6,ik)+ptheta2v(1,5,ik+1))*grad_p(5,k+1)

    div_u = compound_div(1,6,ik)*grad_p(5,k+1)
    tri_p(map_w3(1)+k) = m3_inv(1,1,ik)*(div_u + t_at_p)
  end do
  k = nlayers-1
  tri_p(map_w3(1)+k) = 0.0_r_def

  ! Lower diagonal term
  k = 0
  tri_m(map_w3(1)+k) = 0.0_r_def
  do k = 1,nlayers-1
    ik = (cell-1)*nlayers + k + 1

    t_at_p = t_normalisation(map_wt(1)+k)*p3theta(1,1,ik)&
            *(ptheta2v(1,5,ik)+ptheta2v(2,6,ik-1))*grad_p(6,k-1)

    div_u = compound_div(1,5,ik)*grad_p(6,k-1)
    tri_m(map_w3(1)+k) = m3_inv(1,1,ik)*(div_u + t_at_p)
  end do

end subroutine compute_p0_vert_precon_code

end module compute_p0_vert_precon_kernel_mod
