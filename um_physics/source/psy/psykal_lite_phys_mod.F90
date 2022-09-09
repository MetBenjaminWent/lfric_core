!----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Provides an implementation of the Psy layer for physics

!> @details Contains hand-rolled versions of the Psy layer that can be used for
!> simple testing and development of the scientific code

module psykal_lite_phys_mod

  use constants_mod,         only : i_def, r_def
  use field_mod,             only : field_type, field_proxy_type
  use mesh_mod,              only : mesh_type

  implicit none
  public

contains
  !---------------------------------------------------------------------
  !> LFRic and PSyclone currently do not have a mechanism to loop over a subset
  !> of cells in a horizontal domain. The relevant PSyclone ticket relating to
  !> this is #487.
  !> The orographic drag kernel only needs to be applied to a subset of land
  !> points where the standard deviation of subgrid orography is more than zero.
  !>
  !> invoke_orographic_drag_kernel: Invokes the kernel which calls the UM
  !> orographic drag scheme only on points where the standard deviation
  !> of subgrid orography is more than zero.
  subroutine invoke_orographic_drag_kernel(                              &
                      du_orog_blk, dv_orog_blk, du_orog_gwd, dv_orog_gwd,&
                      dtemp_orog_blk, dtemp_orog_gwd, u_in_w3, v_in_w3,&
                      wetrho_in_w3, theta_in_wth, exner_in_w3, sd_orog,  &
                      grad_xx_orog, grad_xy_orog, grad_yy_orog,          &
                      mr_v, mr_cl, mr_ci,                                &
                      height_w3, height_wth,                             &
                      taux_orog_blk, tauy_orog_blk,                      &
                      taux_orog_gwd, tauy_orog_gwd )

    use orographic_drag_kernel_mod, only: orographic_drag_kernel_code
    use mesh_mod, only: mesh_type
    implicit none

    ! Increments from orographic drag
    type(field_type), intent(inout) :: du_orog_blk, dv_orog_blk, & ! Winds
                                       du_orog_gwd, dv_orog_gwd, &
                                       dtemp_orog_blk, dtemp_orog_gwd ! Temperature

    ! Inputs to orographic drag scheme
    type(field_type), intent(in) :: u_in_w3, v_in_w3,           & ! Winds
                                    wetrho_in_w3, theta_in_wth, & ! Density, Temperature
                                    exner_in_w3,                & ! Exner pressure
                                    sd_orog, grad_xx_orog,      & ! Orography ancils
                                    grad_xy_orog, grad_yy_orog, & !
                                    mr_v, mr_cl, mr_ci,         & ! mixing ratios
                                    height_w3, height_wth         ! Heights
    ! Diagnostics from orographic drag
    ! Stress from ...
    type(field_type), intent(inout) :: taux_orog_blk, tauy_orog_blk, & ! ... orographic flow blocking drag
                                       taux_orog_gwd, tauy_orog_gwd    ! ... orographic gravity wave drag

    integer :: cell

    ! Number of degrees of freedom
    integer :: ndf_w3, undf_w3, ndf_wtheta, undf_wtheta

    ! These are in ANY_DISCONTINUOUS_SPACE_1
    integer :: ndf_adspc1_sd_orog, undf_adspc1_sd_orog

    integer :: nlayers

    type(field_proxy_type) :: du_blk_proxy, dv_blk_proxy,             &
                              du_orog_gwd_proxy, dv_orog_gwd_proxy,   &
                              dtemp_blk_proxy, dtemp_orog_gwd_proxy,  &
                              u_in_w3_proxy, v_in_w3_proxy,           &
                              wetrho_in_w3_proxy, theta_in_wth_proxy, &
                              exner_in_w3_proxy,                      &
                              sd_orog_proxy, grad_xx_orog_proxy,      &
                              grad_xy_orog_proxy, grad_yy_orog_proxy, &
                              mr_v_proxy, mr_cl_proxy, mr_ci_proxy,   &
                              height_w3_proxy, height_wth_proxy,      &
                              taux_blk_proxy, tauy_blk_proxy,         &
                              taux_gwd_proxy, tauy_gwd_proxy

    integer, pointer :: map_adspc1_sd_orog(:,:) => null(), &
                        map_w3(:,:) => null(),                  &
                        map_wtheta(:,:) => null()

    TYPE(mesh_type), pointer :: mesh => null()

    ! Initialise field and/or operator proxies
    du_blk_proxy = du_orog_blk%get_proxy()
    dv_blk_proxy = dv_orog_blk%get_proxy()
    du_orog_gwd_proxy = du_orog_gwd%get_proxy()
    dv_orog_gwd_proxy = dv_orog_gwd%get_proxy()
    dtemp_blk_proxy = dtemp_orog_blk%get_proxy()
    dtemp_orog_gwd_proxy = dtemp_orog_gwd%get_proxy()
    u_in_w3_proxy = u_in_w3%get_proxy()
    v_in_w3_proxy = v_in_w3%get_proxy()
    wetrho_in_w3_proxy = wetrho_in_w3%get_proxy()
    theta_in_wth_proxy = theta_in_wth%get_proxy()
    exner_in_w3_proxy = exner_in_w3%get_proxy()
    sd_orog_proxy = sd_orog%get_proxy()
    grad_xx_orog_proxy = grad_xx_orog%get_proxy()
    grad_xy_orog_proxy = grad_xy_orog%get_proxy()
    grad_yy_orog_proxy = grad_yy_orog%get_proxy()
    mr_v_proxy = mr_v%get_proxy()
    mr_cl_proxy = mr_cl%get_proxy()
    mr_ci_proxy = mr_ci%get_proxy()
    height_w3_proxy = height_w3%get_proxy()
    height_wth_proxy = height_wth%get_proxy()

    taux_blk_proxy = taux_orog_blk%get_proxy()
    tauy_blk_proxy = tauy_orog_blk%get_proxy()
    taux_gwd_proxy = taux_orog_gwd%get_proxy()
    tauy_gwd_proxy = tauy_orog_gwd%get_proxy()

    ! Initialise number of layers
    nlayers = du_blk_proxy%vspace%get_nlayers()

    ! Create a mesh object
    mesh => du_orog_blk%get_mesh()

    ! Look-up dofmaps for each function space
    map_w3 => du_blk_proxy%vspace%get_whole_dofmap()
    map_wtheta => dtemp_blk_proxy%vspace%get_whole_dofmap()
    map_adspc1_sd_orog => sd_orog_proxy%vspace%get_whole_dofmap()

    ! Initialise number of DoFs for w3
    ndf_w3 = du_blk_proxy%vspace%get_ndf()
    undf_w3 = du_blk_proxy%vspace%get_undf()

    ! Initialise number of DoFs for wtheta
    ndf_wtheta = dtemp_blk_proxy%vspace%get_ndf()
    undf_wtheta = dtemp_blk_proxy%vspace%get_undf()

    ! Initialise number of DoFs for adspc1_sd_orog
    ndf_adspc1_sd_orog = sd_orog_proxy%vspace%get_ndf()
    undf_adspc1_sd_orog = sd_orog_proxy%vspace%get_undf()

    ! Loop over cells
    do cell=1,mesh%get_last_edge_cell()
      ! Only call orographic_drag_kernel_code at points where the
      ! standard deviation of the subgrid orography is more than zero.
      if ( sd_orog_proxy%data(map_adspc1_sd_orog(1, cell)) > 0.0_r_def ) then

        call orographic_drag_kernel_code(                                  &
                      nlayers, du_blk_proxy%data, dv_blk_proxy%data,       &
                      du_orog_gwd_proxy%data, dv_orog_gwd_proxy%data,      &
                      dtemp_blk_proxy%data, dtemp_orog_gwd_proxy%data,     &
                      u_in_w3_proxy%data, v_in_w3_proxy%data,              &
                      wetrho_in_w3_proxy%data, theta_in_wth_proxy%data,    &
                      exner_in_w3_proxy%data,                              &
                      sd_orog_proxy%data, grad_xx_orog_proxy%data,         &
                      grad_xy_orog_proxy%data, grad_yy_orog_proxy%data,    &
                      mr_v_proxy%data, mr_cl_proxy%data, mr_ci_proxy%data, &
                      height_w3_proxy%data, height_wth_proxy%data,         &
                      taux_blk_proxy%data, tauy_blk_proxy%data,            &
                      taux_gwd_proxy%data, tauy_gwd_proxy%data,            &
                      ndf_w3, undf_w3, map_w3(:,cell),                     &
                      ndf_wtheta, undf_wtheta, map_wtheta(:,cell),         &
                      ndf_adspc1_sd_orog, undf_adspc1_sd_orog,   &
                      map_adspc1_sd_orog(:,cell) )

      end if ! sd_orog_proxy%data(map_adspc1_sd_orog(1, cell)) > 0.0_r_def

    end do


    ! Set halos dirty/clean for fields modified in the above loop
    call du_blk_proxy%set_dirty()
    call dv_blk_proxy%set_dirty()
    call du_orog_gwd_proxy%set_dirty()
    call dv_orog_gwd_proxy%set_dirty()
    call dtemp_blk_proxy%set_dirty()
    call dtemp_orog_gwd_proxy%set_dirty()
    call taux_blk_proxy%set_dirty()
    call tauy_blk_proxy%set_dirty()
    call taux_gwd_proxy%set_dirty()
    call tauy_gwd_proxy%set_dirty()

  end subroutine invoke_orographic_drag_kernel
  !---------------------------------------------------------------------
  !> Contains the Psy-layer to build the stochastic physics
  !> forcing pattern. At the moment it requires to pass the arrays
  !> stph_spectral_coeffc and stph_spectral_coeffs via the kernel argument.
  !> Psyclone does not recognize arrays in the argument yet
  !> this functionality is being developed in PSyclone ticket 1312
  !> at https://github.com/stfc/PSyclone/issues/1312
  !> Hence this module could be removed once the PSyclone ticket is
  !> completed
  subroutine invoke_spectral_2_cs_kernel_type(fp_spt, longitude, pnm_star, height_wth, &
    stph_spectral_coeffc, stph_spectral_coeffs, &
    spt_level_bottom, spt_level_top, &
    spt_n_max, spt_spectral_dim)

  use spectral_2_cs_kernel_mod, ONLY: spectral_2_cs_code
  use mesh_mod, ONLY: mesh_type

  implicit none

  integer(KIND=i_def), intent(in) :: spt_level_bottom, spt_level_top, spt_n_max, spt_spectral_dim
  type(field_type), intent(in) :: fp_spt, longitude, pnm_star, height_wth
  integer(KIND=i_def) cell
  integer(KIND=i_def) nlayers
  type(field_proxy_type) fp_spt_proxy, longitude_proxy, pnm_star_proxy, height_wth_proxy
  integer(KIND=i_def), pointer :: map_adspc1_longitude(:,:) => null(), map_adspc2_pnm_star(:,:) => null(), &
  &map_wtheta(:,:) => null()
  integer(KIND=i_def) ndf_wtheta, undf_wtheta, ndf_adspc1_longitude, undf_adspc1_longitude, ndf_adspc2_pnm_star, &
  &undf_adspc2_pnm_star
  type(mesh_type), pointer :: mesh => null()

  ! Add arrays stph_spectral_coeffc & stph_spectral_coeffs by hand
  real(kind=r_def), intent(in), dimension(:) :: stph_spectral_coeffc
  real(kind=r_def), intent(in), dimension(:) :: stph_spectral_coeffs
  integer(kind=i_def), parameter :: nranks_array = rank(stph_spectral_coeffc)
  integer(kind=i_def), dimension(nranks_array) :: dims_array

  ! Get the upper bound for each rank of each scalar array
  ! Do dims_array for coeffc and coeffs?
  dims_array = shape(stph_spectral_coeffc)

  !
  ! Initialise field and/or operator proxies
  !
  fp_spt_proxy     = fp_spt%get_proxy()
  longitude_proxy  = longitude%get_proxy()
  pnm_star_proxy   = pnm_star%get_proxy()
  height_wth_proxy = height_wth%get_proxy()
  !
  ! Initialise number of layers
  !
  nlayers = fp_spt_proxy%vspace%get_nlayers()
  !
  ! Create a mesh object
  !
  mesh => fp_spt_proxy%vspace%get_mesh()
  !
  ! Look-up dofmaps for each function space
  !
  map_wtheta           => fp_spt_proxy%vspace%get_whole_dofmap()
  map_adspc1_longitude => longitude_proxy%vspace%get_whole_dofmap()
  map_adspc2_pnm_star  => pnm_star_proxy%vspace%get_whole_dofmap()
  !
  ! Initialise number of DoFs for wtheta
  !
  ndf_wtheta  = fp_spt_proxy%vspace%get_ndf()
  undf_wtheta = fp_spt_proxy%vspace%get_undf()
  !
  ! Initialise number of DoFs for adspc1_longitude
  !
  ndf_adspc1_longitude  = longitude_proxy%vspace%get_ndf()
  undf_adspc1_longitude = longitude_proxy%vspace%get_undf()
  !
  ! Initialise number of DoFs for adspc2_pnm_star
  !
  ndf_adspc2_pnm_star  = pnm_star_proxy%vspace%get_ndf()
  undf_adspc2_pnm_star = pnm_star_proxy%vspace%get_undf()
  !
  ! Call kernels and communication routines
  !
  do cell=1,mesh%get_last_edge_cell()
      call spectral_2_cs_code(nlayers, &
                              ! Add fields
                              fp_spt_proxy%data, longitude_proxy%data, pnm_star_proxy%data, height_wth_proxy%data, &
                              ! Add arrays
                              nranks_array, dims_array, stph_spectral_coeffc, stph_spectral_coeffs, &
                              ! Add SPT scalars
                              spt_level_bottom, spt_level_top, spt_n_max, spt_spectral_dim, &
                              ! Add fields' assoc. space variables
                              ndf_wtheta, undf_wtheta, map_wtheta(:,cell), &
                              ndf_adspc1_longitude, undf_adspc1_longitude, map_adspc1_longitude(:,cell), &
                              ndf_adspc2_pnm_star, undf_adspc2_pnm_star, map_adspc2_pnm_star(:,cell))
  end do
  !
  ! Set halos dirty/clean for fields modified in the above loop
  !
  call fp_spt_proxy%set_dirty()
  !
  !
  end subroutine invoke_spectral_2_cs_kernel_type

end module psykal_lite_phys_mod
