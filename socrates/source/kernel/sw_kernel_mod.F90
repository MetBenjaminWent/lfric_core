!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Interface to Socrates for shortwave fluxes (external illumination)

module sw_kernel_mod

use argument_mod,      only : arg_type,                  &
                              GH_FIELD, GH_INTEGER,      &
                              GH_READ, GH_WRITE,         &
                              GH_READWRITE, CELLS,       &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              ANY_DISCONTINUOUS_SPACE_3
use fs_continuity_mod, only:  W3, Wtheta
use constants_mod,     only : r_def, i_def, radians_to_degrees
use kernel_mod,        only : kernel_type

implicit none

private

public :: sw_kernel_type
public :: sw_code

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, extends(kernel_type) :: sw_kernel_type
  private
  type(arg_type) :: meta_args(38) = (/                             &
    arg_type(GH_FIELD,   GH_WRITE,     Wtheta),                    & ! sw_heating_rate
    arg_type(GH_FIELD,   GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_surf
    arg_type(GH_FIELD,   GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_surf
    arg_type(GH_FIELD,   GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_blue_surf
    arg_type(GH_FIELD,   GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_blue_surf
    arg_type(GH_FIELD,   GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2), & ! sw_up_tile
    arg_type(GH_FIELD,   GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2), & ! sw_up_blue_tile
    arg_type(GH_FIELD,   GH_READWRITE, Wtheta),                    & ! sw_heating_rate_rts
    arg_type(GH_FIELD,   GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_surf_rts
    arg_type(GH_FIELD,   GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_surf_rts
    arg_type(GH_FIELD,   GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_blue_surf_rts
    arg_type(GH_FIELD,   GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_blue_surf_rts
    arg_type(GH_FIELD,   GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! sw_up_tile_rts
    arg_type(GH_FIELD,   GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! sw_up_blue_tile_rts
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! theta
    arg_type(GH_FIELD,   GH_READ,      W3),                        & ! exner
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! exner_in_wth
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! rho_in_wth
    arg_type(GH_FIELD,   GH_READ,      W3),                        & ! height_w3
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! height_wth
    arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! cos_zenith_angle
    arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! lit_fraction
    arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! cos_zenith_angle_rts
    arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! lit_fraction_rts
    arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! stellar_irradiance_rts
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! mv
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! mcl
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! mci
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! area_fraction
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! liquid_fraction
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! ice_fraction
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! ozone
    arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! tile_fraction
    arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_3), & ! tile_sw_direct_albedo
    arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_3), & ! tile_sw_diffuse_albedo
    arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! latitude
    arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! longitude
    arg_type(GH_INTEGER, GH_READ                                )  & ! timestep
    /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: sw_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

! @param[in]    nlayers                 Number of layers
! @param[out]   sw_heating_rate         SW heating rate
! @param[out]   sw_down_surf            SW downward surface flux
! @param[out]   sw_direct_surf          SW unscattered surface flux
! @param[out]   sw_down_blue_surf       SW blue downward surface flux
! @param[out]   sw_direct_blue_surf     SW blue unscattered surface flux
! @param[out]   sw_up_tile              SW upward tiled surface flux
! @param[out]   sw_up_blue_tile         SW blue upward tiled surface flux
! @param[inout] sw_heating_rate_rts     SW heating rate
! @param[inout] sw_down_surf_rts        SW downward surface flux
! @param[inout] sw_direct_surf_rts      SW unscattered surface flux
! @param[inout] sw_down_blue_surf_rts   SW blue downward surface flux
! @param[inout] sw_direct_blue_surf_rts SW blue unscattered surface flux
! @param[inout] sw_up_tile_rts          SW upward tiled surface flux
! @param[inout] sw_up_blue_tile_rts     SW blue upward tiled surface flux
! @param[in]    theta                   Potential temperature field
! @param[in]    exner                   Exner pressure in density space
! @param[in]    exner_in_wth            Exner pressure in wth space
! @param[in]    rho_in_wth              Density in wth space
! @param[in]    height_w3               Height of w3 levels above surface
! @param[in]    height_wth              Height of wth levels above surface
! @param[in]    cos_zenith_angle        Cosine of the stellar zenith angle
! @param[in]    lit_fraction            Lit fraction of the timestep
! @param[in]    cos_zenith_angle_rts    Cosine of the stellar zenith angle
! @param[in]    lit_fraction_rts        Lit fraction of the timestep
! @param[in]    stellar_irradiance_rts  Stellar irradaince at the planet
! @param[in]    mv                      Water vapour field
! @param[in]    mcl                     Cloud liquid field
! @param[in]    mci                     Cloud ice field
! @param[in]    area_fraction           Total cloud area fraction field
! @param[in]    liquid_fraction         Liquid cloud fraction field
! @param[in]    ice_fraction            Ice cloud fraction field
! @param[in]    ozone                   Ozone field
! @param[in]    tile_fraction           Surface tile fractions
! @param[in]    tile_sw_direct_albedo   SW direct tile albedos
! @param[in]    tile_sw_diffuse_albedo  SW diffuse tile albedos
! @param[in]    latitude                Latitude field
! @param[in]    longitude               Longitude field
! @param[in]    timestep                Timestep number
! @param[in]    ndf_wth                 No. DOFs per cell for wth space
! @param[in]    undf_wth                No. unique of DOFs for wth space
! @param[in]    map_wth                 Dofmap for wth space column base cell
! @param[in]    ndf_2d                  No. of DOFs per cell for 2D space
! @param[in]    undf_2d                 No. unique of DOFs for 2D space
! @param[in]    map_2d                  Dofmap for 2D space column base cell
! @param[in]    ndf_tile                Number of DOFs per cell for tiles
! @param[in]    undf_tile               Number of total DOFs for tiles
! @param[in]    map_tile                Dofmap for tile space column base cell
! @param[in]    ndf_w3                  No. of DOFs per cell for w3 space
! @param[in]    undf_w3                 No. unique of DOFs for w3 space
! @param[in]    map_w3                  Dofmap for w3 space column base cell
! @param[in]    ndf_rtile               No. of DOFs per cell for rtile space
! @param[in]    undf_rtile              No. unique of DOFs for rtile space
! @param[in]    map_rtile               Dofmap for rtile space column base cell
subroutine sw_code(nlayers,                          &
                   sw_heating_rate,                  &
                   sw_down_surf,                     &
                   sw_direct_surf,                   &
                   sw_down_blue_surf,                &
                   sw_direct_blue_surf,              &
                   sw_up_tile,                       &
                   sw_up_blue_tile,                  &
                   sw_heating_rate_rts,              &
                   sw_down_surf_rts,                 &
                   sw_direct_surf_rts,               &
                   sw_down_blue_surf_rts,            &
                   sw_direct_blue_surf_rts,          &
                   sw_up_tile_rts,                   &
                   sw_up_blue_tile_rts,              &
                   theta,                            &
                   exner,                            &
                   exner_in_wth,                     &
                   rho_in_wth,                       &
                   height_w3,                        &
                   height_wth,                       &
                   cos_zenith_angle,                 &
                   lit_fraction,                     &
                   cos_zenith_angle_rts,             &
                   lit_fraction_rts,                 &
                   stellar_irradiance_rts,           &
                   mv,                               &
                   mcl,                              &
                   mci,                              &
                   area_fraction,                    &
                   liquid_fraction,                  &
                   ice_fraction,                     &
                   ozone,                            &
                   tile_fraction,                    &
                   tile_sw_direct_albedo,            &
                   tile_sw_diffuse_albedo,           &
                   latitude, longitude,              &
                   timestep,                         &
                   ndf_wth, undf_wth, map_wth,       &
                   ndf_2d, undf_2d, map_2d,          &
                   ndf_tile, undf_tile, map_tile,    &
                   ndf_w3, undf_w3, map_w3,          &
                   ndf_rtile, undf_rtile, map_rtile)

  use well_mixed_gases_config_mod, only: &
    co2_mix_ratio, n2o_mix_ratio, ch4_mix_ratio, o2_mix_ratio
  use radiation_config_mod, only: n_radstep,                  &
    l_planet_grey_surface, planet_albedo, l_rayleigh_sw,      &
    i_cloud_ice_type_sw, i_cloud_liq_type_sw,                 &
    cloud_representation, cloud_overlap, cloud_inhomogeneity, &
    cloud_representation_no_cloud,                            &
    cloud_representation_liquid_and_ice,                      &
    cloud_representation_conv_strat_liq_ice,                  &
    cloud_overlap_maximum_random,                             &
    cloud_overlap_random,                                     &
    cloud_overlap_exponential_random,                         &
    cloud_inhomogeneity_homogeneous,                          &
    cloud_inhomogeneity_scaling,                              &
    cloud_inhomogeneity_mcica,                                &
    cloud_inhomogeneity_cairns,                               &
    cloud_horizontal_rsd, cloud_vertical_decorr
  use set_thermodynamic_mod, only: set_thermodynamic
  use set_cloud_top_mod, only: set_cloud_top
  use jules_control_init_mod, only: n_surf_tile
  use socrates_init_mod, only: n_sw_band
  use socrates_runes, only: runes, ip_source_illuminate,                    &
    ip_cloud_representation_off, ip_cloud_representation_ice_water,         &
    ip_cloud_representation_csiw, ip_overlap_max_random, ip_overlap_random, &
    ip_overlap_exponential_random, ip_inhom_homogeneous, ip_inhom_scaling,  &
    ip_inhom_mcica, ip_inhom_cairns
  use socrates_bones, only: bones
  use xios, only: xios_date, xios_get_current_date, &
    xios_date_get_day_of_year, xios_date_get_second_of_day

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers, timestep
  integer(i_def), intent(in) :: ndf_wth, ndf_w3, ndf_2d
  integer(i_def), intent(in) :: ndf_tile, ndf_rtile
  integer(i_def), intent(in) :: undf_wth, undf_w3, undf_2d
  integer(i_def), intent(in) :: undf_tile, undf_rtile

  integer(i_def), dimension(ndf_wth),   intent(in) :: map_wth
  integer(i_def), dimension(ndf_w3),    intent(in) :: map_w3
  integer(i_def), dimension(ndf_2d),    intent(in) :: map_2d
  integer(i_def), dimension(ndf_tile),  intent(in) :: map_tile
  integer(i_def), dimension(ndf_rtile), intent(in) :: map_rtile

  real(r_def), dimension(undf_wth),  intent(out) :: sw_heating_rate
  real(r_def), dimension(undf_2d),   intent(out) :: sw_down_surf, &
    sw_direct_surf, sw_down_blue_surf, sw_direct_blue_surf
  real(r_def), dimension(undf_tile), intent(out) :: sw_up_tile, sw_up_blue_tile

  real(r_def), dimension(undf_wth),  intent(inout) :: sw_heating_rate_rts
  real(r_def), dimension(undf_2d),   intent(inout) :: sw_down_surf_rts, &
    sw_direct_surf_rts, sw_down_blue_surf_rts, sw_direct_blue_surf_rts
  real(r_def), dimension(undf_tile), intent(inout) :: sw_up_tile_rts, &
    sw_up_blue_tile_rts

  real(r_def), dimension(undf_w3),   intent(in) :: exner, height_w3
  real(r_def), dimension(undf_wth),  intent(in) :: theta, exner_in_wth, &
    rho_in_wth, height_wth, mv, mcl, mci, &
    area_fraction, liquid_fraction, ice_fraction, ozone
  real(r_def), dimension(undf_2d), intent(in) :: &
    cos_zenith_angle, lit_fraction, &
    cos_zenith_angle_rts, lit_fraction_rts, stellar_irradiance_rts
  real(r_def), dimension(undf_tile),  intent(in) :: tile_fraction
  real(r_def), dimension(undf_rtile), intent(in) :: &
    tile_sw_direct_albedo, tile_sw_diffuse_albedo
  real(r_def), dimension(undf_2d), intent(in) :: latitude, longitude

  ! Local variables for the kernel
  type(xios_date) :: datetime
  integer(i_def), parameter :: n_profile = 1
  integer(i_def) :: i_cloud_representation, i_overlap, i_inhom
  integer(i_def) :: rand_seed(n_profile)
  integer(i_def) :: k, n_cloud_layer
  integer(i_def) :: df_rtile, i_tile, i_band
  integer(i_def) :: wth_1, wth_nlayers
  real(r_def), dimension(nlayers) :: &
    ! Heat capacity for each layer
    layer_heat_capacity,             &
    ! Layer pressure and temperature
    p_layer, t_layer,                &
    ! Mass of layer per square metre
    d_mass,                          &
    ! Effective radius of droplets
    liq_dim
  real(r_def), dimension(0:nlayers) :: sw_direct, sw_down, sw_up
  real(r_def), parameter :: dl = 100.0_r_def
    ! Number of unique seeds per degree of latitude / longitude (for MCICA).
    ! 100 per degree equates to approximately a 1km square area of the
    ! globe (for the Earth).

  ! Tiled surface fields
  real(r_def) :: frac_tile(n_profile, n_surf_tile)
  real(r_def) :: albedo_diff_tile(n_profile, n_surf_tile, n_sw_band)
  real(r_def) :: albedo_dir_tile(n_profile, n_surf_tile, n_sw_band)
  real(r_def) :: flux_up_tile_rts(n_surf_tile)
  real(r_def) :: flux_up_blue_tile_rts(n_surf_tile)
  real(r_def) :: flux_up_tile(n_surf_tile)
  real(r_def) :: flux_up_blue_tile(n_surf_tile)


  ! Set wth indexing
  wth_1 = map_wth(1)+1
  wth_nlayers = map_wth(1)+nlayers

  if (mod(timestep-1_i_def, n_radstep) == 0) then
    ! Radiation time-step: full calculation of radiative fluxes

    ! Properties of clouds
    select case (cloud_representation)
    case (cloud_representation_no_cloud)
      i_cloud_representation = ip_cloud_representation_off
    case (cloud_representation_liquid_and_ice)
      i_cloud_representation = ip_cloud_representation_ice_water
    case (cloud_representation_conv_strat_liq_ice)
      i_cloud_representation = ip_cloud_representation_csiw
    case default
      i_cloud_representation = ip_cloud_representation_off
    end select
    select case (cloud_overlap)
    case (cloud_overlap_maximum_random)
      i_overlap = ip_overlap_max_random
    case (cloud_overlap_random)
      i_overlap = ip_overlap_random
    case (cloud_overlap_exponential_random)
      i_overlap = ip_overlap_exponential_random
    case default
      i_overlap = ip_overlap_max_random
    end select
    select case (cloud_inhomogeneity)
    case (cloud_inhomogeneity_homogeneous)
      i_inhom = ip_inhom_homogeneous
    case (cloud_inhomogeneity_scaling)
      i_inhom = ip_inhom_scaling
    case (cloud_inhomogeneity_mcica)
      i_inhom = ip_inhom_mcica
      call xios_get_current_date(datetime)
      ! Generate a unique seed for each gridpoint (actually for each area of
      ! the globe with a size defined by dl) at the given time in seconds.
      rand_seed &
        = int((latitude(map_2d(1))*radians_to_degrees + 90.0_r_def)*dl, i_def) &
        + int(longitude(map_2d(1))*radians_to_degrees*dl, i_def)*180_i_def*dl &
        + ((abs(int(datetime%year, i_def)-2000_i_def)*366_i_def) &
        + int(xios_date_get_day_of_year(datetime), i_def))*86400_i_def &
        + int(xios_date_get_second_of_day(datetime), i_def)
    case (cloud_inhomogeneity_cairns)
      i_inhom = ip_inhom_cairns
    case default
      i_inhom = ip_inhom_homogeneous
    end select

    ! Hardwire droplet effective radius for now
    do k=1, nlayers
      liq_dim(k) = 7.0e-6_r_def
    end do

    ! Use the highest cloud layer for the number of cloud layers
    call set_cloud_top(nlayers, &
      area_fraction(map_wth(1):map_wth(1)+nlayers), n_cloud_layer)

    ! Set up pressures, temperatures, masses and heat capacities
    call set_thermodynamic(nlayers,                &
      exner(map_w3(1):map_w3(1)+nlayers-1),        &
      exner_in_wth(map_wth(1):map_wth(1)+nlayers), &
      theta(map_wth(1):map_wth(1)+nlayers),        &
      rho_in_wth(map_wth(1):map_wth(1)+nlayers),   &
      height_w3(map_w3(1):map_w3(1)+nlayers-1),    &
      height_wth(map_wth(1):map_wth(1)+nlayers),   &
      p_layer, t_layer, d_mass, layer_heat_capacity)

    ! Tile fractions
    do i_tile = 1, n_surf_tile
      frac_tile(1, i_tile) = tile_fraction(map_tile(i_tile))
    end do

    ! Tile albedos
    df_rtile = 0
    do i_band = 1, n_sw_band
      do i_tile = 1, n_surf_tile
        df_rtile = df_rtile + 1
        albedo_diff_tile(1, i_tile, i_band) &
          = tile_sw_diffuse_albedo(map_rtile(df_rtile))
        albedo_dir_tile(1, i_tile, i_band) &
          = tile_sw_direct_albedo(map_rtile(df_rtile))
      end do
    end do

    ! Calculate the SW fluxes (RUN the Edwards-Slingo two-stream solver)
    call runes(n_profile, nlayers,                                             &
      spectrum_name          = 'sw',                                           &
      i_source               = ip_source_illuminate,                           &
      n_cloud_layer          = n_cloud_layer,                                  &
      p_layer_1d             = p_layer,                                        &
      t_layer_1d             = t_layer,                                        &
      mass_1d                = d_mass,                                         &
      density_1d             = rho_in_wth(wth_1:wth_nlayers),                  &
      h2o_1d                 = mv(wth_1:wth_nlayers),                          &
      o3_1d                  = ozone(wth_1:wth_nlayers),                       &
      co2_mix_ratio          = co2_mix_ratio,                                  &
      n2o_mix_ratio          = n2o_mix_ratio,                                  &
      ch4_mix_ratio          = ch4_mix_ratio,                                  &
      o2_mix_ratio           = o2_mix_ratio,                                   &
      cos_zenith_angle       = cos_zenith_angle_rts(map_2d(1):map_2d(1)),      &
      solar_irrad            = stellar_irradiance_rts(map_2d(1):map_2d(1)),    &
      l_grey_albedo          = l_planet_grey_surface,                          &
      grey_albedo            = planet_albedo,                                  &
      l_tile                 = .not.l_planet_grey_surface,                     &
      n_tile                 = n_surf_tile,                                    &
      frac_tile              = frac_tile,                                      &
      albedo_diff_tile       = albedo_diff_tile,                               &
      albedo_dir_tile        = albedo_dir_tile,                                &
      cloud_frac_1d          = area_fraction(wth_1:wth_nlayers),               &
      liq_frac_1d            = liquid_fraction(wth_1:wth_nlayers),             &
      ice_frac_1d            = ice_fraction(wth_1:wth_nlayers),                &
      liq_mmr_1d             = mcl(wth_1:wth_nlayers),                         &
      ice_mmr_1d             = mci(wth_1:wth_nlayers),                         &
      liq_dim_1d             = liq_dim,                                        &
      cloud_horizontal_rsd   = cloud_horizontal_rsd,                           &
      cloud_vertical_decorr  = cloud_vertical_decorr,                          &
      rand_seed              = rand_seed,                                      &
      layer_heat_capacity_1d = layer_heat_capacity,                            &
      l_rayleigh             = l_rayleigh_sw,                                  &
      l_mixing_ratio         = .true.,                                         &
      i_cloud_representation = i_cloud_representation,                         &
      i_overlap              = i_overlap,                                      &
      i_inhom                = i_inhom,                                        &
      i_st_water             = i_cloud_liq_type_sw,                            &
      i_st_ice               = i_cloud_ice_type_sw,                            &
      l_invert               = .true.,                                         &
      flux_direct_1d         = sw_direct,                                      &
      flux_down_1d           = sw_down,                                        &
      flux_up_1d             = sw_up,                                          &
      flux_up_tile_1d        = flux_up_tile_rts,                               &
      flux_up_blue_tile_1d   = flux_up_blue_tile_rts,                          &
      flux_direct_blue_surf  = sw_direct_blue_surf_rts(map_2d(1):map_2d(1)),   &
      flux_down_blue_surf    = sw_down_blue_surf_rts(map_2d(1):map_2d(1)),     &
      heating_rate_1d        = sw_heating_rate_rts(wth_1:wth_nlayers))

    ! Set level 0 increment such that theta increment will equal level 1
    sw_heating_rate_rts(map_wth(1)) = sw_heating_rate_rts(map_wth(1) + 1) &
                                    * exner_in_wth(map_wth(1))            &
                                    / exner_in_wth(map_wth(1) + 1)

    ! Set surface fluxes
    sw_down_surf_rts(map_2d(1)) = sw_down(0)
    sw_direct_surf_rts(map_2d(1)) = sw_direct(0)

    ! Set tiled fluxes
    do i_tile = 1, n_surf_tile
      sw_up_tile_rts(map_tile(i_tile)) = flux_up_tile_rts(i_tile)
      sw_up_blue_tile_rts(map_tile(i_tile)) = flux_up_blue_tile_rts(i_tile)
    end do
  else
    ! Not a radiation time-step

    ! Set tiled fluxes
    do i_tile = 1, n_surf_tile
      flux_up_tile_rts(i_tile) = sw_up_tile_rts(map_tile(i_tile))
      flux_up_blue_tile_rts(i_tile) = sw_up_blue_tile_rts(map_tile(i_tile))
    end do
  end if

  if (n_radstep == 1) then
    ! Radiation timestep = model timestep
    sw_heating_rate(map_wth(1):wth_nlayers) &
      = sw_heating_rate_rts(map_wth(1):wth_nlayers)

    sw_down_surf(map_2d(1))        = sw_down_surf_rts(map_2d(1))
    sw_direct_surf(map_2d(1))      = sw_direct_surf_rts(map_2d(1))
    sw_down_blue_surf(map_2d(1))   = sw_down_blue_surf_rts(map_2d(1))
    sw_direct_blue_surf(map_2d(1)) = sw_direct_blue_surf_rts(map_2d(1))

    do i_tile = 1, n_surf_tile
      sw_up_tile(map_tile(i_tile))      = sw_up_tile_rts(map_tile(i_tile))
      sw_up_blue_tile(map_tile(i_tile)) = sw_up_blue_tile_rts(map_tile(i_tile))
    end do
  else
    ! Corrections to model timestep. The radiative fluxes have been calculated
    ! for a mean sun angle over the radiation timestep and must be converted
    ! to fluxes appropriate for the mean sun angle over the shorter model
    ! timestep.

    ! The bare "bones" of a simple radiative transfer calculation.
    ! Apply the solar zenith angle correction.
    call bones(n_profile, nlayers,                                             &
      n_tile                    = n_surf_tile,                                 &
      l_cos_zen_correction      = .true.,                                      &
      cos_zen_rts               = cos_zenith_angle_rts(map_2d(1):map_2d(1)),   &
      lit_frac_rts              = lit_fraction_rts(map_2d(1):map_2d(1)),       &
      cos_zen_mts               = cos_zenith_angle(map_2d(1):map_2d(1)),       &
      lit_frac_mts              = lit_fraction(map_2d(1):map_2d(1)),           &
      heating_rate_1d_rts       = sw_heating_rate_rts(wth_1:wth_nlayers),      &
      flux_up_tile_1d_rts       = flux_up_tile_rts,                            &
      flux_up_blue_tile_1d_rts  = flux_up_blue_tile_rts,                       &
      flux_direct_surf_rts      = sw_direct_surf_rts(map_2d(1):map_2d(1)),     &
      flux_down_surf_rts        = sw_down_surf_rts(map_2d(1):map_2d(1)),       &
      flux_direct_blue_surf_rts = sw_direct_blue_surf_rts(map_2d(1):map_2d(1)),&
      flux_down_blue_surf_rts   = sw_down_blue_surf_rts(map_2d(1):map_2d(1)),  &
      heating_rate_1d_mts       = sw_heating_rate(wth_1:wth_nlayers),          &
      flux_up_tile_1d_mts       = flux_up_tile,                                &
      flux_up_blue_tile_1d_mts  = flux_up_blue_tile,                           &
      flux_direct_surf_mts      = sw_direct_surf(map_2d(1):map_2d(1)),         &
      flux_down_surf_mts        = sw_down_surf(map_2d(1):map_2d(1)),           &
      flux_direct_blue_surf_mts = sw_direct_blue_surf(map_2d(1):map_2d(1)),    &
      flux_down_blue_surf_mts   = sw_down_blue_surf(map_2d(1):map_2d(1)))

    ! Set level 0 increment such that theta increment will equal level 1
    sw_heating_rate(map_wth(1)) = sw_heating_rate(map_wth(1) + 1) &
                                * exner_in_wth(map_wth(1))        &
                                / exner_in_wth(map_wth(1) + 1)

    ! Set tiled fluxes
    do i_tile = 1, n_surf_tile
      sw_up_tile(map_tile(i_tile)) = flux_up_tile(i_tile)
      sw_up_blue_tile(map_tile(i_tile)) = flux_up_blue_tile(i_tile)
    end do
  end if

end subroutine sw_code

end module sw_kernel_mod
