!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Processes diagnostics for bl_exp_alg

module bl_exp_diags_mod

  use constants_mod,       only: l_def
  use field_mod,           only: field_type
  use integer_field_mod,   only: integer_field_type
  use io_config_mod,       only: subroutine_timers
  use timer_mod,           only: timer

  implicit none

  private

  ! Logical indicating whether diagnostics are requested
  logical( l_def ) :: gross_prim_prod_flag

  public :: initialise_diags_for_bl_exp
  public :: output_diags_for_bl_exp

contains

  !> @brief Initialise fields for locally-computed diagnostics

  subroutine initialise_diags_for_bl_exp(zh, gross_prim_prod)

    implicit none

    type( field_type ), intent(in)    :: zh
    type( field_type ), intent(inout) :: gross_prim_prod

    if ( subroutine_timers ) call timer("bl_exp_diags")

    gross_prim_prod_flag = .true.
    call zh%copy_field_properties(gross_prim_prod)

    if ( subroutine_timers ) call timer("bl_exp_diags")

  end subroutine initialise_diags_for_bl_exp

  !> @brief Output diagnostics from bl_exp_alg
  subroutine output_diags_for_bl_exp(ntml, cumulus, bl_type_ind, wvar, dsldzm, &
                                     gradrinr, rhokh_bl, tke_bl, dtrdz_tq_bl,  &
                                     rdz_tq_bl, zhsc, level_ent, level_ent_dsc,&
                                     ent_we_lim, ent_t_frac, ent_zrzi,         &
                                     ent_we_lim_dsc, ent_t_frac_dsc,           &
                                     ent_zrzi_dsc,                             &
                                     tile_fraction, z0m_tile, z0m,             &
                                     gross_prim_prod, net_prim_prod, gc_tile,  &
                                     ustar,                                    &
                                     dust_flux)

    implicit none

    ! Prognostic fields to output
    type( field_type ), intent(in)    :: ntml, cumulus, bl_type_ind, wvar, &
                                         dsldzm, gradrinr, rhokh_bl, tke_bl, &
                                         dtrdz_tq_bl, rdz_tq_bl, zhsc, &
                                         ent_we_lim, ent_t_frac, ent_zrzi, &
                                         ent_we_lim_dsc, ent_t_frac_dsc, &
                                         ent_zrzi_dsc
    type(integer_field_type), intent(in) :: level_ent, level_ent_dsc
    type( field_type ), intent(in)    :: tile_fraction, z0m_tile, z0m, &
                                         gross_prim_prod, net_prim_prod, &
                                         gc_tile, ustar
    type( field_type ), intent(in)    :: dust_flux

    if ( subroutine_timers ) call timer("bl_exp_diags")

    ! Prognostic fields
    call ntml%write_field('turbulence__ntml')
    call cumulus%write_field('turbulence__cumulus')
    call bl_type_ind%write_field('turbulence__bl_type_ind')
    call wvar%write_field('turbulence__wvar')
    call dsldzm%write_field('turbulence__dsldzm')
    call gradrinr%write_field('turbulence__gradrinr')
    call rhokh_bl%write_field('turbulence__rhokh') 
    call tke_bl%write_field('turbulence__tke') 
    call dtrdz_tq_bl%write_field('turbulence__dtrdz_tq') 
    call rdz_tq_bl%write_field('turbulence__rdz_tq') 
    call zhsc%write_field('turbulence__zhsc') 
    call level_ent%write_field('turbulence__level_ent') 
    call level_ent_dsc%write_field('turbulence__level_ent_dsc') 
    call ent_we_lim%write_field('turbulence__ent_we_lim') 
    call ent_t_frac%write_field('turbulence__ent_t_frac') 
    call ent_zrzi%write_field('turbulence__ent_zrzi') 
    call ent_we_lim_dsc%write_field('turbulence__ent_we_lim_dsc') 
    call ent_t_frac_dsc%write_field('turbulence__ent_t_frac_dsc') 
    call ent_zrzi_dsc%write_field('turbulence__ent_zrzi_dsc') 

    call tile_fraction%write_field('surface__tile_fraction')
    call z0m_tile%write_field('surface__z0m_tile')
    call z0m%write_field('surface__z0m')
    call net_prim_prod%write_field('surface__net_prim_prod')
    call gc_tile%write_field('surface__gc_tile') 
    call ustar%write_field('surface__ustar') 

    call dust_flux%write_field('aerosol__dust_flux') 

    ! Diagnostics computed in the kernel
    if (gross_prim_prod_flag) &
         call gross_prim_prod%write_field('surface__gross_prim_prod')


    if ( subroutine_timers ) call timer("bl_exp_diags")

  end subroutine output_diags_for_bl_exp
end module bl_exp_diags_mod
