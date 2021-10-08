!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Initialise, define and finalise the linearisation state.

module linear_model_data_mod

  use constants_mod,                  only : i_def, r_def, str_def
  use clock_mod,                      only : clock_type
  use gungho_model_data_mod,          only : model_data_type
  use pure_abstract_field_mod,        only : pure_abstract_field_type
  use field_mod,                      only : field_type
  use field_collection_mod,           only : field_collection_type
  use finite_element_config_mod,      only : element_order
  use init_gungho_prognostics_alg_mod,only : init_gungho_prognostics_alg
  use mr_indices_mod,                 only : nummr, &
                                             mr_names
  use moist_dyn_mod,                  only : num_moist_factors
  use moist_dyn_factors_alg_mod,      only : moist_dyn_factors_alg
  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W2, W3, WTheta
  use linear_config_mod,              only : ls_option,         &
                                             ls_option_analytical
  use linear_data_algorithm_mod,      only : linear_copy_model_to_ls, &
                                             linear_init_pert_random
  use log_mod,                        only : log_event,         &
                                             log_scratch_space, &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_ERROR

  implicit none

  public linear_create_ls,        &
         linear_setup_ls_field,   &
         linear_init_ls,          &
         linear_init_pert

contains

  !> @brief   Create the fields in the ls fields field collection.
  !> @details At present, this only includes the preparation for an
  !!          analytical definition. But this could be extended to include
  !!          the preparation for reading ls fields from a file,
  !!          with a time axis.
  !> @param[inout] model_data The working data set for a model run
  !> @param[in]    mesh_id The identifier given to the current 3d mesh
  !> @param[in]    twod_mesh_id The identifier given to the current 2d mesh
  subroutine linear_create_ls( model_data, mesh_id, &
                               twod_mesh_id )

    implicit none

    type( model_data_type ), target, intent(inout) :: model_data
    integer(i_def),                  intent(in)    :: mesh_id
    integer(i_def),                  intent(in)    :: twod_mesh_id

    type( field_collection_type ), pointer :: depository => null()
    type( field_collection_type ), pointer :: ls_fields => null()
    type( field_type ),            pointer :: ls_mr(:) => null()
    type( field_type ),            pointer :: ls_moist_dyn(:) => null()

    integer(i_def)     :: imr
    character(str_def) :: name
    character(str_def) :: moist_dyn_name

    depository => model_data%depository
    ls_fields => model_data%ls_fields
    ls_mr => model_data%ls_mr
    ls_moist_dyn => model_data%ls_moist_dyn

    write(log_scratch_space,'(A,A)') "Create ls fields: "// &
          "Setting up ls field collection"
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    ls_fields = field_collection_type(name='ls_fields')

    if ( ls_option == ls_option_analytical ) then

      call linear_setup_ls_field( &
           "ls_rho", depository, ls_fields, W3, mesh_id )
      call linear_setup_ls_field( &
           "ls_exner", depository, ls_fields, W3, mesh_id )
      call linear_setup_ls_field( &
           "ls_theta", depository, ls_fields, Wtheta, mesh_id )
      call linear_setup_ls_field( &
           "ls_u", depository, ls_fields, W2, mesh_id )

      do imr = 1,nummr
        name = trim('ls_' // adjustl(mr_names(imr)) )
        call linear_setup_ls_field( &
             name, depository, ls_fields, Wtheta, mesh_id, mr=ls_mr, imr=imr )
      end do

      do imr = 1, num_moist_factors
        write(moist_dyn_name, "(A12, I1)") "ls_moist_dyn", imr
        name = trim(moist_dyn_name)
        call linear_setup_ls_field( &
             name, depository, ls_fields, Wtheta, mesh_id, mr=ls_moist_dyn, &
             imr=imr )
      end do

    else
      call log_event( "LS setup not available for requested ls_option ", &
           LOG_LEVEL_ERROR)
    end if

  end subroutine linear_create_ls

  !> @brief Create a linearisation state field and add to the depository.
  !> @param[in]    name        The field name
  !> @param[inout] depository  The depository field collection
  !> @param[inout] ls_fields   The linearisation state field collection
  !> @param[in]    fs          The function space of the field
  !> @param[in]    mesh_id     The identifier given to the current 3d mesh
  !> @param[inout] mr          An array of fields that hold the linearisation
  !!                           state moisture mixing ratios
  !> @param[in]    imr         The index of the moisture mixing ratio array
  subroutine linear_setup_ls_field( name, depository, ls_fields, fs, mesh_id, mr, imr )

    implicit none

    character(*),                       intent(in)     :: name
    type( field_collection_type ),      intent(inout)  :: depository
    type( field_collection_type ),      intent(inout)  :: ls_fields
    integer (i_def),                    intent(in)     :: fs
    integer(i_def),                     intent(in)     :: mesh_id
    type(field_type), optional,         intent(inout)  :: mr(:)
    integer(i_def), optional,           intent(in)     :: imr

    type(function_space_type),       pointer :: field_space => null()
    class(pure_abstract_field_type), pointer :: abs_fld_ptr => null()
    type(field_type)                         :: new_field

    write(log_scratch_space,'(3A,I6)') &
         "Creating new field for ", trim(name)
    call log_event(log_scratch_space,LOG_LEVEL_INFO)

    field_space => function_space_collection%get_fs( &
                                        mesh_id, element_order, fs )

    ! Initialise the field and add to the depository.
    if ( present(imr) ) then
      call mr(imr)%initialise( field_space, name=trim(name) )
      call depository%add_field(mr(imr))
    else
      call new_field%initialise( field_space, name=trim(name) )
      call depository%add_field(new_field)
    end if

    ! Add the field pointer to the target field collection
    abs_fld_ptr => depository%get_field(name)
    call ls_fields%add_reference_to_field(abs_fld_ptr)

    ! Nullify pointers
    abs_fld_ptr => null()
    field_space => null()

  end subroutine linear_setup_ls_field

  !> @brief   Define the linearisation state values.
  !> @details At the present, these can only be defined from an analytical
  !!          solution.
  !> @param[in]    mesh_id      The identifier given to the current 3d mesh
  !> @param[in]    twod_mesh_id The identifier given to the current 2d mesh
  !> @param[inout] model_data   The working data set for a model run
  !> @param[in]    clock        The model time
  subroutine linear_init_ls( mesh_id, twod_mesh_id, model_data, clock )

    use gungho_step_mod,                only : gungho_step

    implicit none

    integer(i_def),                  intent(in)    :: mesh_id
    integer(i_def),                  intent(in)    :: twod_mesh_id
    type( model_data_type ), target, intent(inout) :: model_data
    class(clock_type),               intent(in)    :: clock

    integer :: i
    type( field_type ), pointer                    :: ls_field => null()

    if(ls_option == ls_option_analytical) then

      ! Procedure to define the linearisation state from an analytical field:
      ! 1. Define the analytical field in the gungho prognostics. (This
      !    is done in initialise_model_data in tl_test_driver).
      ! 2. Evolve the prognostic fields using gungho_step -
      !    this avoids the linearisation state being zero.
      ! 3. Copy the prognostic fields to the linearisation fields, and set
      !    the prognostic fields to zero.

      ! Evolve the prognostic fields.
      do i=1,10
        call gungho_step( mesh_id,      &
                          twod_mesh_id, &
                          model_data,   &
                          clock )
      end do

      ! Copy the prognostic fields to the LS and then zero the prognostics.
      call linear_copy_model_to_ls( model_data )

    else
     call log_event('This ls_option not available', LOG_LEVEL_ERROR)
    end if

    ! Print the min and max values of the linearisation fields.
    ls_field => model_data%ls_fields%get_field("ls_u")
    call ls_field%log_minmax(LOG_LEVEL_INFO,'ls_u')

    ls_field => model_data%ls_fields%get_field("ls_rho")
    call ls_field%log_minmax(LOG_LEVEL_INFO,'ls_rho')

    ls_field => model_data%ls_fields%get_field("ls_exner")
    call ls_field%log_minmax(LOG_LEVEL_INFO,'ls_exner')

    ls_field => model_data%ls_fields%get_field("ls_theta")
    call ls_field%log_minmax(LOG_LEVEL_INFO,'ls_theta')

  end subroutine linear_init_ls

  !> @brief   Define the initial perturbation values.
  !> @details Define the initial perturbation - currently from random data
  !> @param[in]    mesh_id      The identifier given to the current 3d mesh
  !> @param[in]    twod_mesh_id The identifier given to the current 2d mesh
  !> @param[inout] model_data   The working data set for a model run
  subroutine linear_init_pert( mesh_id, twod_mesh_id, model_data )

    implicit none

    integer(i_def),                  intent(in)    :: mesh_id
    integer(i_def),                  intent(in)    :: twod_mesh_id
    type( model_data_type ), target, intent(inout) :: model_data

    call linear_init_pert_random( mesh_id,      &
                                  twod_mesh_id, &
                                  model_data )

  end subroutine linear_init_pert

end module linear_model_data_mod
