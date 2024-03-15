!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>  @brief Support for configuring input/output files
module lfric_xios_metafile_mod

  use constants_mod,                 only: str_def, i_def, l_def
  use log_mod,                       only: log_event,                         &
                                           log_level_error,                   &
                                           log_level_info,                    &
                                           log_level_warning
  use xios,                          only: xios_file,                         &
                                           xios_field,                        &
                                           xios_get_handle,                   &
                                           xios_add_child,                    &
                                           xios_set_attr,                     &
                                           xios_is_defined_field_attr,        &
                                           xios_get_field_attr,               &
                                           xios_set_field_attr,               &
                                           xios_is_defined_fieldgroup_attr,   &
                                           xios_get_fieldgroup_attr
  use lfric_xios_diag_mod,           only: field_is_valid,                    &
                                           get_field_grid_ref,                &
                                           get_field_domain_ref,              &
                                           get_field_axis_ref
  implicit none

  !> @brief Wrap XIOS file handle
  type :: metafile_type
    type(xios_file) :: handle
  contains
    private
    procedure, public :: init => metafile_init
    procedure, public :: get_handle => metafile_get_handle

    ! destructor - here to avoid gnu compiler bug
    final :: metafile_destructor
  end type metafile_type

private

public :: metafile_type, add_field

contains

  !> @brief Support for legacy checkpoints
  !> @details This affects lbc and gungho prognostics, which have traditionally been
  !> using the legacy checkpoint domains checkpoint_W3, checkpoint_W2 etc.
  !> @param[in] field    XIOS field object
  !> @param[in] field_id XIOS fielld id
  subroutine handle_legacy_fields(field, field_id)
    implicit none

    type(xios_field), intent(in) :: field
    character(*), intent(in) :: field_id

    character(20), parameter :: Wtheta = 'checkpoint_Wtheta'
    character(20), parameter :: W3 = 'checkpoint_W3'
    character(20), parameter :: W2 = 'checkpoint_W2'

    character(20), parameter :: wtheta_fields(7) = [character(20) :: &
      'theta', 'm_v', 'm_cl', 'm_r', 'm_ci', 'm_s', 'm_g']
    character(20), parameter :: w3_fields(3) =  [character(20) ::  &
     'rho', 'exner', 'ageofair']
    character(20), parameter :: w2_fields(1) = [character(20) :: &
      'u']
    character(20), parameter :: lbc_wtheta_fields(7) = [character(20) ::  &
      'lbc_theta', 'lbc_m_v', 'lbc_m_cl', 'lbc_m_r', 'lbc_m_ci', 'lbc_m_s', 'lbc_m_g']
    character(20), parameter :: lbc_w3_fields(2) = [character(20) :: &
      'lbc_rho', 'lbc_exner']
    character(20), parameter :: lbc_w2_fields(3) = [character(20) :: &
      'lbc_u', 'boundary_u_diff', 'boundary_u_driving']

    character(20) :: domain_id

    if (any(wtheta_fields == field_id)) then
      domain_id = Wtheta
    else if (any(w3_fields == field_id)) then
      domain_id = W3
    else if (any(w2_fields == field_id)) then
      domain_id = W2
    else if (any(lbc_wtheta_fields == field_id)) then
      domain_id = Wtheta
    else if (any(lbc_w3_fields == field_id)) then
      domain_id = W3
    else if (any(lbc_w2_fields == field_id)) then
      domain_id = W2
    else
      domain_id = ''
      call log_event('unexpected legacy field: ' // trim(field_id), log_level_error)
    end if

    call xios_set_attr(field, domain_ref=domain_id)
  end subroutine handle_legacy_fields

  !> @brief Get file handle from XIOS
  !> @param[inout] self  Metafile object
  !> @param[in] file_id  XIOS id of file
  subroutine metafile_init(self, file_id)
    implicit none
    class(metafile_type), intent(in out) :: self

    character(*), intent(in) :: file_id

    call xios_get_handle(file_id, self%handle)
  end subroutine metafile_init

  !> @brief Accessor for handle
  !> @param[inout] self  Metafile object
  !> @param[out] handle  XIOS file handle
  subroutine metafile_get_handle(self, handle)
    implicit none
    class(metafile_type), intent(in) :: self
    type(xios_file), intent(out) :: handle

    handle = self%handle
  end subroutine metafile_get_handle

  !> @brief Destructor of metafile object
  !> @param[inout] self  Metafile object
  subroutine metafile_destructor(self)
    implicit none
    type(metafile_type), intent(inout) :: self
    ! empty
  end subroutine metafile_destructor

  !> @brief Get field precision from XIOS
  !> @param[in] field_id  XIOS id of field
  !> @param[in] dflt      Default precision
  !> @param[in] group_id  XIOS id of fieldgroup
  !> @return              Precision returned
  function get_field_precision(field_id, dflt, group_id) result(prec)
    implicit none
    character(*), intent(in) :: field_id
    integer(i_def), intent(in) :: dflt
    character(*), optional, intent(in) :: group_id

    integer(i_def) :: prec
    logical(l_def) :: has_prec

    if (field_is_valid(field_id)) then
      call xios_is_defined_field_attr(field_id, prec=has_prec)
      if (has_prec) then
        call xios_get_field_attr(field_id, prec=prec)
        return
      end if
    end if

    if (present(group_id)) then
      call xios_is_defined_fieldgroup_attr(group_id, prec=has_prec)
      if (has_prec) then
        call xios_get_fieldgroup_attr(group_id, prec=prec)
        return
      end if
    end if

    prec = dflt
  end function get_field_precision

  !> @brief Add copy of dictionary field to file
  !> @param[in] metafile       XIOS file wrapper
  !> @param[in] dict_field_id  XIOS id of dictionary field to be copied
  !> @param[in] prefix         ID prefix to be used, .e.g, "checkpoint_"
  !> @param[in] operation      XIOS field operation, e.g., "once"
  !> @param[in] id_as_name     Use dictionary field ID as field name?
  !> @param[in] legacy         Use legacy checkpointing domain?
  subroutine add_field(metafile, dict_field_id, prefix, operation, id_as_name, legacy)
    implicit none
    type(metafile_type), intent(in) :: metafile
    character(*), intent(in) :: dict_field_id
    character(*), intent(in) :: prefix
    character(*), intent(in) :: operation
    logical(l_def), optional, intent(in) :: id_as_name
    logical(l_def), optional, intent(in) :: legacy

    character(20), parameter :: lfric_dict = 'lfric_dictionary'
    integer(i_def), parameter :: dflt_prec = 8

    type(xios_file)    :: file
    character(str_def) :: field_id
    character(str_def) :: field_name
    type(xios_field)   :: field
    character(str_def) :: grid_ref
    character(str_def) :: domain_ref
    character(str_def) :: axis_ref
    integer(i_def)     :: prec
    logical(l_def)     :: use_id_as_name
    logical(l_def)     :: use_legacy

    use_id_as_name = .false.
    if (present(id_as_name)) use_id_as_name = id_as_name

    use_legacy = .false.
    if (present(legacy)) use_legacy = legacy

    call metafile%get_handle(file)

    field_id = trim(prefix) // dict_field_id

    if (field_is_valid(field_id)) then
      ! old style checkpoint configuration - enable field
      call xios_set_field_attr(field_id, enabled=.true.)
    else
      ! new style - add field to checkpoint file
      call xios_add_child(file, field, field_id)
      call log_event('adding checkpoint field ' // trim(field_id), log_level_info)
      if (.not. field_is_valid(field_id)) &
        call log_event('internal error: added field invalid', log_level_error)

      ! copy name and precision from dictionary field
       if (use_id_as_name) then
        field_name = dict_field_id ! correct for checkpointing
      else
        call xios_get_field_attr(dict_field_id, name=field_name)
        if (field_name /= dict_field_id) &
          call log_event('internal error - mismatch: ' // trim(field_name) &
            // ' vs ' // trim(dict_field_id), log_level_warning)
      end if

      prec = get_field_precision(dict_field_id, dflt_prec, lfric_dict)

      call xios_set_attr(field, name=field_name, prec=prec, operation=operation)

      if (use_legacy) then
        call handle_legacy_fields(field, dict_field_id)
      else
        grid_ref = get_field_grid_ref(dict_field_id)
        if (grid_ref /= '') then
          call xios_set_attr(field, grid_ref=grid_ref)
        else
          domain_ref = get_field_domain_ref(dict_field_id)
          axis_ref = get_field_axis_ref(dict_field_id)
          if (domain_ref /= '') call xios_set_attr(field, domain_ref=domain_ref)
          if (axis_ref /= '') call xios_set_attr(field, axis_ref=axis_ref)
        end if
      end if
    end if

  end subroutine add_field

end module lfric_xios_metafile_mod
