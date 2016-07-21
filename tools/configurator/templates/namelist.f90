{#- This is the skeleton of the namelist loading module.                   -#}
{#- The Jinja templating library is used to insert the actual code.        -#}
!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
!> Manages the {{listname}} namelist.
!>
module {{listname}}_config_mod

  use constants_mod, only : {{kindlist | join( ', ' )}}
  use log_mod,       only : log_event, log_scratch_space, LOG_LEVEL_ERROR

  implicit none

  private
  public ::
{%- if enumerations %}
{{- ' ' }}
{%-   for enumeration in enumerations.keys() | sort %}
{{-enumeration}}_from_key, key_from_{{enumeration}}{{','}}
{%-   endfor %} &
{{'           '}}
{%- endif -%}
{{' '}}read_{{listname}}_namelist, {{listname}}_is_loadable, {{listname}}_is_loaded

{%- if enumerations %}
{{-'\n'}}
{%-   for enumeration, pairs in enumerations | dictsort %}
{%-     for pair in pairs %}
  integer(i_native), public, parameter :: {{listname}}_{{enumeration}}_{{pair.key}} = {{pair.value}}
{%-     endfor %}
{%-   endfor %}
{%- endif %}

{%- if parameters %}
{{-'\n'}}
{%-   for name, ftype in parameters | dictsort %}
  {{ftype.typex}}({{ftype.kind}}), public, protected :: {{name}}
{%-   endfor %}
{%- endif %}

  logical :: namelist_loaded = .false.

{%- if enumerations %}
{{-'\n'}}
{%-   for enumeration, pairs in enumerations | dictsort %}
  character(str_short), parameter :: {{enumeration}}_key({{pairs | length()}}) &
{%      set indent = '          = [character(len=str_short) :: ' %}
{%-       for pair in pairs %}
{%-         if not loop.first %}, &{{'\n'}}{% endif %}
{{- indent }}'{{ pair.key }}'
{%-       set indent = ' ' * indent | length() %}
{%-     endfor %}]
{%-   endfor %}
{%- endif %}

contains

{%- for enumeration, pairs in enumerations | dictsort %}

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> \param[in] key Enumeration key.
  !>
  integer(i_native) function {{enumeration}}_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    key_index = 1
    do
      if (trim({{enumeration}}_key(key_index)) == trim(key)) then
        {{enumeration}}_from_key = key_index + {{listname}}_{{enumeration}}_{{pairs[0]['key']}} - 1
        return
      else
        key_index = key_index + 1
        if (key_index > ubound({{enumeration}}_key, 1)) then
          write( log_scratch_space, &
                 '("Key ''", A, "'' not recognised for {{listname}} {{enumeration}}")' ) key
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function {{enumeration}}_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> \param[in] value Enumeration value.
  !>
  character(str_short) function key_from_{{enumeration}}( value )

    implicit none

    integer(i_native), intent(in) :: value

    integer(i_native) :: key_index

    key_index = value - {{listname}}_{{enumeration}}_{{pairs[0]['key']}} + 1
    if (key_index < lbound({{enumeration}}_key, 1) &
        .or. key_index > ubound({{enumeration}}_key, 1)) then
      write( log_scratch_space, &
             '("Value ", I0, " is not in {{listname}} {{enumeration}}")' ) value
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    key_from_{{enumeration}} = {{enumeration}}_key( key_index )

  end function key_from_{{enumeration}}
{%- endfor %}

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> \param [in] file_unit Unit number of the file to read from.
  !>
  subroutine read_{{listname}}_namelist( file_unit )
    implicit none
    integer(i_native), intent(in) :: file_unit
    call read_namelist( file_unit
{%- if enumerations -%}
, {{enumerations.keys() | sort | join( ', ' )}}
{%- endif -%}
                      {{' '}})
  end subroutine read_{{listname}}_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit
{%- if enumerations -%}
, {{ enumerations.keys() | sort | decorate( 'dummy_' ) | join( ', ' ) }}
{%- endif -%}
                          {{' '}})

{%- if constants %}

    use constants_mod, only : {{constants | join( ', ' )}}
{%- endif %}

    implicit none

    integer(i_native), intent(in)  :: file_unit
{%- if enumerations %}
{%-   for enumeration, pairs in enumerations.iteritems() %}
    integer(i_native), intent(out) :: dummy_{{enumeration}}
{%-   endfor %}
{{-'\n'}}
{%- for enumeration in enumerations.keys() | sort %}
    character(str_short) :: {{enumeration}}
{%- endfor %}
{%- endif %}

    namelist /{{listname}}/ {{ variables.keys() | sort | join( ', ' ) }}

    integer(i_native) :: condition

    read( file_unit, nml={{listname}}, iostat=condition, iomsg=log_scratch_space )
    if (condition /= 0) then
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

{%- if enumerations %}
{{-'\n'}}
{%-   for enumeration in enumerations.keys() | sort %}
    dummy_{{enumeration}} = {{enumeration}}_from_key( {{enumeration}} )
{%-   endfor %}
{%- endif %}

    namelist_loaded = .true.
{% for name, code in initialisation.iteritems() %}
    {{name}} = {{code[0]}}
{% endfor %}
  end subroutine read_namelist

  !> Can this namelist be loaded?
  !>
  !> \return True if it is possible to load the namelist.
  !>
  function {{listname}}_is_loadable()

    implicit none

    logical :: {{listname}}_is_loadable

    {{listname}}_is_loadable = .not. namelist_loaded

  end function {{listname}}_is_loadable

  !> Has this namelist been loaded?
  !>
  !> \return True if the namelist has been loaded.
  !>
  function {{listname}}_is_loaded()

    implicit none

    logical :: {{listname}}_is_loaded

    {{listname}}_is_loaded = namelist_loaded

  end function {{listname}}_is_loaded

end module {{listname}}_config_mod
