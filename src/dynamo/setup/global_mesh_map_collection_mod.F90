!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
!
!>
!> @brief   Holds and manages the multiple global mesh maps
!>
!> @details A container which holds a collection of global mesh maps
!>          It will handle the creation and storing of requested global maps
!>          from a single source global mesh to a number of appropriately
!>          specified different global meshes.
!
module global_mesh_map_collection_mod

  use constants_mod,       only: i_def
  use linked_list_mod,     only: linked_list_type, linked_list_item_type
  use global_mesh_map_mod, only: global_mesh_map_type
  use log_mod,             only: log_event, log_scratch_space                  &
                               , LOG_LEVEL_TRACE, LOG_LEVEL_ERROR

  implicit none

  private

  type, public :: global_mesh_map_collection_type
    private

    !> Linked list of global_mesh_map_type objects.
    type(linked_list_type) :: global_mesh_map_list

    !> The id of the global_mesh object that is the source for all
    !> maps in this collection.
    integer(i_def)         :: source_id

    !> The number of cells in the global_mesh object that is the
    !> source for all maps in this collection.
    integer(i_def)         :: source_ncells

  contains

    !> @brief  Returns the id of the source global mesh object
    procedure, public :: get_source_id

    !> @brief  Returns the number of cells in the source global mesh object
    procedure, public :: get_source_ncells

    !> @brief     Adds a global_mesh_map object to the collection
    !> @param[in] target_id         Integer id of the global mesh object which
    !>                              this maps the source global mesh object to.
    !> @param[in] global_mesh_map   Integer array [::] which provides global
    !>                              cell ids of the target global mesh object
    !>                              which overlaps with source global mesh
    !>                              cells.
    !>                              This connectivity should be in the
    !>                              dimensions of [number of target cells for
    !>                              each source cell, number of source cells]
    procedure, public :: add_global_mesh_map

    !> @brief     Returns pointer to global_mesh_map object which maps global
    !>            cell ids in the target global mesh to global cell ids in the
    !>            source global mesh.
    !> @param[in] target_id         Integer id of the target global mesh object
    !>                              of requested global_mesh_map object.
    procedure, public :: get_global_mesh_map

    procedure, public :: clear
    final             :: global_mesh_map_collection_destructor

  end type global_mesh_map_collection_type

  interface global_mesh_map_collection_type
    module procedure global_mesh_map_collection_constructor
  end interface

contains

!==============================================================================
!> @brief  Constructs the mesh map collection object
!> @return self The constructed global_mesh_map_collection object
function global_mesh_map_collection_constructor(source_id, source_ncells)      &
     result(self)

  implicit none

  type(global_mesh_map_collection_type) :: self
  integer(i_def) :: source_id
  integer(i_def) :: source_ncells

  self%global_mesh_map_list = linked_list_type()
  self%source_id     = source_id
  self%source_ncells = source_ncells

end function global_mesh_map_collection_constructor


function get_source_id(self) result (source_id)

  implicit none

  class(global_mesh_map_collection_type), intent(in) :: self
  integer(i_def) :: source_id

  source_id = self%source_id

end function get_source_id

function get_source_ncells(self) result (source_ncells)

  implicit none

  class(global_mesh_map_collection_type), intent(in) :: self
  integer(i_def) :: source_ncells

  source_ncells = self%source_ncells

end function get_source_ncells

!==============================================================================

subroutine add_global_mesh_map(self, target_id, intermesh_gid_map)

  implicit none

  class(global_mesh_map_collection_type), intent(inout) :: self
  integer(i_def), intent(in) :: target_id
  integer(i_def), intent(in) :: intermesh_gid_map(:,:)

  type(global_mesh_map_type) :: global_mesh_map

  integer(i_def) :: global_mesh_map_id

  global_mesh_map_id = 1000*self%source_id + target_id

  if (self%global_mesh_map_list%item_exists(global_mesh_map_id)) then
    ! Do nothing as map already exists
    write(log_scratch_space, '(A,I0,A)')                                       &
        'Skipping task: Global mesh map (id: ', global_mesh_map_id,            &
        ') already exists.'
    call log_event(log_scratch_space, LOG_LEVEL_TRACE)
    return
  end if

  ! The intermesh_gid_map is dimensioned:
  ! ( n mapped cells from target per source cell, ncells in source )
  ! Check that the intermesh_gid_map is valid for this source 
  if (size(intermesh_gid_map,2) /= self%source_ncells) then
    ! Do nothing as map already exists
    write(log_scratch_space, '(A,I0,A,I0,A)')                                  &
        'Invalid intermesh GID mapping: Number of source cells in GID map (',  &
        size(intermesh_gid_map,2),                                             &
        ') does not match number of source cells for this map collection (',   &
        self%source_ncells,')'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    return
  end if

  global_mesh_map = global_mesh_map_type( self%source_id, target_id            &
                                        , intermesh_gid_map )

  call self%global_mesh_map_list%insert_item( global_mesh_map )

  return
end subroutine add_global_mesh_map

!==============================================================================

function get_global_mesh_map( self, target_id ) result (global_mesh_map)

  implicit none

  class(global_mesh_map_collection_type), intent(in) :: self
  integer(i_def), intent(in) :: target_id

  type(global_mesh_map_type), pointer :: global_mesh_map
  integer(i_def) :: global_mesh_map_id

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type),pointer :: loop => null()

  global_mesh_map_id = 1000*self%source_id + target_id

  if (self%global_mesh_map_list%item_exists(global_mesh_map_id)) then

    loop => self%global_mesh_map_list%get_head()

    do
      ! Loop over list looking of correct item id
      if ( global_mesh_map_id == loop%payload%get_id() ) then

        select type(m => loop%payload)
        type is (global_mesh_map_type)
          global_mesh_map => m
        end select
        exit

      end if

      loop => loop%next
    end do

  else

    ! Requested map does not exist between this collections source
    ! global_mesh and the target global mesh.
    nullify(global_mesh_map)
    write(log_scratch_space, '(A,I0,A)')                                       &
        'Requested global mesh map (id: ', global_mesh_map_id,                 &
        ') does not exist.'
    call log_event(log_scratch_space, LOG_LEVEL_TRACE)
    return

  end if

end function get_global_mesh_map

!==============================================================================

!> Clear all items from the linked list of global_mesh_maps
subroutine clear(self)

  implicit none

  class(global_mesh_map_collection_type), intent(inout) :: self

  call self%global_mesh_map_list%clear()

  return
end subroutine clear

!==============================================================================

! Global mesh map collection destructor
subroutine global_mesh_map_collection_destructor(self)

  implicit none

  type (global_mesh_map_collection_type), intent(inout) :: self

  call self%clear()

  return
end subroutine global_mesh_map_collection_destructor

end module global_mesh_map_collection_mod
