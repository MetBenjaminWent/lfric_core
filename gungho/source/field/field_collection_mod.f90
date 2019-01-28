!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Holds and manages fields in a collection
!>
!> @details A container that holds a collection of fields. Fields that are
!>          presented to the field_collection through the add_field() method are
!>          copied, so when the original goes out of scope, the copy in the
!>          field_collection will continue to be maintained.
!
module field_collection_mod

  use constants_mod,      only: i_def, l_def, str_def
  use field_mod,          only: field_parent_type, &
                                field_type, &
                                field_pointer_type
  use log_mod,            only: log_event, log_scratch_space, &
                                LOG_LEVEL_ERROR, LOG_LEVEL_INFO
  use linked_list_mod,    only: linked_list_type, &
                                linked_list_item_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Type that holds a collection of fields in a linked list
  !-----------------------------------------------------------------------------
  type, public :: field_collection_type

    private
    !> The name of the field collection if provided. 
    character(str_def)     :: name = 'unnamed_collection'

    !> A linked list of the fields contained within the collection
    type(linked_list_type) :: field_list

  contains
    procedure, public :: add_field
    procedure, public :: add_reference_to_field
    procedure, public :: get_field
    procedure, public :: get_iterator
    procedure, public :: clear
    final             :: field_collection_destructor
  end type field_collection_type
  !-----------------------------------------------------------------------------

  interface field_collection_type
    module procedure field_collection_constructor
  end interface

  !-----------------------------------------------------------------------------
  ! Type that iterates through a field collection
  !-----------------------------------------------------------------------------
  type, public :: field_collection_iterator_type

    private
    !> A pointer to the field within the collection that will be
    !> the next to be returned 
    type(linked_list_item_type), pointer :: current

  contains
    procedure, public :: next
    procedure, public :: has_next
  end type field_collection_iterator_type

  interface field_collection_iterator_type
    module procedure field_collection_iterator_constructor
  end interface

contains

!> Constructor for a field collection
!> @param [in] name The name given to the collection
function field_collection_constructor(name) result(self)

  implicit none

  character(*), intent(in), optional :: name

  type(field_collection_type) :: self

  self%field_list = linked_list_type()
  if (present(name))self%name = trim(name)

end function field_collection_constructor

!> Constructor for a field collection iterator
!> @param [in] collection The collection to iterate over
function field_collection_iterator_constructor(collection) result(self)

  implicit none

  type(field_collection_type) :: collection
  type(field_collection_iterator_type) :: self

  self%current => collection%field_list%get_head()

  if(.not.associated(self%current))then
    write(log_scratch_space, '(2A)') &
       'Cannot create an iterator on an empty field collection: ', &
        trim(collection%name)
    call log_event( log_scratch_space, LOG_LEVEL_ERROR)
  end if

end function field_collection_iterator_constructor

!> Adds a field to the collection. The field maintained in the collection will
!> either be a copy of the original or a field pointer object containing a 
!> pointer to a field held elsewhere..
!> @param [in] field The field that is to be copied into the collection or a
!>                   field pointer object that is to be stored in the collection
subroutine add_field(self, field)

  implicit none

  class(field_collection_type), intent(inout) :: self
  class(field_parent_type), intent(in) :: field

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type), pointer :: loop => null()

  ! start at the head of the mesh collection linked list
  loop => self%field_list%get_head()

  do
    ! If list is empty or we've got to the end of list and we didn't find the
    ! field, then we can add it and go home
    if ( .not. associated(loop) ) then
      call self%field_list%insert_item( field )
      exit
    end if

    ! otherwise if we already have the field, then exit with error
    select type(listfield => loop%payload)
      type is (field_type)
        select type(infield => field)
          type is (field_type)
            if ( trim(infield%get_name()) == &
                                   trim(listfield%get_name()) ) then
              write(log_scratch_space, '(4A)') &
                 'Field [', trim(infield%get_name()), &
                 '] already exists in field collection: ', trim(self%name)
              call log_event( log_scratch_space, LOG_LEVEL_ERROR)
            end if
        end select
      type is (field_pointer_type)
        select type(infield => field)
          type is (field_pointer_type)
            if ( trim(infield%field_ptr%get_name()) == &
                                   trim(listfield%field_ptr%get_name()) ) then
              write(log_scratch_space, '(4A)') &
                 'Field [', trim(infield%field_ptr%get_name()), &
                 '] already exists in field collection: ', trim(self%name)
              call log_event( log_scratch_space, LOG_LEVEL_ERROR)
            end if
        end select
    end select

    loop => loop%next
  end do

end subroutine add_field

!> Adds a pointer to a field to the collection. The pointer will point to a
!> field held elsewhere. If that field is destroyed - the pointer will become
!> an orphan
!> @param [in] field_ptr A pointer to a field that is to be referenced in the
!>                       collection
! The routine accepts a pointer to a field. It packages it up into a
! field_pointer object and calls the routine to add this to the collection
subroutine add_reference_to_field(self, field_ptr)

  implicit none

  class(field_collection_type), intent(inout) :: self
  type(field_type), pointer, intent(in) :: field_ptr

  type(field_pointer_type) :: field_pointer

  ! Create a field pointer object that just contains a field pointer
  field_pointer = field_pointer_type( field_ptr )

  call self%add_field( field_pointer )

end subroutine add_reference_to_field

!> Access a field from the collection
!> @param [in] field_name The name of the field to be accessed
!> @return field Pointer to the field that is extracted
function get_field(self, field_name) result(field)

  implicit none

  class(field_collection_type), intent(inout) :: self

  character(*), intent(in) :: field_name
  type(field_type), pointer :: field

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type), pointer :: loop => null()

  ! start at the head of the mesh collection linked list
  loop => self%field_list%get_head()

  do
    ! If list is empty or we're at the end of list and we didn't find the
    ! field, fail with an error
    if ( .not. associated(loop) ) then
      write(log_scratch_space, '(4A)') 'No field [', trim(field_name), &
         '] in field collection: ', trim(self%name)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR)
    end if
    ! otherwise search list for the name of field we want

    ! 'cast' to the field_type 
    select type(listfield => loop%payload)
      type is (field_type)
      if ( trim(field_name) == trim(listfield%get_name()) ) then
          field => listfield
          exit
      end if
      type is (field_pointer_type)
      if ( trim(field_name) == trim(listfield%field_ptr%get_name()) ) then
          field => listfield%field_ptr
          exit
      end if
    end select

    loop => loop%next
  end do

end function get_field

!> Returns an iterator on the field collection
function get_iterator(self) result(iterator)

  implicit none

  class(field_collection_type), intent(inout) :: self
  type(field_collection_iterator_type) :: iterator

  iterator=field_collection_iterator_type(self)

end function get_iterator

!> Clears all items from the field collection linked list
subroutine clear(self)

  implicit none

  class(field_collection_type), intent(inout) :: self

  call self%field_list%clear()

  return
end subroutine clear

!> Destructor for the field collection
subroutine field_collection_destructor(self)

  implicit none

  type (field_collection_type), intent(inout) :: self

  call self%clear()

  return
end subroutine field_collection_destructor

!> Returns the next field form the collection
!> @return field Pointer to the field that is next in the collection
function next(self) result (field)

  implicit none

  class(field_collection_iterator_type), intent(inout), target :: self
  type(field_type), pointer :: field

  ! Extract a pointer to the current field in the collection
  select type(listfield => self%current%payload)
    type is (field_type)
      field => listfield
    type is (field_pointer_type)
      field => listfield%field_ptr
  end select
  ! Move the current field pointer onto the next field in the collection
  self%current => self%current%next

end function next

!> Checks if there are any further fields in the collection being iterated over
!> @return next Returns true if there is another field in the collection,
!>              and  false if there isn't.
function has_next(self) result(next)

  implicit none

  class(field_collection_iterator_type), intent(in) :: self
  logical(l_def) :: next

  next = .true.
  if(.not.associated(self%current)) next = .false.

end function has_next

end module field_collection_mod
