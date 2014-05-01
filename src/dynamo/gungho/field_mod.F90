!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!> @brief A module providing field related classes.
!>
!> @detail Both the full fat field representation and a light weight proxy are
!> defined here.

module field_mod

  use constants_mod,            only : dp
  use function_space_mod,       only : function_space_type
  use gaussian_quadrature_mod,  only : gaussian_quadrature_type, ngp

  implicit none

  private
  public  :: field_data_from_proxy

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  !> Algorithm layer representation of a field.
  !>
  !> This is a proxy to the actual field description with no accessor methods.
  !>
  type, public :: field_type

    private

    type( field_data_type ), pointer :: real_field => null()

  end type field_type

  interface field_type

    module procedure field_constructor

  end interface field_type

  !> Psy layer representation of a field.
  !>
  !> Objects of this type hold all the data of the field ready for kernels to
  !> work on.
  !>
  type, public :: field_data_type
    private
    !> The number of layers 
    integer :: nlayers
    !> The number of unique degrees of freedom
    integer :: undf
    !> Each field has a pointer to the function space on which it lives
    type( function_space_type ), pointer, public :: vspace => null( )
    !> Each field has a pointer to the gaussian quadrature rule which will be
    !! used to integrate over its values
    type( gaussian_quadrature_type ), pointer, public &
                                      :: gaussian_quadrature => null( )
    !> Allocatable array of type real which holds the values of the field
    real(kind=dp), allocatable, public :: data( : )

  contains
    !> Wrapper to <code>get_ncell</code> from the underlying function space
    !! @param[in] self The calling field
    !! @return An integer, the number of cells in a single layer (columns)
    procedure, public :: get_ncell

    !> Accessor function to get the number of layers
    !! @param[in] self the calling field
    !! @return An integer, the number of layers
    procedure, public :: get_nlayers

    !> Create a new light-weight proxy for this object.
    !>
    !> This is used by the algorithm layer to refer to the object without
    !> having access to its methods.
    !>
    !> @return The proxy object.
    !>
    procedure, public :: new_proxy
  end type field_data_type

  interface field_data_type

    module procedure field_data_constructor

  end interface

contains

  !> Obtain the actual field data object from a proxy object.
  !>
  !> This is a non-type-bound procedure (class method) used by the Psy layer.
  !>
  !> @param [in] proxy The proxy object to dereference.
  !> @return The dereferenced data object.
  function field_data_from_proxy( proxy ) result( data )

    implicit none

    type( field_type ), intent( in ) :: proxy
    type( field_data_type ), pointer :: data

    data => proxy%real_field

  end function field_data_from_proxy

  !---------------------------------------------------------------------------
  ! Constructors
  !---------------------------------------------------------------------------
  !> Construct a <code>field_type</code> object.
  !>
  !> @param [in] field_data the data object the new object should reference.
  !> @return The new object.
  !>
  function field_constructor( field_data ) result( new_field )

    implicit none

    type( field_data_type ), pointer, intent( in ) :: field_data

    type( field_type ) :: new_field

    new_field%real_field => field_data

  end function field_constructor

  !> Construct a <code>field_data_type</code> object.
  !>
  !> @param [in] vector_space the function space that the field lives on
  !> @param [in] gq the gaussian quadrature rule
  !> @param [in] num_layers integer number of layers for the field
  !> @return self the field
  !>
  function field_data_constructor( vector_space, &
                                   gq,           &
                                   num_layers) result(self)

    type(function_space_type), target, intent(in) :: vector_space
    type(gaussian_quadrature_type), target, intent(in) :: gq
    integer, intent(in) :: num_layers

    type(field_data_type), target :: self

    self%vspace => vector_space
    self%gaussian_quadrature => gq
    self%nlayers = num_layers
    self%undf = self%vspace%get_undf()

    ! allocate the array in memory
    allocate(self%data(self%undf))

  end function field_data_constructor

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  !> Wrapper to <code>get_ncell</code> from the underlying function space
  !! @param[in] self The calling field
  !! @return An integer, the number of cells in a single layer (columns)
  integer function get_ncell(self)
    class(field_data_type) :: self
    get_ncell=self%vspace%get_ncell()
    return
  end function get_ncell

  !> Accessor function to get the number of layers
  !! @param[in] self the calling field
  !! @return An integer, the number of layers
  integer function get_nlayers(self)
    class(field_data_type) :: self
    get_nlayers=self%nlayers
    return
  end function get_nlayers

  !> Create a new light-weight proxy for this object.
  !>
  !> This is used by the algorithm layer to refer to the object without
  !> having access to its methods.
  !>
  !> @return The proxy object.
  !>
  type( field_type ) function new_proxy( self ) result( proxy )

    implicit none

    class( field_data_type ), target, intent( in ) :: self

    proxy = field_type( self )

  end function new_proxy

end module field_mod
