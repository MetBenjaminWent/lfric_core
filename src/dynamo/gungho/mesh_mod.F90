!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Defines various mesh details.

!> @details information about the mesh is held in  here.
module mesh_mod
implicit none

integer :: num_cells
integer :: num_layers
integer :: element_order
logical :: l_spherical

integer :: v_unique_dofs(4,2)
integer :: v_dof_entity(4,0:3)    

end module mesh_mod

