!!<summary>The bcs algorithm requires matrix operations such as appending, deleting
!!finding the intersections and complements of data sets. Since these
!!methods are pretty general, it makes sense to put them in their own module.</summary>
!!<comments date="27 Jun 2013" author="Conrad W. Rosenbrock">
!!- Changed array reallocations to use move_alloc instead of normal copy/overwrite.</comments>
module matrix_sets
  use num_types
  implicit none
  public intersection, complement, appendto, delete

  !!<summary>Appends an element to the end if a list. Accepts either real or
  !!integer lists and values.</summary>
  INTERFACE appendto
     module procedure appendto_one_d_list_real, appendto_one_d_list_integer
  END INTERFACE appendto
  !!<summary>Removes an element at a specified index from the list. Real or
  !!integer lists supported</summary>
  INTERFACE delete
     module procedure delete_element_from_one_d_list_real,delete_element_from_one_d_list_integer
  END INTERFACE delete
contains
  !!<summary>Determines the intersection of the specified lists. Returns the values
  !!that are common to both lists in the intersect parameter. Comparison is made using
  !!the built-in fortran operator @CREF[.eq.]. Zero values are not considered common even
  !!if they are in both lists.</summary>
  !!<parameter name="listone, listtwo">The two lists to compare values for and take the
  !!intersection of.</parameter>
  !!<parameter name="intersect">The resultant list of indices. Will be overwritten if
  !!already allocated.</parameter>
  SUBROUTINE intersection(listone, listtwo, intersect)
    integer, intent(in) :: listone(:), listtwo(:)
    integer, allocatable, intent(inout) :: intersect(:)

    !!<local name="N,M">The size of the lists 1 and 2 respectively.</local>
    !!<local name="indices">Counter for the number of indices common to both lists.</local>
    !!<local name="numberofnonzero">The final number of non-zero values common to both lists.</local>
    integer, allocatable ::templist(:)
    integer i, N, M, indices, numberofnonzero
    
    N = size(listone, 1)
    M = size(listtwo, 1)
  
    if (N > M) then
       allocate(templist(N))
    else
       allocate(templist(M))
    endif

    templist = 0
    indices = 1

    do i=1, N
       if (any(listtwo .eq. listone(i))) then
          templist(indices) = listone(i)
          indices = indices + 1
       endif
    enddo

    !Some of the elements in listone could have been zero and thus common to both lists.
    !But we are only interested in the non-zero values.
    numberofnonzero = count(templist > 0)
    if (allocated(intersect)) deallocate(intersect)
    allocate(intersect(numberofnonzero))
    intersect(1:numberofnonzero) = templist(1:numberofnonzero)
  END SUBROUTINE intersection

  !!<summary>Determines the complement of the two lists. Returns the values
  !!of those elements in the first list that are not in the second list.</summary>
  !!<parameter name="listone, listtwo" regular="true">The two lists of integer values to
  !!take the complement of.</parameter>
  !!<parameter name="comp">The resulting list of values that are only in list one.</parameter>
  SUBROUTINE complement(listone,listtwo,comp)
    integer, intent(in) :: listone(:), listtwo(:)
    integer, allocatable, intent(inout) :: comp(:)
    
    !!<local name="N,M">The size of the lists 1 and 2 respectively.</local>
    !!<local name="indices">Counter for the number of values exclusive to list one.</local>
    !!<local name="numberofnonzero">The final number of non-zero values unique to list one.</local>
    integer, allocatable ::templist(:)
    integer i, N, M, indices,numberofnonzero

    N =size(listone,1)
    M = size(listtwo,1)

    if (N > M) then
       allocate(templist(N))
    else
       allocate(templist(M))
    endif

    templist = 0
    indices = 1
    do i=1, N
       if (any(listtwo .eq. listone(i))) then
          cycle
       else
          templist(indices) = listone(i)
          indices = indices + 1
       endif
    enddo

    numberofnonzero = count(templist > 0)
    if (allocated(comp)) deallocate(comp)
    allocate(comp(numberofnonzero))
    comp(1:numberofnonzero) = templist(1:numberofnonzero)
  END SUBROUTINE complement

  !!<summary>Removes the element at the specified index from the real-valued list.</summary>
  !!<parameter name="list">The list of real values to delete from.</parameter>
  !!<parameter name="deleteindex">The index whose value should be removed from the list.</parameter>
  SUBROUTINE delete_element_from_one_d_list_real(list,deleteindex)
    real(dp), allocatable, intent(inout) :: list(:)
    integer, intent(in) :: deleteindex
    real(dp), allocatable :: templist(:)

    allocate(templist(size(list,1)-1))
    if (deleteindex .gt. 1) then
       templist(1:deleteindex-1) = list(1:deleteindex-1)
       templist(deleteindex:size(list,1)-1) = list(deleteindex+1:size(list, 1))
    else
       templist(:) = list(2:size(list,1))
    end if    
    call move_alloc(templist, list)
  END SUBROUTINE delete_element_from_one_d_list_real

  !!<summary>Removes the element at the specified index from the integer-valued list.</summary>
  !!<parameter name="list">The list of real values to delete from.</parameter>
  !!<parameter name="deleteindex">The index whose value should be removed from the list.</parameter>
  SUBROUTINE delete_element_from_one_d_list_integer(list, deleteindex)
    integer, allocatable, intent(inout) :: list(:)
    integer, intent(in) :: deleteindex
    integer, allocatable :: templist(:)

    allocate(templist(size(list,1)-1))
    if (deleteindex .gt. 1) then
       templist(1:deleteindex-1) = list(1:deleteindex-1)
       templist(deleteindex:size(list,1)-1) = list(deleteindex+1:size(list, 1))
    else
       templist(:) = list(2:size(list,1))
    end if    
    call move_alloc(templist, list)
  END SUBROUTINE delete_element_from_one_d_list_integer

  !!<summary>Adds the specified value to the end of the real-valued list.</summary>
  !!<parameter name="list">The real-valued list to append to.</parameter>
  !!<parameter name="appendage">The extra value to append to the list.</parameter>
  SUBROUTINE appendto_one_d_list_real(list,appendage)
    real(dp), allocatable, intent(inout) :: list(:)
    real(dp), intent(in) :: appendage
    real(dp), allocatable :: templist(:)

    allocate(templist(size(list,1)+1))
    templist(1:size(list,1)) = list
    templist(size(list,1)+1) = appendage

    call move_alloc(templist, list)
  END SUBROUTINE appendto_one_d_list_real

  !!<summary>Adds the specified value to the end of the integer-valued list.</summary>
  !!<parameter name="list">The integer-valued list to append to.</parameter>
  !!<parameter name="appendage">The extra value to append to the list.</parameter>
  SUBROUTINE appendto_one_d_list_integer(list, appendage)
    integer, allocatable, intent(inout) :: list(:)
    integer, intent(in) :: appendage
    integer, allocatable :: templist(:)

    allocate(templist(size(list,1)+1))
    templist(1:size(list,1)) = list
    templist(size(list,1)+1) = appendage

    call move_alloc(templist, list)
  END SUBROUTINE appendto_one_d_list_integer
end module matrix_sets
