! #define DEBUG
module f_util
  implicit none

  public select_triplet
  public periodic_distance
  public match_structure
  public tr_cl_sc

contains

  ! translate cluster to supercell
  subroutine tr_cl_sc(nAsc, nCell, nA, npt, scijkl, sc_mat, inv_sc_mat, sc_ref, apos, cluster, tr_cl)
    implicit none
    integer, intent(in) :: nAsc, nCell, nA, npt
    integer(8), intent(in) :: scijkl(4, nAsc), sc_ref(3, nCell), cluster(4, npt)
    real(8), intent(in) :: sc_mat(3,3), inv_sc_mat(3,3), apos(3, nA)
    integer(8), intent(out) :: tr_cl(npt, nCell)
    !f2py intent(in) scijkl, sc_mat(3,3), inv_sc_mat, sc_ref, apos, cluster
    !f2py intent(out) tr_cl
    !f2py integer, intent(hide),depend(scijkl) :: nAsc=shape(scijkl,1)
    !f2py integer, intent(hide),depend(sc_ref) :: nCell=shape(sc_ref,1)
    !f2py integer, intent(hide),depend(apos) :: nA=shape(apos,1)
    !f2py integer, intent(hide),depend(cluster) :: npt=shape(cluster,1)

    integer i, j, iA, found
    integer ct(4)
    real(8) ct_f(3)
!    if (npt >8) then
!       print *, 'ERROR recompile with ct(4,maxnpt) maxnpt>=', npt
!       stop
!    end if

 !   ct(4,1:npt) = cluster(4,1:npt)
    do i=1, nCell
       !  [get_structure_ordering((coords+ijk[None,:]).dot(sc.inv_sc_mat), sc.frac_coords,0) for ijk in sc.ijk_ref]
       do j=1, npt
          ct(4) = cluster(4, j)
          ct(1:3) = cluster(1:3, j) + sc_ref(1:3, i)
          ct_f = matmul(inv_sc_mat, ct(1:3)+ apos(:,ct(4)+1))+1e-12
          ct_f = ct_f - floor(ct_f)
          ct_f = matmul(sc_mat, ct_f)- apos(:,ct(4)+1)
          ct(1:3) = nint(ct_f(1:3))
          found=-1
          do iA=1, nAsc
!             print *, 'i,j, iA;  clus;  scref;  ct;  scijkl'
!             print '(4I4)', i,j,iA, found, cluster(:, j), sc_ref(:,i), -99, ct, scijkl(:, iA)
             if (scijkl(1, iA) /= ct(1)) cycle
             if (scijkl(2, iA) /= ct(2)) cycle
             if (scijkl(3, iA) /= ct(3)) cycle
             if (scijkl(4, iA) /= ct(4)) cycle
             found = iA
!             print *, 'debug found', found
             exit
          end do
          if (found <0) then
             print *, 'ERROR translating cluster', cluster
             stop
          end if
          tr_cl(j, i) = found
       end do
    end do
  end subroutine tr_cl_sc

  
  ! latt(:,1-3) the three lattice vectors
  ! p1 is the "perpect structure"
  ! p2 is the structure to be analyzed
  ! n_i: number of atoms
  ! tp_i: type of atoms
  ! pi:   for each atom in p2, its matching (position only) atom index in p1. -1 for no matching
  ! stat: for each atom in p2, its matching status:
  !     M=matched, m=matched, tol_m<d<=tol_v
  !     S=substituted (matched but different type)
  !     s=substituted (matched, tol_m<d<=tol_v but different type)
  !     I=interstitial (nothing within tol_v)
  ! dist:  for each atom in p2, the distance to the matched atom in p1
  ! pi1:   for each atom in p1, its matching (position only) atom index in p2. -1 for no matching

  function periodic_distance(latt, dr12)
    implicit none
    real(8) periodic_distance
    integer, parameter  :: dim=3
    real(8), intent(in) :: latt(dim, 3), dr12(dim)

    integer i
    real(8) dr_all(dim, 27), dr(dim)
    real(8), dimension(dim, 27) :: img0
    img0= reshape((/-1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,0   ,  1,  1,  1,&
      & 1, 1, 1, 1, 1, 1, -1, -1, -1, 0, 0, 0,  1,1   ,  1, -1, -1, -1,  0,  0,  0,  1,  1,  1, -1, -1,&
      & -1,  0,  0,  0,1   ,  1,  1, -1,  0,  1, -1,  0,  1, -1,  0,  1, -1,  0,  1, -1,  0,1   , -1,  0,&
      & 1, -1, 0, 1, -1, 0,  1, -1,  0,  1/), shape(img0))

    dr = dr12
    dr = dr - nint(dr)
    do i=1,27
       dr_all(:,i)=dr+img0(:,i)
    end do
    periodic_distance = minval(norm2(matmul(latt, dr_all),1))
  end function periodic_distance


  subroutine match_structure(latt, n1, p1, tp1, n2, p2, tp2, tol_m, tol_v, st1, pi1, di1, st2, pi2, di2)
    implicit none
    integer, parameter  :: dim=3
    real(8), intent(in) :: latt(3, 3)
    integer, intent(in) :: n1, n2
    real(8), intent(in) :: p1(3, n1), p2(3, n2), tol_m, tol_v
    integer, intent(in) :: tp1(n1), tp2(n2)
    character, intent(out):: st1(n1), st2(n2)
    integer, intent(out):: pi1(n1), pi2(n2)
    real(8), intent(out):: di1(n1), di2(n2)
    !f2py intent(in) latt, p1, tp1, p2, tp2, tol_m, tol_v
    !f2py intent(out) st1, pi1, di1, st2, pi2, di2
    !f2py integer, intent(hide), depend(p1) :: n1=shape(p1, 1)
    !f2py integer, intent(hide), depend(p2) :: n2=shape(p2, 1)

    real(8) dr, min_dr, tmp
    real(8), allocatable :: d_table(:, :)
    integer i, j, i_min, j_min, iloop

    allocate(d_table(n2, n1))
!    !$ OMP PARALLEL DO PRIVATE(fcrd)
    do i=1, n1
       do j=1, n2
          d_table(j,i)= periodic_distance(latt, p1(:,i)- p2(:,j))
       end do
    end do
!    !$ OMP END PARALLEL DO

    pi1=-1
    pi2=-1
    di1=1e20
    di2=1e20
    st1='V'
    st2='I'
    do iloop=1, min(n1, n2)
       min_dr = 1E20
       do i=1, n1
          if (pi1(i)>=0) cycle
          do j=1, n2
             if (pi2(j)>=0) cycle
             if (d_table(j,i)<min_dr) then
                i_min = i
                j_min = j
                min_dr = d_table(j,i)
             end if
          end do
       end do
       if (min_dr<= tol_v) then
          pi1(i_min) = j_min-1
          pi2(j_min) = i_min-1
          di1(i_min) = min_dr
          di2(j_min) = min_dr
          if (tp1(i_min)==tp2(j_min)) then
             st1(i_min)='M'
             st2(j_min)='M'
          else
             st1(i_min)='S'
             st2(j_min)='S'
          end if
          if (min_dr> tol_m) then
             st1(i_min)=char(ichar(st1(i_min))+32)
             st2(j_min)=char(ichar(st2(j_min))+32)
          end if
       else
          exit
       end if
    end do
    deallocate(d_table)
  end subroutine match_structure


  subroutine select_triplet(poscar, rmat, cutoff, nlat, natom)
    implicit none
    integer, intent(in) :: nlat, natom
    double precision, intent(in) ::  cutoff
    double precision, intent(in) ::  poscar(natom,3)
    double precision, intent(in) ::  rmat(3,3)
    !f2py intent(in) :: nlat
    !f2py intent(in) :: poscar, rmat, cutoff 
    !f2py integer, intent(hide),depend(poscar) :: natom=shape(poscar,0)

!!!f2py intent(in, out) :: dataout

    ! local var
    integer i, j, k,  ntotal, nlatset
    integer iatom1, iatom2, iatom3
    integer ilat1, ilat2, ilat3
    integer tmp1, tmp2
    double precision radius
    double precision clus(3,4)
    double precision pos1(3), pos2(3), pos3(3)
    integer lat1(3), lat2(3), lat3(3)

    !out
    integer, allocatable :: dataout(:,:)

    !inside
    double precision, allocatable :: latset(:,:)


    !generate whole lattice site
    nlatset=(2*nlat+1)*(2*nlat+1)*(2*nlat+1)
    allocate(latset(nlatset,3))
    tmp1=0
    do i=-nlat, nlat
       do j=-nlat, nlat
          do k=-nlat, nlat
             tmp1=tmp1+1
             latset(tmp1,1)=i
             latset(tmp1,2)=j
             latset(tmp1,3)=k
          end do
       end do
    end do
    print *, "Total length of nlatset: ", tmp1
    do i=1, nlatset
       print '(3I5)', int(latset(i,1)),int(latset(i,2)),int(latset(i,3))
    end do

    !loop over lattice site and atomic site
    ntotal=0d0
    do iatom1=1, natom
       do ilat2=1, nlatset
          do iatom2=1, natom
             do ilat3=1, nlatset
                do iatom3=1, natom
                   call getrad(poscar(iatom1,:), latset(ilat2,:), poscar(iatom2,:), &
                        latset(ilat3,:), poscar(iatom3,:), rmat, radius)
                   if (radius <= cutoff) then
                      ntotal=ntotal+1
                      print *, "radius: ", radius
                      print 100, iatom1, ilat2, iatom2, ilat3, iatom3
                      print *
                   end if
                end do
             end do
          end do
       end do
    end do
    print *, "Total selected clusters: ", ntotal
100 format("lat1: ", i8,  "   lat2: ", 2i8, "   lat3: ", 2i8)

    !write to output file
    allocate(dataout(ntotal,12))
    ntotal=0d0
    dataout=0d0
    do iatom1=1, natom
       do ilat2=1, nlatset
          do iatom2=1, natom
             do ilat3=1, nlatset
                do iatom3=1, natom
                   call getrad(poscar(iatom1,:), latset(ilat2,:), poscar(iatom2,:), &
                        latset(ilat3,:), poscar(iatom3,:), rmat, radius)
                   if (radius <= cutoff) then
                      ntotal=ntotal+1
                      dataout(ntotal, 4)=iatom1-1d0
                      dataout(ntotal, 5:7)=latset(ilat2,:)
                      dataout(ntotal, 8)=iatom2-1d0
                      dataout(ntotal, 9:11)=latset(ilat3,:)
                      dataout(ntotal, 12)=iatom3-1d0
                      print "(3(4I5,';'))", dataout(ntotal,:)
                   end if
                end do
             end do
          end do
       end do
    end do
    print *, "Total selected clusters: ", ntotal

    !write to file  
    open(33,file='triplet-selection',form='formatted',status='unknown')
    do i=1, ntotal
       write(33, "(12I5)") dataout(i,:)
    end do
    close(33)
    return
  end subroutine select_triplet


  subroutine getrad(pos1, lat2, pos2, lat3, pos3, rmat, radius)
    double precision pos1(3), pos2(3), pos3(3)
    double precision lat1(3), lat2(3), lat3(3)
    double precision rads(3), tmp(3)
    double precision rmat(3,3)  
    double precision clus(3,3)

    double precision radius, maxrad
    integer i, j, k

    lat1=0d0
    clus(1,:)=lat1(:)+pos1(:)
    clus(2,:)=lat2(:)+pos2(:)
    clus(3,:)=lat3(:)+pos3(:)

    !transform to direct coordinate
    clus=matmul(clus, rmat)

    !first radius
    tmp=clus(1,:)-clus(2,:)
    rads(1)=(dot_product(tmp,tmp))**(0.5)
    maxrad=rads(1)
    !second radius
    tmp=clus(1,:)-clus(3,:)
    rads(2)=(dot_product(tmp,tmp))**(0.5)
    if (rads(2)>maxrad) then
       maxrad=rads(2)
    end if
    !third radius
    tmp=clus(2,:)-clus(3,:)
    rads(3)=(dot_product(tmp,tmp))**(0.5)
    if (rads(3)>maxrad) then
       maxrad=rads(3)
    end if

    radius=maxrad
  end subroutine getrad


  subroutine writemtx(Alist, totNF, Ncorr)
    integer eff, i, j
    integer totNF, Ncorr
    double precision Alist(totNF, Ncorr)

    eff=0d0
    do i=1, totNF
       do j=1, Ncorr
          if (abs(Alist(i,j))>1d-10) then
             eff=eff+1
          end if
       end do
    end do

    open(40,file='Amatrix.mtx',form='formatted',status='unknown')
    write(40, *) "%%MatrixMarket matrix coordinate real general"
    write(40, *) "%Created by Wolfram Mathematica 9.0 : www.wolfram.com"
    print *, "To be written: "
    print *, totNF, Ncorr, eff
    write(40, "(3I16)") totNF, Ncorr, eff
    do i=1, totNF
       do j=1, Ncorr
          if (abs(Alist(i,j))>1d-10) then
             write(40, "(2I5, E16.8)") i, j, Alist(i,j)
          end if
       end do
    end do

    close(40)

    return

  end subroutine writemtx
end module f_util
