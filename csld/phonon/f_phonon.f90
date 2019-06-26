! #define DEBUG
module f_phonon
  implicit none

  public init
  public init_nac
  public init_dpcor
  public get_dispersion
  public get_eig_e_vec
  public get_dos_new
  public calc_thermal
  public get_fcm_dipole
  public get_dm

  real(8), private, parameter :: PI      = 3.14159265358979323846
  real(8), private, parameter :: UNITTHZ = 15.6333046
  real(8), private, parameter :: BohrInAngstrom = 0.52917721067
  real(8), private, parameter :: HartreeIneV =  27.21138602
  !  real, public, parameter :: THztoEV =  0.0041356668
  !  real, public, parameter :: THztoinverCM =  33.3564227583
  real(8), private, parameter :: kBbyeV = 8.61733238E-5
  real(8), private, parameter :: atomic2eVAng =  HartreeIneV*BohrInAngstrom


  !
  ! NOTE: all real space coordinates (incl. apos, disp) are cartesian
  ! Therefore reciprocal space coordinates should also be be cartesian
  !
  complex(kind=8), allocatable, private  :: dm(:,:,:)
  real(8), private :: avec(3,3), qvec(3,3)
  integer, private :: natoms, nfcm, dim_dm
  real(8), allocatable, private  :: mass(:)
  real(8), allocatable, private :: apos(:,:)
  real(8), allocatable, private :: disp(:,:)
  real(8), allocatable, private :: fcm(:,:,:)
  integer, allocatable, private :: typ(:,:)
  ! non-analytic correction parameters
  integer, private :: nac_method = -1, Ncell_mixsp
  real(8), private :: diel(3,3), rho2, pcVol, pref_nac, etafac=2.0
  integer, allocatable, private :: idx_mixsp(:)
  real(8), allocatable, private :: Zstar(:,:,:), Zscalar(:)
  ! force constants for long range ASR corrections
  integer, private :: nfcm_c
  real(8), allocatable, private :: disp_c(:,:)
  real(8), allocatable, private :: fcm_c(:,:,:)
  integer, allocatable, private :: typ_c(:,:)

  ! for matrix diagonalization
  integer, private :: lw
  real(8), private, allocatable :: RWORK(:), W(:)
  complex(kind=8), private, allocatable :: WORK(:)
  complex(kind=8), allocatable, private  :: dmtmp(:,:)

contains


!!!
  !! Initialize lattice
  !!
  subroutine init(a_lat, ma, x_f, na, ijk, tp, fc, nf)
    real(8), intent(in) :: a_lat(3,3)
    integer, intent(in) :: na, nf
    real(8), intent(in) :: ma(na)
    real(8), intent(in) :: x_f(na, 3)
    integer, intent(in) :: ijk(nf, 3)
    integer, intent(in) :: tp(nf, 2)
    real(8), intent(in) :: fc(nf, 3, 3)
    !f2py intent(in) a_lat, ma, x_f, ijk, tp, fc
    !f2py integer, intent(hide),depend(ma) :: na=shape(ma,0)  
    !f2py integer, intent(hide),depend(tp) :: nf=shape(tp,0)  

    ! local variables

    natoms=na
    dim_dm = 3*natoms
    lw= MAX(1,2*dim_dm-1)
    nfcm=nf
    if (allocated(mass)) deallocate(mass)
    if (allocated(apos)) deallocate(apos)
    if (allocated(disp)) deallocate(disp)
    if (allocated(typ)) deallocate(typ)
    if (allocated(fcm)) deallocate(fcm)
    allocate(mass(natoms), apos(3,natoms), disp(nfcm, 3), fcm(nfcm, 3, 3), typ(nfcm, 2))


    avec= transpose(a_lat)
    ! NOTE: DMdipole always has column-first
    ! so the lattice vectors are avec(:,1), avec(:,2), avec(:,3)
    call lattc(avec, qvec)
    mass(:)= ma(:)
    apos= transpose(matmul(x_f, a_lat))
    fcm(:,:,:) = fc(:,:,:)
    ! index starting at 1, not 0 in python
    typ(:,:)= tp(:,:) + 1
    !
    ! NOTE: we adopt the phase convention that drops the explicit atomic coordinates in the density matrix
    !
    !    do i=1, nfcm
    !       disp(i, :) = x(typ(i, 1), :) - x(typ(i, 2), :)
    !    end do
    disp = matmul(ijk, a_lat)

    !    print *, 'debug natoms, nfcm', natoms, nfcm, na, nf
    !    print *, 'debug mass ', mass
    !    print *, 'debug avec ', avec
    !    print *, 'debug x ', x
    !    print *, 'debug disp ', disp
    !    print *, 'debug typ ', typ
    !    print *, 'debug fcm ', fcm

    if (allocated(WORK)) deallocate(WORK)
    if (allocated(RWORK)) deallocate(RWORK)
    if (allocated(W)) deallocate(W)
    allocate(WORK(lw), RWORK(MAX(1,3*dim_dm-2)), W(dim_dm))
    if (allocated(dmtmp)) deallocate(dmtmp)
    allocate(dmtmp(3*natoms,3*natoms))
  end subroutine init


  ! Init force constants for long range ASR corrections
  subroutine init_dpcor(ijk, tp, fc, nf)
    integer, intent(in) :: nf
    integer, intent(in) :: ijk(nf, 3)
    integer, intent(in) :: tp(nf, 2)
    real(8), intent(in) :: fc(nf, 3, 3)
    !f2py intent(in) ijk, tp, fc
    !f2py integer, intent(hide),depend(tp) :: nf=shape(tp,0)  

    nfcm_c = nf
    if (allocated(disp_c)) deallocate(disp_c)
    if (allocated(fcm_c)) deallocate(fcm_c)
    if (allocated(typ_c)) deallocate(typ_c)
    allocate(disp_c(nfcm_c, 3), fcm_c(nfcm_c, 3, 3), typ_c(nfcm_c, 2))
    fcm_c(:,:,:) = fc(:,:,:)
    typ_c(:,:)= tp(:,:) + 1
    disp_c = matmul(ijk, transpose(avec))
  end subroutine init_dpcor

  !
  ! Initialize dielectric tensor and Born effective charge
  !
  subroutine init_nac(method, d_tensor, born_chg, na, r_i, v_i, ijk_idx, nijk, etafac_in)
    integer, intent(in) :: method
    real(8), intent(in) :: d_tensor(3,3), r_i, v_i
    integer, intent(in) :: na, nijk, ijk_idx(nijk)
    real(8), intent(in) :: born_chg(na, 3, 3)
    real(8), intent(in) :: etafac_in
    !f2py intent(in) method, d_tensor, born_chg, r_i, v_i, ijk_idx
    !f2py integer, intent(hide),depend(born_chg) :: na=shape(born_chg,0)
    !f2py integer, intent(hide),depend(ijk_idx) :: nijk=shape(ijk_idx,0)
    integer i

    ! should be 0, 1 or 2
    nac_method = method
    if ((nac_method < -1) .or. (nac_method > 2)) then
       print *, 'expecting NAC method -1,0,1,2  but found', nac_method
    end if
    if (natoms /= na) then
       print *, 'ERROR num of atoms in born chg ', na, ' != ', natoms
    end if
    if (allocated(Zstar)) deallocate(Zstar)
    if (allocated(Zscalar)) deallocate(Zscalar)
    if (allocated(idx_mixsp)) deallocate(idx_mixsp)
    Ncell_mixsp = nijk
    allocate(Zstar(3, 3, natoms), idx_mixsp(Ncell_mixsp), Zscalar(natoms))

    diel= d_tensor
    do i=1, natoms
       Zstar(:,:,i) = born_chg(i,:,:)
    end do
    Zscalar(:) = Zstar(1,1,:)
    idx_mixsp(:) = ijk_idx(:) +1
    rho2 = r_i*r_i
    !       print *, 'debug v_i pcVol ', v_i, pcVol, ' r_i', r_i
    pcVol = v_i
    pref_nac = 4*PI*atomic2eVAng/pcVol
    etafac = etafac_in

    !    print *, 'debug method ', method
    !    print *, 'debug diel ', d_tensor
    !    print *, 'debug born1 ', born_chg(1,:,:)
    !    print *, 'debug born2 ', born_chg(2,:,:)
    !    print *, 'debug rho ', r_i
    !    print *, 'debug pcVol ', pcVol
  end subroutine init_NAC


  ! return DM to external call
  subroutine get_dm(kpts, nk, nb, dm_out)
    integer, intent(in)   :: nk, nb
    real(8), intent(in)   :: kpts(nk, 3)
    complex(kind=8),intent(out)   :: dm_out(nb, nb, nk)
    !f2py intent(in)  kpts, nb
    !f2py intent(out) dm_out
    !f2py integer, intent(hide),depend(kpts) :: nk=shape(kpts,0)

    if (nb /= dim_dm) then
       print *, 'ERROR get_dm num_band', dim_dm, " !=", nb
       STOP
    end if  
    call calc_dm(kpts, nk)
    dm_out(:,:,:) = dm(:,:,:)
    deallocate(dm)
  end subroutine get_dm


  !
  ! Common subroutine to calculate dynamical matrix
  !
  subroutine calc_dm(kpts, nk, skip_alloc)
    integer, intent(in)   :: nk
    real(8), intent(in):: kpts(nk, 3)
    logical, optional, intent(in)   :: skip_alloc
    !f2py intent(in) kpts
    !f2py integer, intent(hide),depend(kpts) :: nk=shape(kpts,0)  

    ! local variables
    integer ik, ipair, ix, jx, i, j, iA, jA
    real(8) tmp_sum, q0(3), kz1(3), kz2(3), q02, q0epsq0
    complex(kind=8) tmp, tmp1, mix_nac

    if (.not. present(skip_alloc)) then
       if (AllOCATED(dm)) then
          deallocate(dm)
       end if
       allocate(dm(dim_dm, dim_dm, nk))
    end if
    dm(:,:,:) = (0., 0.)
    do ik=1, nk
       q0(:) = kpts(ik, :)
       if ((nac_method >= 0).and.(nac_method <=2)) then
          q02 = dot_product(q0, q0)
          if ((q02 < 10E-30) .and.(nac_method /= 0)) then
             ! this is Gamma
             ! we need to know the direction of the symmetry line
             if ((ik>=2).and.(dot_product(kpts(ik-1,:),kpts(ik-1,:))>10E-12)) then
                q0(:) = kpts(ik-1,:)*1E-4
             else if ((ik<=nk-1).and.(dot_product(kpts(ik+1,:),kpts(ik+1,:))>10E-12)) then
                q0(:) =kpts(ik+1,:)*1E-4
             else
                ! Cannot determine the direction of the Gamma point
              !  q0(1)=0
              !  q0(2)=0
              !  q0(3)=1E-5
             end if
             q02 = dot_product(q0, q0)
          end if

          if (nac_method /= 0) then
             q0epsq0 = dot_product(matmul(q0, diel), q0)
             if (nac_method == 1) then
                mix_nac = exp(-q02/rho2)
             else
                mix_nac = cmplx(0,kind=8)
                do ipair = 1, Ncell_mixsp
                   tmp_sum = dot_product(q0, disp(idx_mixsp(ipair),:))
                   mix_nac = mix_nac + cmplx(cos(tmp_sum), sin(tmp_sum),kind=8)
                end do
                mix_nac = mix_nac/Ncell_mixsp
             end if
             !print *, 'debug q3=', q0, ' mix=', mix_nac

             do iA=1, natoms
                do i=1, 3
                   kz1(i) = sum(q0(:) * Zstar(:, i, iA))
                end do
                do jA=1, natoms
                   do i=1, 3
                      kz2(i) = sum(q0(:) * Zstar(:, i, jA))
                   end do
                   tmp = pref_nac*mix_nac/(sqrt(mass(iA)*mass(jA))* q0epsq0)
                   !                print *, iA, jA, tmp
                   do ix=1, 3
                      do jx=1, 3
                         dm(3*(iA-1)+ix, 3*(jA-1)+jx,ik) =dm(3*(iA-1)+ix, 3*(jA-1)+jx,ik) &
                              +tmp*kz1(ix)*kz2(jx)
                         dmtmp(3*(iA-1)+ix, 3*(jA-1)+jx) = tmp*kz1(ix)*kz2(jx)
                      end do
                   end do
                end do
             end do
          else 
             !   write(*,'(A,3F10.5)') 'DEBUG q=', q0
!             call calc_dm_dipole(q0, dmtmp, natoms)
!             dm(:,:,ik)= dm(:,:,ik)+ dmtmp
             call calc_dm_dipole(q0, dm(:,:,ik), natoms)
          end if
       end if
       call calc_dm_sr(q0, dm(:,:,ik), natoms, nfcm, disp, fcm, typ)
       ! write(*,*) nac_method
       ! write(*,'(12f7.3)') (dmtmp(:,iA), iA=1,3*natoms)
       ! if (dot_product(q0,q0)<1E-14) then
       !    write(*,'(12f7.3)') (dm(:,iA,ik), iA=1,3*natoms)
       ! end if
#ifdef DEBUG
       open (unit = 1, file = "dmat_all.dat")
       do iA=1,natoms*3
          do jA=1,natoms*3
             write(1, *) dm(iA, jA, ik)
          end do
       end do
       close(1)
#endif
    end do

  end subroutine calc_dm


  ! Short-range part using force constants
  subroutine calc_dm_sr(q, dmat, nat, nf, dis, fc, tp)
    implicit none
    real(8), intent(in) :: q(3)
    integer, intent(in) :: nat, nf
    complex(kind=8), intent(out) :: dmat(3*nat,3*nat)
    real(8), intent(in) :: dis(nf,3)
    real(8), intent(in) :: fc(nf,3,3)
    integer, intent(in) :: tp(nf,2)
    
    ! local variables
    integer ipair, ix, jx
    real(8) tmp_sum
    complex(kind=8) tmp, fac
    logical diff_pos

    do ipair=1, nf
       tmp_sum = dot_product(q, dis(ipair,:))
       diff_pos= any(abs(dis(ipair,:))>1E-12)
       fac = cmplx(cos(tmp_sum),sin(tmp_sum),kind=8)/sqrt(mass(tp(ipair,1))* mass(tp(ipair,2)))
       do ix=1, 3
          do jx=1, 3
             tmp=fc(ipair, ix, jx)* fac
             dmat(3*(tp(ipair,1)-1)+ix, 3*(tp(ipair,2)-1)+jx) = &
                  dmat(3*(tp(ipair,1)-1)+ix, 3*(tp(ipair,2)-1)+jx)+tmp
             !print *, ix, jx, tmp, dmat(3*(tp(ipair,1)-1)+ix,3*(tp(ipair,2)-1)+jx) 
             if (diff_pos .OR.(tp(ipair,1) /= tp(ipair,2))) then
                dmat(3*(tp(ipair,2)-1)+jx, 3*(tp(ipair,1)-1)+ix) = &
                     dmat(3*(tp(ipair,2)-1)+jx, 3*(tp(ipair,1)-1)+ix) + conjg(tmp)
             end if
             !print *,'a', ix, jx, tmp,dmat(3*(tp(ipair,1)-1)+ix,3*(tp(ipair,2)-1)+jx)
          end do
       end do
    end do
  end subroutine calc_dm_sr



  !
  ! calculate eigen energies
  !
  subroutine get_dispersion(kpts, nk, connect_pol, eigE, nb)
    integer, intent(in)   :: nk, nb
    real(8), intent(in)   :: kpts(nk, 3)
    integer, intent(in)   :: connect_pol
    real(8),intent(out)   :: eigE(nb, nk)
    !f2py intent(in)  kpts, connect_pol, nb
    !f2py intent(out) eigE
    !f2py integer, intent(hide),depend(kpts) :: nk=shape(kpts,0)  

    ! local variables
    integer ik, info, i

    !   print *, 'DEBUG get_dispersion nband=', nband, ' nk=', nk

    call calc_dm(kpts, nk)

    do ik=1, nk
       !       print *, 'DM at ', ik, dm(ik,:,:)
       call ZHEEV('N', 'L', dim_dm, dm(:,:,ik), dim_dm, W, WORK, lw, RWORK, info)
       if (info/=0) then
          print *, 'ERROR ZHEEV failed: ', info
          print *, 'debug dm ', dm(:, :, ik)
          print *, 'debug eig ', W(:)
       end if
       !       print *, 'W ', W(:)
       do i=1, dim_dm
          eigE(i, ik) = sign(sqrt(abs(W(i))), W(i)) * UNITTHZ
       end do
    end do
    !    print *, 'f90 eigE', eigE
    deallocate(dm)

  end subroutine get_dispersion





  subroutine get_eig_e_vec(kgrid, nk, nb, eigE, eigV)
    integer, intent(in)   :: nk, nb
    real(8), intent(in):: kgrid(nk, 3)
    real(8), intent(out):: eigE(nb, nk)
    complex(kind=8), intent(out):: eigV(nb, nb, nk)
    !f2py intent(in)  kgrid, nb
    !f2py intent(out) eigE, eigV
    !f2py integer, intent(hide),depend(kgrid) :: nk=shape(kgrid,0)

    ! local variables
    integer ik, i, info

    if (nb /= dim_dm) then
       print *, 'input nb of get_eig_e_vec should be', dim_dm
       stop
    end if
    if (AllOCATED(dm)) then
       deallocate(dm)
    end if
    allocate(dm(dim_dm, dim_dm, 1))
    do ik=1, nk
       call calc_dm(kgrid(ik,:), 1, .true.)
       !       print *, 'DM at ', ik, dm(ik,:,:)
       call ZHEEV('V', 'L', dim_dm, dm(:,:,1), dim_dm, W, WORK, lw, RWORK, info)
       if (info/=0) then
          print *, 'ERROR ZHEEV failed: ', info
       end if
       !       print *, 'W ', W(:)
       do i=1, dim_dm
          eigE(i, ik) = sign(sqrt(abs(W(i))), W(i)) * UNITTHZ
       end do
       eigV(:, :, ik) = dm(:,:,1)
    end do

    deallocate(dm)
  end subroutine get_eig_e_vec

  subroutine set_dos_en(eigE, density, width, nk, nen, dim)
    implicit none
    integer, intent(in) :: nk, nen, dim
    real(8), intent(out):: density(nen, 2), width
    real(8), intent(in) :: eigE(dim, nk)
    !local var
    real(8) emin, emax, padding
    integer i

    emin = minval(eigE)
    emax = maxval(eigE)
    padding = (emax-emin)/nen
    emin = emin-padding
    emax = emax+padding
    !    print *, 'emin, emax, width', emin, emax, width
    width = (emax - emin)/(nen-1)
    do i=1, nen
       density(i, 1) = emin + (i-1)*width
       density(i, 2) = 0
    end do
  end subroutine set_dos_en



  ! thermal properties. 
  ! Assuming energy in eV
  subroutine calc_thermal(den, nen, Tlist, nT, dat)
    integer, intent(in)   :: nT, nen
    real(8), intent(in):: den(nen, 2)
    real(8), intent(in):: Tlist(nT)
    real(8), intent(out):: dat(nT, 5)
    !f2py intent(in)  den, range
    !f2py intent(out) dat
    !f2py integer, intent(hide),depend(den) :: nen=shape(den,0)
    !f2py integer, intent(hide),depend(Tlist) :: nT=shape(Tlist,0)

    ! local variables
    integer i, idos
    real(8) t, en, efree, entropy, cv, dos, nph, expVal, hf, kBT, ewidth

    if (nen<=1) then
       print *, 'WARNING: need at least two DOS data points'
    end if
    ewidth= den(2,1)- den(1,1)

    do i=1, nT
       t = Tlist(i)
       en=0.0
       efree=0.0
       cv=0.0
       do idos= 1, nen
          hf = den(idos, 1)
          dos = den(idos, 2)
          if (dos < 1E-16) cycle
          if (hf <0) then
             print *, 'WARNING: imaginary phonon encountered and ignored'
             cycle
          end if
          kBT= kBbyeV*t
          expVal= exp(hf/kBT)
          nph=1.0/(expVal-1.0)

          en=en+ hf * (0.5+ nph) *dos
          ! heat capacity in kB
          cv=cv+ (hf*nph/kBT)**2 *expVal *dos;
          efree=efree+ (0.5*hf+ kBT*log(1.0-1.0/expVal))*dos
          ! print *, 'debug en,efree=', en, efree
       end do
       ! entropy in kB
       entropy=(en-efree)/kBT
       dat(i, 1) = t
       dat(i, 2) = en*ewidth
       dat(i, 3) = efree*ewidth
       dat(i, 4) = entropy*ewidth
       dat(i, 5) = cv*ewidth
    end do


  end subroutine calc_thermal




  ! density of states
  ! energy in THz, so epsilon should also be in THz
  ! ismear= -1 using the tetrahedron method
  !       = 0 (Gauss), 1 (Lorentzian) smearing
  subroutine get_dos_new(do_pdos,ndiv, kgrid, nkin, nen, ismear, epsilon, natom, density, pdos)
    integer, intent(in)   :: nkin, nen, ismear, natom, ndiv(3)
    logical, intent(in) :: do_pdos
    real(8), intent(in) :: kgrid(nkin, 3),epsilon
    real(8), TARGET, intent(out):: density(nen, 2)
    real(8), intent(out):: pdos(natom, nen)
    !f2py intent(in)  do_pdos, ndiv, kgrid, nen, ismear, natom,epsilon
    !f2py intent(out) density, pdos
    !f2py integer, intent(hide),depend(kgrid) :: nkin=shape(kgrid,0)

    ! local variables
    real(8), allocatable :: eigE(:,:), proj(:,:,:)  ! eigE(iband, ik), proj(iatom, iband, ik)
    real(8) width, fac, smear_wt
    integer nk, ik, i, j, k, info
    complex(kind=8),allocatable :: eigv(:,:,:)
    real(8), pointer     :: XE(:), Y(:) ! energies for DOS plot, and DOS values
    ! tetrahedron var
    integer Ntet, Nptk, Nkmx, Ntmx
    integer, allocatable :: Idef(:,:)
    real(8), allocatable  :: Ptk(:,:), intdos(:), intdosp(:,:)




    if (ismear==-1) then
       ! tetrahedron method
       Nkmx = ndiv(1)*ndiv(2)*ndiv(3)*8
       Ntmx = Nkmx*6
       allocate(Idef(5, Ntmx), Ptk(4, Nkmx), intdos(nen),intdosp(natom, nen))
       call SETK06(ndiv(1),ndiv(2),ndiv(3),avec(:,1),avec(:,2),avec(:,3),Ptk,Nptk,Idef,Ntet,Nkmx,Ntmx)
       nk = Nptk
       allocate(eigE(dim_dm, nk), proj(natom,dim_dm,nk), eigv(dim_dm,dim_dm,nk))
       call get_eig_e_vec(transpose(Ptk(1:3,1:nk)), nk, dim_dm, eigE, eigV)
    else
       nk = nkin
       allocate(eigE(dim_dm, nk), proj(natom,dim_dm,nk), eigv(dim_dm,dim_dm,nk))
       call get_eig_e_vec(kgrid, nk, dim_dm, eigE, eigV)
    end if

    do ik=1,nk
       do j=1, dim_dm
          do k=1, natom
             proj(k, j, ik)= sum(abs(eigV(3*(k-1)+1:3*k, j, ik))**2)
          end do
       end do
    end do
    deallocate(eigv)
    call set_dos_en(eigE, density, width, nk, nen, dim_dm)
    XE => density(:,1)
    Y => density(:,2)
    if (do_pdos) pdos=0

    if (ismear==-1) then
       ! tetrahedron method
       if (do_pdos) then
          call PDSTET(proj,natom,dim_dm,natom,eigE,dim_dm,dim_dm,Idef,    &
               & Ntet,XE,nen,Y,intdos,pdos,intdosp,natom)
       else
          call DOSTET(eigE,dim_dm,dim_dm,Idef,Ntet,XE,nen,Y,intdos)
       end if
    else
       ! smearing
       do ik=1, nk
          do i=1, nen
             do j=1, dim_dm
                if (ismear==1) then
                   ! Lorentzian
                   smear_wt = 1/(((eigE(j,ik)-XE(i))/epsilon)**2 + 1)
                else
                   ! Gauss, ismear=0
                   smear_wt = exp(-0.5* ((eigE(j,ik)-XE(i))/epsilon)**2)
                end if
                do k=1, natom
                   pdos(k, i)=pdos(k, i)+ smear_wt* proj(k, j, ik)
                end do
                Y(i)=Y(i)+ smear_wt
             end do
          end do
       end do
    end if

    do i=1, natom
       fac = 3/(width* sum(pdos(i,:)))
       pdos(i,:)= pdos(i,:)*fac
    end do
    deallocate(eigE, proj)
    if (allocated(Idef)) deallocate(Idef)
    if (allocated(Ptk)) deallocate(Ptk)
    if (allocated(intdos)) deallocate(intdos)
    if (allocated(intdosp)) deallocate(intdosp)
  end subroutine get_dos_new



  ! wrapper for DM_dipole_dipole
  ! using stored cell info
  subroutine calc_dm_dipole(q, dmat, nat)
    implicit none
    real(8), intent(in) :: q(3)
    integer, intent(in):: nat
    complex(kind=8), intent(out) :: dmat(3*nat,3*nat)

    !  complex(kind=8), pointer :: Amat(3,natoms,3,natoms)
    real(8)  errlim
    integer iA, jA, ix, jx

    if (nat /= natoms) STOP 'ERROR nat != natoms'
    errlim = 1d-18
    call get_fcm_dipole(avec, qvec, errlim, apos, nat, Zstar, diel, q, dmat)
    do iA=1, nat
       do ix=1,3
          do jA=1, nat 
             do jx=1,3
                dmat(3*(iA-1)+ix, 3*(jA-1)+jx) =dmat(3*(iA-1)+ix, 3*(jA-1)+jx)/sqrt(mass(iA)*mass(jA))
             end do
          end do
       end do
    end do

    ! Long-range ASR correction using force constants fcm_c
    if (allocated(fcm_c)) then
#ifdef DEBUG
       dmtmp(:,:)=dmat(:,:)
       open (unit = 1, file = "dmat_lr.dat")
       do ix=1,nat*3
          do jx=1,nat*3
             write(1, *) dmtmp(ix, jx)
          end do
       end do
       close(1)
#endif
       call calc_dm_sr(q, dmat, natoms, nfcm_c, disp_c, fcm_c, typ_c)
#ifdef DEBUG
       dmtmp= dmat-dmtmp
       open (unit = 1, file = "dmat_dpcor.dat")
       do ix=1,nat*3
          do jx=1,nat*3
             write(1, *) dmtmp(ix, jx)
          end do
       end do
       close(1)
#endif
    end if
  end subroutine calc_dm_dipole


  ! stand-alone wrapper for DM_dipole_dipole
  ! for calculation of dipole DM of supercell
  ! all input passedly directly, NOT from the module
  subroutine get_fcm_dipole(rb, qb, errlim, apos, nat, Zeff, eps, qph, dmat)
    implicit none
    integer, intent(IN)                    :: nat
    real(8), dimension(3,3), intent(IN)     :: rb, qb, eps
    real(8), intent(IN)                     :: errlim
    real(8), dimension(3,nat), intent(IN)   :: apos
    real(8), dimension(3,3,nat), intent(IN) :: Zeff
    real(8), dimension(3), intent(IN)       :: qph
    complex(kind=8), TARGET, dimension(3*nat,3*nat), intent(OUT) :: dmat
    !f2py intent(in) rb, qb, errlim, apos, Zeff, eps, qph
    !f2py intent(out) dmat
    !f2py integer, intent(hide),depend(apos) :: nat=shape(apos,1)  

    real(8)  eta
!    integer trp
    !  complex(kind=8), pointer :: Amat(3,natoms,3,natoms)

    !  Amat=>dmat
    dmat = (0d0, 0d0)
    eta= 4d0
    call find_Ewald_eta_screened(rb, errlim, nat, apos, eps, eta)
!    open (unit=99, file='trp-fac.txt', status='old', action='read')
!    read(99, *), trp, etafac
!    close(99)
!    print *, 'calling DM_dipole_dipole eta=', eta
    call DM_dipole_dipole(rb, qb, eta, errlim, apos, nat, Zeff, eps, qph/(2*PI),dmat)
!    print *, 'calling DMstandard'
!    print *, 'debug trp=', trp
!    print *, 'qph=',qph, ' max(abs(D-D*))', maxval(abs(dmat-conjg(dmat)))
!    if (trp .EQ. 1) then
!       print *, 'DM_dp conjugated'
!       dmat = conjg(dmat)
!    end if
    call DMstandard(-qph/(2*PI), dmat, nat, apos)
    dmat = dmat * atomic2eVAng
  end subroutine get_fcm_dipole



  ! wrapper to get ionic forces Fewald
  subroutine calc_force_ion(rb, errlim, atp, nat, Zh, frc)
    implicit none
    integer, intent(in) :: nat
    real(8), intent(in)  :: rb(3,3), atp(3,nat), Zh(nat), errlim
    real(8), intent(out) :: frc(3,nat)

    real(8) eta

    call find_Ewald_eta_screened(avec, errlim, nat, apos, diel, eta)
    call Fewald(rb, eta, errlim, atp, nat, Zh, frc)
    frc = frc * atomic2eVAng
    !         call Eion(rb, eta, errlim, apos(:,:,i), Natoms, Zbare, ecc)
    !         ecc=ecc*Har2eV/ang2bohr
  end subroutine calc_force_ion

end module f_phonon
