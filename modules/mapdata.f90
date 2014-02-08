MODULE mapdata
  !
  use constants,  only : dp
  !
  implicit none
  !
  real(dp), dimension(1:3, 1:3) :: scell
  real(dp), dimension(1:3, 1:3) :: ksi
  real(dp), allocatable :: rv(:, :)
  integer, allocatable :: nn(:), iiorb(:, :)
  integer norb_uc, norb_sc
  !
 CONTAINS
  !
  SUBROUTINE init_mapping
    !
    use para,       only : inode, para_sync
    use constants,  only : dp, fin, stdout
    use wanndata,   only : norb
    !
    implicit none
    !
    integer ii, t1
    !
    if (inode.eq.0) then
      !
      open(unit=fin, file="unfold.map")
      do ii=1, 3
        read(fin, *) scell(:, ii)
      enddo
      !  Transformation matrix from UC lattice to SC lattice : A'=scell*A 
      read(fin, *) norb_uc, norb_sc
      !  Number of orbitals in single Unit cell; number of orbitals in supercell
      if (norb.ne.norb_sc) then
        write(stdout, *) "  !!! ERROR: The hamiltonian dimension is inconsistent with mapping file"
        stop
      endif
      !
      write(stdout, *) "   #   : Unit cell to supercell transformation"
      do ii=1, 3
        write(stdout, '(A,I1,A,I2,A,I2,A,I2,A)') "    #   : A'", ii, "=", int(scell(1, ii)), "A1 +", int(scell(2, ii)), "A2 +", int(scell(3, ii)), "A3"
      enddo
      !
      write(stdout, *) "   #   : === Norb_uc ===  :  === Norb_sc ==="
      write(stdout, '(A3,10X,I5,15X,I5)') "#", norb_uc, norb_sc
      !
    endif
    !
    call para_sync(scell, 3, 3)
    call para_sync(norb_uc)
    call para_sync(norb_sc)
    !
    allocate(rv(1:3, 1:norb_sc))
    allocate(nn(1:norb_uc))
    allocate(iiorb(1:norb_sc, 1:norb_uc))
    !
    !  then BZ from UC to SC : B'=scell^-1*A, so sigma=scell^-1
    !     But, then kpt coordinates translate according to ksi=scell
    !
    ksi(:,:)=scell(:,:)
    nn(:)=0
    iiorb(:,:)=-1
    if (inode.eq.0) then
      !
      do ii=1, norb_sc
        read(fin, *) t1, rv(:,ii)   ! NC orb index; r(N)
        nn(t1)=nn(t1)+1
        iiorb(nn(t1), t1)=ii
      enddo
      close(unit=fin)
      !
      write(stdout, *) "  # :==== MAPPING OF PRIMITIVE CELL ORBITALS TO SUPERCELL ===="
      do ii=1, norb_uc
        write(stdout, *) "  # ", ii, ":", nn(ii), ":",(iiorb(t1, ii), t1=1, nn(ii))
      enddo
      write(stdout, *) "  # :========================================================="
      !
    endif
    !
    call para_sync(rv, 3, norb_sc)
    call para_sync(nn, norb_uc)
    call para_sync(iiorb, norb_sc, norb_uc)
    !
  END SUBROUTINE
  !
  SUBROUTINE finalize_mapping
    !
    implicit none
    !
    if ( allocated(rv) ) deallocate(rv)
    if ( allocated(nn) ) deallocate(nn)
    if ( allocated(iiorb) ) deallocate(iiorb)
    !
  END SUBROUTINE
  !
END MODULE
