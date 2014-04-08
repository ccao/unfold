MODULE input
  !
  use constants
  !
  implicit none
  !
  character(len=80) seed
  real(dp), dimension(1:3,1:3) :: scell      ! Super Cell transformation matrix (from UC to SC)
  real(dp) elow, ehigh, de
  real(dp) beta
  character(len=3) fmode                    ! Hamiltonian file mode: n: original wannier90 output r: reduced format
  integer nkpt                               ! Number of output kpoints (in UC)
  real(dp), allocatable :: kvec(:, :)        ! Output K-points
  !
CONTAINS
  !
  ! #1: seed
  ! #2: elow ehigh de
  ! #3: mode (0: single kpt; 1: band structure; 2: plane-cut; 3: whole space)
  ! #4:  if (mode==0) kvec; if (mode==1) nksec nk_per_sec; otherwise nkx nky nkz
  !
  SUBROUTINE read_input()
    !
    use constants
    use para
    use specdata,   only : nen, ene, spec
    !
    implicit none
    !
    integer ii, iksec, ik
    integer nksec, nk_per_sec
    integer nkx, nky, nkz
    integer ikx, iky, ikz
    integer mode
    !
    real(dp), allocatable :: kbnd_vec(:, :)
    !
    if (inode.eq.0) then
      !
      open(unit=fin, file="unfold.inp")
      !
      read(fin, *) seed, fmode
      !
      read(fin, *) elow, ehigh, de
      read(fin, *) beta
      !
      nen=(ehigh-elow)/de+1
      de=(ehigh-elow)/(nen-1)
      !
      read(fin, *) mode
      !
      write(stdout, '(1A,1A20)') " #: Seed name:", seed
      write(stdout, '(1A,F12.8,1A,F12.8,1A)'), " #:    energy range [", elow, ",", ehigh, "]"
      write(stdout, '(1A,F9.6)') " #:   with beta:", beta
      !
    endif
    !
    CALL para_sync(scell, 3, 3)
    CALL para_sync(elow)
    CALL para_sync(ehigh)
    CALL para_sync(de)
    CALL para_sync(beta)
    CALL para_sync(nen)
    CALL para_sync(mode)
    !
    allocate(ene(1:nen))
    allocate(spec(1:nen))
    !
    do ii=1, nen
      ene(ii)=elow+(ii-1.d0)*de
    enddo
    !
    select case (mode)
      case (0)
        nkpt=1
        allocate(kvec(1:3, 1:1))
        if (inode.eq.0) then
          read(fin, *) kvec(:, 1)
          write(stdout, *) "   #   : single point mode"
          write(stdout, '(1A,F12.8,1A,F12.8,1A,F12.8,1A)') "   #   @ (", kvec(1,1), ",", kvec(2,1), ",", kvec(3,1), ")"
        endif
        CALL para_sync(kvec, 3, nkpt)
      case (1)
        if (inode.eq.0) read(fin, *) nksec, nk_per_sec
        CALL para_sync(nksec)
        CALL para_sync(nk_per_sec)
        nkpt=(nksec-1)*nk_per_sec+1
        allocate(kvec(1:3, nkpt))
        if (inode.eq.0) then
          write(stdout, *) "   #   : Band structure mode"
          write(stdout, '(1A,I5,1A)') "   #   :  with total of ", nkpt, " K-points"
          allocate(kbnd_vec(1:3, 1:nksec))
          do ii=1, nksec
            read(fin, *) kbnd_vec(:, ii)
          enddo
          do ii=1, nksec-1
            do ik=1, nk_per_sec
              kvec(:, (ii-1)*nk_per_sec+ik)=kbnd_vec(:, ii)+(kbnd_vec(:, ii+1)-kbnd_vec(:, ii))*(ik-1)/nk_per_sec
            enddo
          enddo
          kvec(:, nkpt)=kbnd_vec(:, nksec)
          deallocate(kbnd_vec)
        endif
        CALL para_sync(kvec, 3, nkpt)
      case (2)
        if (inode.eq.0) read(fin, *) nkx, nky
        CALL para_sync(nkx)
        CALL para_sync(nky)
        nkpt=nkx*nky
        allocate(kvec(1:3, 1:nkpt))
        do ikx=1, nkx
          do iky=1, nky
            ik=(iky-1)*nkx+ikx
            kvec(1, ik)=(ikx-1.d0)/nkx
            kvec(2, ik)=(iky-1.d0)/nky
            kvec(3, ik)=0.d0
          enddo
        enddo
      case (3)
        if (inode.eq.0) read(fin, *) nkx, nky, nkz
        CALL para_sync(nkx)
        CALL para_sync(nky)
        CALL para_sync(nkz)
        nkpt=nkx*nky*nkz
        allocate(kvec(1:3, 1:nkpt))
        do ikx=1, nkx
          do iky=1, nky
            do ikz=1, nkz
              ik=(ikz-1)*nkx*nky+(iky-1)*nkx+ikx
              kvec(1, ik)=(ikx-1.d0)/nkx
              kvec(2, ik)=(iky-1.d0)/nky
              kvec(3, ik)=(ikz-1.d0)/nkz
            enddo
          enddo
        enddo
    end select
    !
    close(unit=fin)
    !
  END SUBROUTINE
  !
  SUBROUTINE finalize_input
    !
    implicit none
    !
    if (allocated(kvec)) deallocate(kvec)
    !
  END SUBROUTINE
END MODULE
