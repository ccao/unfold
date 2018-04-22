!
program unfold
  !
  use para,        only : init_para, finalize_para, inode
  use constants,   only : dp, fout
  use wanndata,    only : read_ham, finalize_wann, read_reduced_ham
  use mapdata,     only : init_mapping, finalize_mapping
  use specdata,    only : nen, init_spec, finalize_spec
  use input,       only : kvec, nkpt, elow, ehigh, seed, finalize_input, read_input, fmode
  !
  implicit none
  !
  integer ik
  !
  call init_para
  call read_input
  if(fmode(1:1).eq.'r') then
    call read_reduced_ham(seed)
  else
    call read_ham(seed)
  endif
  call init_spec
  call init_mapping
  !
  if (inode.eq.0) then
    open(unit=fout, file="unfold.dat")
    write(fout, '(2I10,2F24.16)') nkpt, nen, elow, ehigh
  endif
  !
  do ik=1, nkpt
    !
    if (inode.eq.0) write(fout, '(3F16.8)') kvec(:, ik)
    call calc_spectrum_k(kvec(:,ik))
    call output_spectrum
    !
  enddo ! ik
  !
  close(unit=fout)
  !
  call finalize_input
  call finalize_spec
  call finalize_wann
  call finalize_mapping
  call finalize_para
  !
end program
