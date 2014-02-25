!
program unfold
  !
  use para,        only : init_para, finalize_para
  use constants,   only : dp, stdout
  use wanndata,    only : read_ham, finalize_wann
  use mapdata,     only : init_mapping, finalize_mapping
  use specdata,    only : nen, init_spec, finalize_spec
  use input,       only : kvec, nkpt, seed, finalize_input, read_input
  !
  implicit none
  !
  integer ik
  !
  call init_para
  call read_input
  call read_ham(seed)
  call init_spec
  call init_mapping
  !
  write(stdout, '(2I10)') nkpt, nen
  !
  do ik=1, nkpt
    !
    write(stdout, '(3F16.8)') kvec(:, ik)
    call calc_spectrum_k(kvec(:,ik))
    call output_spectrum
    !
  enddo ! ik
  !
  call finalize_input
  call finalize_spec
  call finalize_wann
  call finalize_mapping
  call finalize_para
  !
end program
