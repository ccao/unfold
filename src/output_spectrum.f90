SUBROUTINE output_spectrum
  !
  use para,      only : inode
  use constants, only : stdout
  use specdata,  only : spec, nen
  use input,     only : nkpt
  !
  implicit none
  !
  integer ii
  !
  if (inode.eq.0) then
    if (nkpt.eq.1) then
      write(stdout, '(1F12.8)') (spec(ii), ii=1, nen)
    else
      write(stdout, '(10F16.9)') (spec(ii), ii=1, nen)
    endif
  endif
  !
END SUBROUTINE
