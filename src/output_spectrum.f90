SUBROUTINE output_spectrum
  !
  use para,      only : inode
  use constants, only : fout
  use specdata,  only : spec, nen
  use input,     only : nkpt
  !
  implicit none
  !
  integer ii
  !
  if (inode.eq.0) then
    if (nkpt.eq.1) then
      write(fout, '(1F12.8)') (spec(ii), ii=1, nen)
    else
      write(fout, '(10F16.9)') (spec(ii), ii=1, nen)
    endif
  endif
  !
END SUBROUTINE
