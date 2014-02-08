include 'lapack.f90'

SUBROUTINE calc_spectrum_k(kv_uc)
  !
  use lapack95,    only : heev
  use constants,   only : dp, cmplx_0, cmplx_i, twopi
  use mapdata,     only : ksi, rv, norb_uc, nn, iiorb 
  use specdata,    only : spec, nen, ene, eig, work
  use wanndata,    only : norb, nrpt, rvec, weight, ham
  !
  implicit none
  !
  real(dp), dimension(1:3) :: kv_uc  ! (unfolded) K-point in UC's IBZ
  !
  real(dp) kv(1:3)    ! Folded K-point in SC's IBZ
  !
  real(dp) ratio, delta, ibz
  complex(dp) fact
  !
  integer ii, jorb, iorb, iorb_uc, info
  !
  spec(:)=0.d0
  !
  work(:,:)=cmplx_0
  !
  do ii=1, 3
    kv(ii)=ibz(SUM(ksi(:, ii)*kv_uc(:)))
  enddo
  !
  do ii=1, nrpt
    ratio=SUM(kv(:)*rvec(:, ii))
    fact=exp(-cmplx_i*twopi*ratio)/weight(ii)
    work(:,:)=work(:,:)+fact*ham(:,:,ii)
  enddo
  !
  call heev(work, eig, 'V', 'U', info)
  !
  ! Unfolding according to:
  !
  ! A_k(\omega)=\Sum_n A_kn
  ! A_kn(\omega)=\Sum_{KJ} |<kn|KJ>|^2 A_{KJ}(\omega)
  !   and 
  !     <kn|KJ>=\sqrt(L/l)\Sum_N exp(-i k r(N)) \delta_{n,[N]}\delta_{[k],K} <KN|KJ>
  !    where k & K denotes kpoints in UC BZ & SC BZ respectively
  !          n & N denotes band index in UC & SC respectively
  !         [k] means folded kpoints in SC BZ corresponding to k 
  !         [N] means unfolded band index in UC corresponding to N
  !
  do jorb=1, norb           ! loop over J  : jorb ==> J
    !
    ratio=0.d0              ! ratio : ==> \Sum_n |<kn|KJ>|^2
    do iorb_uc=1, norb_uc   ! loop over n  : iorb_uc ==> n
      !
      fact=cmplx_0          ! fact: <kn|KJ>
      do ii=1, nn(iorb_uc)  ! loop over all N : trick: only N's that [N]=n
        iorb=iiorb(ii, iorb_uc) ! reverse finding N
        fact=fact+exp(-cmplx_i*twopi*SUM(kv_uc(:)*rv(:, iorb)))*work(iorb, jorb)
      enddo ! ii
      ratio=ratio+conjg(fact)*fact
      !
    enddo ! iorb_uc
    !
    do ii=1, nen
      spec(ii)=spec(ii)+ratio*delta(ene(ii)-eig(jorb))
    enddo ! ii
    !
  enddo ! jorb
  !
END SUBROUTINE

real(dp) FUNCTION delta(en)
  !
  use constants, only : dp
  use input,     only : beta
  !
  implicit none
  !
  real(dp) en
  !
  delta=beta**2/(en**2+beta**2)
  !
  return
  !
END FUNCTION

real(dp) FUNCTION ibz(x)
  !
  use constants, only : dp
  !
  implicit none
  !
  real(dp) x
  ibz=x-floor(x)
  return
  !
END FUNCTION
