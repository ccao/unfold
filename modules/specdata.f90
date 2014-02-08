MODULE specdata
  !
  use constants
  !
  implicit none
  !
  integer nen
  real(dp), allocatable :: ene(:)
  real(dp), allocatable :: spec(:)
  real(dp), allocatable :: eig(:)
  complex(dp), allocatable :: work(:, :)
  !
CONTAINS
  !
  SUBROUTINE init_spec
    !
    use wanndata,  only : norb
    !
    implicit none
    !
    allocate(eig(1:norb))
    allocate(work(1:norb, 1:norb))
    !
  END SUBROUTINE

  SUBROUTINE finalize_spec
    !
    implicit none
    !
    if (allocated(ene)) deallocate(ene)
    if (allocated(spec)) deallocate(spec)
    if (allocated(eig)) deallocate(eig)
    if (allocated(work)) deallocate(work)
    !
  END SUBROUTINE
  !
END MODULE
