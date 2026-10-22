      module A1_gef_mod
      implicit none

      real, allocatable :: Etab(:,:)
      real, allocatable :: Jtab(:,:)
      real, allocatable :: Ytab(:)
      real, allocatable :: Emean(:)
      real, allocatable :: dEmean(:)

      contains

      subroutine gef_allocate(n)
      integer, intent(in) :: n

      call gef_deallocate()

      if (n <= 0) return

      allocate(Etab(n,0:1000))
      allocate(Jtab(n,1:100))
      allocate(Ytab(n))
      allocate(Emean(n))
      allocate(dEmean(n))

      Etab = 0.
      Jtab = 0.
      Ytab = 0.
      Emean = 0.
      dEmean = 0.
      end subroutine gef_allocate

      subroutine gef_deallocate()

      if (allocated(Etab)) deallocate(Etab)
      if (allocated(Jtab)) deallocate(Jtab)
      if (allocated(Ytab)) deallocate(Ytab)
      if (allocated(Emean)) deallocate(Emean)
      if (allocated(dEmean)) deallocate(dEmean)

      end subroutine gef_deallocate

      end module A1_gef_mod
