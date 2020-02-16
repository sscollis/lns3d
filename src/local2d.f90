!=============================================================================!
        module local2d
!
!  Local variables for two-dimensional analysis
!
!=============================================================================!
        use global

        real, allocatable :: g1v(:,:,:), g2v(:,:,:)
        real, allocatable :: g11v(:,:,:), g12v(:,:,:), g22v(:,:,:)
        !$sgi distribute g1v(*,*,block), g2v(*,*,block)
        !$sgi distribute g11v(*,*,block), g12v(*,*,block), g22v(*,*,block)

        real :: gu(3,3), grho(3), gt(3), gp(3)
                             
        real :: divu, g1divu, g2divu, g3divu
        
        real :: rho, u1, u2, u3, t, p, rhoinv

        real :: mu,     dmu,    d2mu
        real :: lm,     dlm,    d2lm
        real :: con,    dcon,   d2con
        real :: g1mu,   g2mu,   g3mu
        real :: g1lm,   g2lm,   g3lm
        real :: g1con,  g2con,  g3con
        real :: g1dmu,  g2dmu,  g3dmu
        real :: g1dlm,  g2dlm,  g3dlm
        real :: g1dcon, g2dcon, g3dcon

        real :: S1jj, S2jj, S3jj, Lapt

        real :: S(3,3)

        real, allocatable :: Ah(:,:,:,:), Bh(:,:,:,:), Dh(:,:,:,:)
        
        real, allocatable :: Vh11(:,:,:), Vh12(:,:,:), Vh22(:,:,:)

        contains

!=============================================================================!
        subroutine mlocal2d
!
!  Allocate variables for 3-d analysis
!
!=============================================================================!
        implicit none
        integer :: i, j, idof, ier
        character*80 :: code='mLocal2$'
!=============================================================================!
        allocate( g1v(ndof,nx,ny), g2v(ndof,nx,ny), g11v(ndof,nx,ny), &
                  g12v(ndof,nx,ny), g22v(ndof,nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory$')

        !$doacross local(i,idof)
        !$omp parallel do private(i,idof)
        do j = 1, ny
          do i = 1, nx
            do idof = 1, ndof
              g1v(idof,i,j) = zero
              g2v(idof,i,j) = zero
              g11v(idof,i,j) = zero
              g12v(idof,i,j) = zero
              g22v(idof,i,j) = zero
            end do
          end do
        end do

!.... the following variables are only needed for linear or implicit problems

        if (linear.eq.1 .or. impl.ne.0) then
          allocate( Ah(ndof,ndof,nx,ny), Bh(ndof,ndof,nx,ny),   &
                    Dh(ndof,ndof,nx,ny), STAT=ier )
          if (ier .ne. 0) call error(code,'Insufficient Memory$')
          allocate( Vh11(6,nx,ny), Vh12(6,nx,ny), Vh22(6,nx,ny), STAT=ier )
          if (ier .ne. 0) call error(code,'Insufficient Memory$')
        end if

        end subroutine mlocal2d

        end module local2d
!=============================================================================!
