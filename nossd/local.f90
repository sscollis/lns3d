!=============================================================================!
        module local
!
!  Local variables
!
!=============================================================================!

        real, allocatable :: g1v(:,:), g2v(:,:)
        real, allocatable :: g11v(:,:), g12v(:,:), g22v(:,:)
        
        real, allocatable :: g1vl(:), g2vl(:)
        real, allocatable :: g11vl(:), g12vl(:), g22vl(:)
        
        real, allocatable :: Ah(:,:,:), Bh(:,:,:), Dh(:,:,:)
        
        real, allocatable :: Vh11(:,:), Vh12(:,:), Vh22(:,:)

        real, allocatable :: ABhi(:,:), Dhi(:,:,:)

!=============================================================================!

        real, allocatable :: gu(:,:,:), grho(:,:),              &
                             gt(:,:),   gp(:,:)
                             
        real, allocatable :: divu(:),   g1divu(:),              &
                             g2divu(:), g3divu(:)
        
        real, allocatable :: rho(:), u1(:), u2(:),              &
                             u3(:),  t(:),  p(:),               &
                             rhoinv(:)

        real, allocatable :: mu(:),     dmu(:),    d2mu(:)
        real, allocatable :: lm(:),     dlm(:),    d2lm(:)
        real, allocatable :: con(:),    dcon(:),   d2con(:)
        real, allocatable :: g1mu(:),   g2mu(:),   g3mu(:)
        real, allocatable :: g1lm(:),   g2lm(:),   g3lm(:)
        real, allocatable :: g1con(:),  g2con(:),  g3con(:)
        real, allocatable :: g1dmu(:),  g2dmu(:),  g3dmu(:)
        real, allocatable :: g1dlm(:),  g2dlm(:),  g3dlm(:)
        real, allocatable :: g1dcon(:), g2dcon(:), g3dcon(:)

        real, allocatable :: S1jj(:),  S2jj(:),  S3jj(:),       &
                             S(:,:,:), Lapt(:)

        end module local
!=============================================================================!
