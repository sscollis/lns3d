!============================================================================!
        program mkpot 
!  
!  THIS CODE IS FOR STREET'S INCOMPRESSIBLE SOLVER!!!!!
!
!  Extract boundary data based on pot field.  Note that the pressure is 
!  assumed to be ndof=5, and the pressure gradient in xi is output.
!
!  Also note that the velocity components are in cartisian frame
!  
!  Author:  Scott Collis
!
!  Revised: 7-23-97
!============================================================================!
        use const
        implicit none

!.... flow data

        real, allocatable :: v(:,:,:)
        
        real, allocatable :: g1v(:,:,:), g2v(:,:,:)

!.... mesh

        real, allocatable :: x(:,:), y(:,:), xi(:), eta(:), s(:), z(:)
        real :: dxi, deta, dz

!.... metrics

        real, allocatable :: m1(:,:),  m2(:,:),  n1(:,:),  n2(:,:),   &
                             m11(:,:), m12(:,:), m22(:,:),            &
                             n11(:,:), n12(:,:), n22(:,:)
        real :: m1l, m2l, bn1, bn2

!.... parameters

        real    :: Ma, Re, Pr, cv, time, gamma, tmp
        integer :: lstep, nx, ny, nz, ndof

        integer :: i, j, k, idof, itmp
        integer :: nx2, ny2, nz2

        character(80) file 
        
        integer :: optx=-1, opty=-1
        logical :: xper=.false., yper=.false.
        logical :: lsym=.true. , rsym=.false., bsym=.false., tsym=.false.
!============================================================================!

!.... read the restart file

        file = 'output.R.1'
        open(unit=10, file=file, form='unformatted')
        read(10) lstep, time, nx, ny, nz, ndof, Re, Ma, Pr, gamma, cv
        allocate( v(ny,nx,ndof) )
        read(10) v
        close(10)

!.... read in the grid file

        allocate( x(ny,nx), y(ny,nx), xi(nx), eta(ny), s(nx) )
        open(unit=10,file='grid.dat',form='unformatted')
        read(10) nx2, ny2, nz2
        if (nx2.ne.nx .or. ny2.ne.ny .or. nz2.ne.nz) then
          write(*,*) 'Grid and data file dimensions do not match'
          call exit(1)
        end if
        read(10) (((x(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((y(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((   tmp, i = 1, nx), j = 1, ny), k = 1, nz)
        close(10)

!.... make the xi grid
        
        dxi = one / float(nx-1)
        do i = 1, nx
          xi(i) = real(i-1) * dxi
        end do

!.... make the eta grid

        deta = one / float(ny-1)
        do j = 1, ny
          eta(j) = float(j-1) * deta
        end do

!.... allocate storage for metrics

        allocate (m1(ny,nx),  m2(ny,nx),  n1(ny,nx),  n2(ny,nx), &
                  m11(ny,nx), m12(ny,nx), m22(ny,nx),            &
                  n11(ny,nx), n12(ny,nx), n22(ny,nx) )

!.... read in the metric file

        open (unit=10,file='metric.dat',form='unformatted', &
              status='old')
        read(10) m1, m2, n1, n2, m11, m12, m22, n11, n12, n22
        close(10)
        
!.... Compute first derivatives of field in the mapped space

        allocate( g1v(ny,nx,ndof), g2v(ny,nx,ndof) )
        call grad(ndof, nx, ny, v, g1v, g2v, dxi, deta, optx, opty, &
                  xper, yper, lsym, rsym, bsym, tsym)

!.... write out boundary values

        open(20,file='top.pot',form='formatted')
        do i = 1, nx
          write(20,"(5(1pe21.14,1x))") (v(ny,i,idof), idof = 1, ndof)
        end do
        close(20)

        open(20,file='right.pot',form='formatted')
        do j = 1, ny
          write(20,"(5(1pe21.14,1x))") (v(j,nx,idof), idof = 1, ndof)
        end do
        close(20)

        open(20,file='pg.pot',form='unformatted')
        write(20) (((v(j,i,idof), j=1, ny), i=1, nx), idof = 2, ndof)
        write(20) ((g1v(j,i,5), j=1, ny), i=1, nx)
        close(20)

        call exit(0)
        end

!=============================================================================!
        module stencil
!
!  Finite Difference coefficient of various orders
!
!=============================================================================!
!.... First derivatives
!=============================================================================!

!.... fourth order central difference ( 1 2 x 4 5 )

        real, parameter :: ga1 =  8.333333333333333333333E-02
        real, parameter :: ga2 = -6.666666666666666666667E-01
        real, parameter :: ga3 =  6.666666666666666666667E-01
        real, parameter :: ga4 = -8.333333333333333333333E-02

!.... fourth order one-sided ( x 2 3 4 5 ) 

        real, parameter :: gc1 = -2.083333333333333333333E+00
        real, parameter :: gc2 =  4.000000000000000000000E+00
        real, parameter :: gc3 = -3.000000000000000000000E+00
        real, parameter :: gc4 =  1.333333333333333333333E+00
        real, parameter :: gc5 = -2.500000000000000000000E-01

!.... fourth order biased difference ( 1 x 2 3 4 5 )

        real, parameter :: gb1 = -2.500000000000000000000E-01
        real, parameter :: gb2 = -8.333333333333333333333E-01
        real, parameter :: gb3 =  1.500000000000000000000E+00
        real, parameter :: gb4 = -5.000000000000000000000E-01
        real, parameter :: gb5 =  8.333333333333333333333E-02

!.... sixth order one-sided ( x 2 3 4 5 6 7 )

        real, parameter :: ge1 = -2.450000000000000000000E+00
        real, parameter :: ge2 =  6.000000000000000000000E+00
        real, parameter :: ge3 = -7.500000000000000000000E+00
        real, parameter :: ge4 =  6.666666666666666666667E+00
        real, parameter :: ge5 = -3.750000000000000000000E+00
        real, parameter :: ge6 =  1.200000000000000000000E+00
        real, parameter :: ge7 = -1.666666666666666666667E-01

!.... sixth order one-pt biased ( 1 x 3 4 5 6 7 ) 

        real, parameter :: gf1 = -1.666666666666666666667E-01
        real, parameter :: gf2 = -1.283333333333333333333E+00
        real, parameter :: gf3 =  2.500000000000000000000E+00 
        real, parameter :: gf4 = -1.666666666666666666667E+00
        real, parameter :: gf5 =  8.333333333333333333333E-01
        real, parameter :: gf6 = -2.500000000000000000000E-01
        real, parameter :: gf7 =  3.333333333333333333333E-02

!==============================================================================
!.... The following stencils are Carpenter's stable and third
!.... order accurate boundary treatment for the explicit fourth
!.... order interior scheme.
!==============================================================================

!.... third order one-sided [Carpenter] ( x 2 3 4 5 6) 
  
        real, parameter :: gg1 = -1.8760320556207377229
        real, parameter :: gg2 =  3.185383225577892867
        real, parameter :: gg3 = -1.8145456794375275725
        real, parameter :: gg4 =  0.5916582410526027442
        real, parameter :: gg5 = -0.10105206800050562464
        real, parameter :: gg6 =  0.014588336428275308769

!.... third order biased [Carpenter] ( 1 x 3 4 5 ) 
  
        real, parameter :: gh1 = -0.38425423267792540204
        real, parameter :: gh2 = -0.29063894776734868107
        real, parameter :: gh3 =  0.6717647845153154114
        real, parameter :: gh4 =  0.07108165983739987271
        real, parameter :: gh5 = -0.07363071876172424507
        real, parameter :: gh6 =  0.00567745485428304409

!.... third order biased [Carpenter] ( 1 2 x 4 5 ) 
  
        real, parameter :: gi1 =  0.18288527868682620658
        real, parameter :: gi2 = -1.0800147541745551643
        real, parameter :: gi3 =  0.6578728964966252582
        real, parameter :: gi4 =  0.17761704868919314564
        real, parameter :: gi5 =  0.0767798363958275586
        real, parameter :: gi6 = -0.015140306093917004671

!.... third order biased [Carpenter] ( 1 2 3 x 5 ) 
  
        real, parameter :: gj1 = -0.03418371033652918578
        real, parameter :: gj2 =  0.22482902574010312173
        real, parameter :: gj3 = -0.8908123329284539625
        real, parameter :: gj4 =  0.16529994771003501478
        real, parameter :: gj5 =  0.6134395520875252998
        real, parameter :: gj6 = -0.07857248227268028806

!=============================================================================!
!.... Second derivatives
!=============================================================================!

!.... fourth order central difference

        real, parameter :: da1 = -8.333333333333333333333E-02
        real, parameter :: da2 =  1.333333333333333333333E+00
        real, parameter :: da3 = -2.500000000000000000000E+00
        real, parameter :: da4 =  1.333333333333333333333E+00
        real, parameter :: da5 = -8.333333333333333333333E-02
        
!.... third order biased difference

        real, parameter :: db1 =  9.166666666666666666667E-01
        real, parameter :: db2 = -1.666666666666666666667E+00
        real, parameter :: db3 =  5.000000000000000000000E-01
        real, parameter :: db4 =  3.333333333333333333333E-01
        real, parameter :: db5 = -8.333333333333333333333E-02
        
!.... fourth order one-sided (assumes f'=0) [not smooth]

        real, parameter :: dc1 = -5.763888888888888888889E+00
        real, parameter :: dc2 =  8.000000000000000000000E+00
        real, parameter :: dc3 = -3.000000000000000000000E+00
        real, parameter :: dc4 =  8.888888888888888888889E-01
        real, parameter :: dc5 = -1.250000000000000000000E-01
        
!.... third order one-sided

        real, parameter :: dd1 =  2.916666666666666666667E+00
        real, parameter :: dd2 = -8.666666666666666666667E+00
        real, parameter :: dd3 =  9.500000000000000000000E+00
        real, parameter :: dd4 = -4.666666666666666666667E+00
        real, parameter :: dd5 =  9.166666666666666666667E-01
        
        end module stencil

!=============================================================================!
        subroutine grad( ndof, nx, ny, v, g1v, g2v, dx, dy, optx, opty,&
                         xper, yper, lsym, rsym, bsym, tsym)
!
!  Take the gradient of a 2-D field.
!  Updated to sixth-order accurate differencing on the interior with
!  the option for optimized fourth-order differencing.
!
!  Revised: 6-28-95
!
!=============================================================================!
        use stencil
        implicit none
        
        integer :: ndof, nx, ny, optx, opty
        logical :: xper, yper
        logical :: lsym, rsym, bsym, tsym
        real    :: v(ny,nx,ndof), g1v(ny,nx,ndof), g2v(ny,nx,ndof)
        real    :: dx, dy
        
        real, parameter :: one = 1.0, pt5 = 0.5
        real dxinv, dyinv
        real a, b, c, w
        real gx1, gx2, gx3, gx4, gx5, gx6
        real gy1, gy2, gy3, gy4, gy5, gy6
        
        integer :: idof, isign
!=============================================================================!

        dxinv  = one / dx
        dyinv  = one / dy

!.... seven point stencil in x

        if (optx.eq.0) then
          c = 1.0 / 60.0
        else if (optx.eq.-1) then
          c = 0.0
        else
          w = 2.0 * 3.1415926535897932385e+0 / 12.0
          c = (w/2.0 - 2.0*Sin(w)/3.0 + Sin(2.0*w)/12.0) / &
              (5.0*Sin(w) - 4.0*Sin(2.0*w) + Sin(3.0*w))
        end if
        
        a = (2.0 + 15.0 * c) / 3.0
        b = -(1.0 + 48.0 * c) / 12.0

        gx1 =  -c * dxinv
        gx2 =  -b * dxinv
        gx3 =  -a * dxinv
        gx4 =   a * dxinv
        gx5 =   b * dxinv
        gx6 =   c * dxinv

!.... seven point stencil in y

        if (opty.eq.0) then
          c = 1.0 / 60.0
        else if (opty.eq.-1) then
          c = 0.0
        else
          w = 2.0 * 3.1415926535897932385e+0 / 12.0
          c = (w/2.0 - 2.0*Sin(w)/3.0 + Sin(2.0*w)/12.0) / &
              (5.0*Sin(w) - 4.0*Sin(2.0*w) + Sin(3.0*w))
        end if
        
        a = (2.0 + 15.0 * c) / 3.0
        b = -(1.0 + 48.0 * c) / 12.0

        gy1 =  -c * dyinv
        gy2 =  -b * dyinv
        gy3 =  -a * dyinv
        gy4 =   a * dyinv
        gy5 =   b * dyinv
        gy6 =   c * dyinv

!=============================================================================!
!.... compute the gradient in x
!=============================================================================!

        if (xper) then

          g1v(:,1,:)        = ( gx1 * v(:,nx-3,:)       + &
                                gx2 * v(:,nx-2,:)       + &
                                gx3 * v(:,nx-1,:)       + &
                                gx4 * v(:,2,:)          + &
                                gx5 * v(:,3,:)          + &
                                gx6 * v(:,4,:)  ) 
  
          g1v(:,2,:)        = ( gx1 * v(:,nx-2,:)       + &
                                gx2 * v(:,nx-1,:)       + &
                                gx3 * v(:,1,:)          + &
                                gx4 * v(:,3,:)          + &
                                gx5 * v(:,4,:)          + &
                                gx6 * v(:,5,:)  ) 
  
          g1v(:,3,:)        = ( gx1 * v(:,nx-1,:)       + &
                                gx2 * v(:,1,:)          + &
                                gx3 * v(:,2,:)          + &
                                gx4 * v(:,4,:)          + &
                                gx5 * v(:,5,:)          + &
                                gx6 * v(:,6,:)  ) 
  
          g1v(:,nx-2,:)     = ( gx1 * v(:,nx-5,:)       + &
                                gx2 * v(:,nx-4,:)       + &
                                gx3 * v(:,nx-3,:)       + &
                                gx4 * v(:,nx-1,:)       + &
                                gx5 * v(:,1,:)          + &
                                gx6 * v(:,2,:)  ) 
  
          g1v(:,nx-1,:)     = ( gx1 * v(:,nx-4,:)       + &
                                gx2 * v(:,nx-3,:)       + &
                                gx3 * v(:,nx-2,:)       + &
                                gx4 * v(:,1,:)          + &
                                gx5 * v(:,2,:)          + &
                                gx6 * v(:,3,:)  ) 
  
          g1v(:,nx,:) = g1v(:,1,:)
          
        else
        
!         g1v(:,1,:)      = ( v(:,2,:) - v(:,1,:) ) * dxinv
          g1v(:,1,:)      = ( gc1 * v(:,1,:)  + &
                              gc2 * v(:,2,:)  + &
                              gc3 * v(:,3,:)  + &
                              gc4 * v(:,4,:)  + &
                              gc5 * v(:,5,:)  ) * dxinv
!         g1v(:,1,:)      = ( ge1 * v(:,1,:)  + &
!                             ge2 * v(:,2,:)  + &
!                             ge3 * v(:,3,:)  + &
!                             ge4 * v(:,4,:)  + &
!                             ge5 * v(:,5,:)  + &
!                             ge6 * v(:,6,:)  + &
!                             ge7 * v(:,7,:)  ) * dxinv
  
!         g1v(:,2,:)      = ( v(:,3,:) - v(:,1,:) ) * pt5 * dxinv
          g1v(:,2,:)      = ( gb1 * v(:,1,:)  + &
                              gb2 * v(:,2,:)  + &
                              gb3 * v(:,3,:)  + &
                              gb4 * v(:,4,:)  + &
                              gb5 * v(:,5,:)  ) * dxinv
!         g1v(:,2,:)      = ( gf1 * v(:,1,:)  + &
!                             gf2 * v(:,2,:)  + &
!                             gf3 * v(:,3,:)  + &
!                             gf4 * v(:,4,:)  + &
!                             gf5 * v(:,5,:)  + &
!                             gf6 * v(:,6,:)  + &
!                             gf7 * v(:,7,:)  ) * dxinv

          g1v(:,3,:)      = ( ga1 * v(:,1,:)  + &
                              ga2 * v(:,2,:)  + &
                              ga3 * v(:,4,:)  + &
                              ga4 * v(:,5,:)  ) * dxinv
          
          g1v(:,nx-2,:)   = ( ga1 * v(:,nx-4,:)  + &
                              ga2 * v(:,nx-3,:)  + &
                              ga3 * v(:,nx-1,:)  + &
                              ga4 * v(:,nx  ,:)  ) * dxinv
  
          g1v(:,nx-1,:)  = -( gb1 * v(:,nx  ,:)  + &
                              gb2 * v(:,nx-1,:)  + &
                              gb3 * v(:,nx-2,:)  + &
                              gb4 * v(:,nx-3,:)  + &
                              gb5 * v(:,nx-4,:)  ) * dxinv
  
          g1v(:,nx,:)    = -( gc1 * v(:,nx  ,:)  + &
                              gc2 * v(:,nx-1,:)  + &
                              gc3 * v(:,nx-2,:)  + &
                              gc4 * v(:,nx-3,:)  + &
                              gc5 * v(:,nx-4,:)  ) * dxinv

        end if

!.... implement the symmetric conditions

        if (lsym) then
          do idof = 1, ndof
            if (idof .eq. 3) then
              isign = -one
            else
              isign = one
            end if
  
            g1v(:,1,idof)     = ( isign * gx1 * v(:,4,idof)     + &
                                  isign * gx2 * v(:,3,idof)     + &
                                  isign * gx3 * v(:,2,idof)     + &
                                  gx4 * v(:,2,idof)             + &
                                  gx5 * v(:,3,idof)             + &
                                  gx6 * v(:,4,idof)             ) 
    
            g1v(:,2,idof)     = ( isign * gx1 * v(:,3,idof)     + &
                                  isign * gx2 * v(:,2,idof)     + &
                                  gx3 * v(:,1,idof)             + &
                                  gx4 * v(:,3,idof)             + &
                                  gx5 * v(:,4,idof)             + &
                                  gx6 * v(:,5,idof)             ) 
    
            g1v(:,3,idof)     = ( isign * gx1 * v(:,2,idof)     + &
                                  gx2 * v(:,1,idof)             + &
                                  gx3 * v(:,2,idof)             + &
                                  gx4 * v(:,4,idof)             + &
                                  gx5 * v(:,5,idof)             + &
                                  gx6 * v(:,6,idof)             ) 
          end do
        end if
        
!.... interior

        g1v(:,4:nx-3,:)  = ( gx1 * v(:,1:nx-6,:)        + &
                             gx2 * v(:,2:nx-5,:)        + &
                             gx3 * v(:,3:nx-4,:)        + &
                             gx4 * v(:,5:nx-2,:)        + &
                             gx5 * v(:,6:nx-1,:)        + &
                             gx6 * v(:,7:nx  ,:)        ) 

!=============================================================================!
!.... compute the gradient in y
!=============================================================================!

        if (yper) then
        
          g2v(1,:,:)       = ( gy1 * v(ny-3,:,:)        + &
                                gy2 * v(ny-2,:,:)       + &
                                gy3 * v(ny-1,:,:)       + &
                                gy4 * v(2,:,:)          + &
                                gy5 * v(3,:,:)          + &
                                gy6 * v(4,:,:)          ) 
  
          g2v(2,:,:)       = ( gy1 * v(ny-2,:,:)        + &
                                gy2 * v(ny-1,:,:)       + &
                                gy3 * v(1,:,:)          + &
                                gy4 * v(3,:,:)          + &
                                gy5 * v(4,:,:)          + &
                                gy6 * v(5,:,:)  ) 
  
          g2v(3,:,:)       = ( gy1 * v(ny-1,:,:)        + &
                                gy2 * v(1,:,:)          + &
                                gy3 * v(2,:,:)          + &
                                gy4 * v(4,:,:)          + &
                                gy5 * v(5,:,:)          + &
                                gy6 * v(6,:,:)  ) 
  
          g2v(ny-2,:,:)    = ( gy1 * v(ny-5,:,:)        + &
                                gy2 * v(ny-4,:,:)       + &
                                gy3 * v(ny-3,:,:)       + &
                                gy4 * v(ny-1,:,:)       + &
                                gy5 * v(1,:,:)          + &
                                gy6 * v(2,:,:)  ) 
  
          g2v(ny-1,:,:)    = ( gy1 * v(ny-4,:,:)        + &
                                gy2 * v(ny-3,:,:)       + &
                                gy3 * v(ny-2,:,:)       + &
                                gy4 * v(1,:,:)          + &
                                gy5 * v(2,:,:)          + &
                                gy6 * v(3,:,:)  ) 
  
          g2v(ny,:,:) = g2v(1,:,:)

        else

          g2v(1,:,:)     = ( gc1 * v(1,:,:)     + &
                             gc2 * v(2,:,:)     + &
                             gc3 * v(3,:,:)     + &
                             gc4 * v(4,:,:)     + &
                             gc5 * v(5,:,:)  ) * dyinv

          g2v(2,:,:)     = ( gb1 * v(1,:,:)     + &
                             gb2 * v(2,:,:)     + &
                             gb3 * v(3,:,:)     + &
                             gb4 * v(4,:,:)     + &
                             gb5 * v(5,:,:)  ) * dyinv
  
          g2v(3,:,:)     = ( ga1 * v(1,:,:)     + &
                             ga2 * v(2,:,:)     + &
                             ga3 * v(4,:,:)     + &
                             ga4 * v(5,:,:)  ) * dyinv
  
          g2v(ny-2,:,:)  = ( ga1 * v(ny-4,:,:)  + &
                             ga2 * v(ny-3,:,:)  + &
                             ga3 * v(ny-1,:,:)  + &
                             ga4 * v(ny  ,:,:)  ) * dyinv
  
          g2v(ny-1,:,:) = -( gb1 * v(ny  ,:,:)  + &
                             gb2 * v(ny-1,:,:)  + &
                             gb3 * v(ny-2,:,:)   + &
                             gb4 * v(ny-3,:,:)  + &
                             gb5 * v(ny-4,:,:)  ) * dyinv
  
          g2v(ny,:,:)   = -( gc1 * v(ny  ,:,:)  + &
                             gc2 * v(ny-1,:,:)  + &
                             gc3 * v(ny-2,:,:)  + &
                             gc4 * v(ny-3,:,:)  + &
                             gc5 * v(ny-4,:,:)  ) * dyinv
                            
        end if
        
!.... interior

        g2v(4:ny-3,:,:) = ( gy1 * v(1:ny-6,:,:) + &
                            gy2 * v(2:ny-5,:,:) + &
                            gy3 * v(3:ny-4,:,:) + &
                            gy4 * v(5:ny-2,:,:) + &
                            gy5 * v(6:ny-1,:,:) + &
                            gy6 * v(7:ny  ,:,:) ) 
        
        return
        end
