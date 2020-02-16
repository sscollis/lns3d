!=============================================================================!
        subroutine genbump( g1v, g2v, g11v, g12v, g22v )
!
!  generate the wall bump
!
!=============================================================================!
        use global
        use bump
        implicit none
        
        real :: g1v(ny,nx,ndof), g2v(ny,nx,ndof)
        real :: g11v(ny,nx,ndof), g12v(ny,nx,ndof), g22v(ny,nx,ndof)
        character*80 ::code='GenBump$'

        integer :: i
        real :: s(nx)
        
        real :: bmu=0.7, bsigma=0.05
        integer :: btype
!=============================================================================!
        allocate( hw(nx) )
            
        call error('genbump$','Must update indices$')

!.... read the mean and standard deviation from the bump.dat file

        open(30,file='bump.inp',err=1000)
        read(30,*,err=1010,end=1010) bmu, bsigma, btype
        close(30)

!.... specify the bump type

        if (btype .eq. 0) then
          s = x(1,:)
        else if (btype .eq. 1) then
          s = Sqrt(x(1,:) + two*x(1,:)**2)/Sqrt(two) + &
              Log(one + four*x(1,:) + two**(onept5)*Sqrt(x(1,:) + &
              two*x(1,:)**2)) / four
        else
          call error(code,'Illegal value of btype in bump.inp$')
        end if

!.... Normalized Gaussian bump

        hw(:) = exp( -pt5 * ( ( s - bmu ) / bsigma )**2 )
        
        allocate( u1p(nx), u2p(nx), u3p(nx), tp(nx), tpp(nx) )
        allocate( u1w(nx), u2w(nx), u3w(nx), tw(nx), twp(nx) )

!.... note there is a negative sign since the normal points outside the domain

        u1p = -( bnb(:,1) * g1v(1,:,2) + bnb(:,2) * g2v(1,:,2) )
        u2p = -( bnb(:,1) * g1v(1,:,3) + bnb(:,2) * g2v(1,:,3) )
        u3p = -( bnb(:,1) * g1v(1,:,4) + bnb(:,2) * g2v(1,:,4) )
        tp  = -( bnb(:,1) * g1v(1,:,5) + bnb(:,2) * g2v(1,:,5) )
        tpp = bnb(:,1) * (bnb(:,1) * g11v(1,:,5) + bnb(:,2) * g12v(1,:,5)) + &
              bnb(:,2) * (bnb(:,1) * g12v(1,:,5) + bnb(:,2) * g22v(1,:,5))

        u1w = -u1p * hw
        u2w = -u2p * hw
        u3w = -u3p * hw
        tw  = -tp  * hw
        twp = -tpp * hw

!.... save the bump profile 

        open(30,file='bump.dat')
        do i = 1, nx
          write(30,10) s(i), hw(i), real(u1w(i)), real(u2w(i)), &
                       real(u3w(i)), real(tw(i)), real(twp(i))
        end do
        close(30)
        
        return
  10    format( 8(1pe16.7e4,1x) )
1000    call error(code,'Error opening bump.inp$')
1010    call error(code,'Error reading bump.inp$')
        end
