c---------------------------------------------------------------------
        subroutine splcrv( npts, s, x, y, bsx, bsy, nknot, 
     &                     xknot, yknot, korder, ncoef )
c---------------------------------------------------------------------
        implicit real (a-h,o-z)

        dimension s(npts), x(npts), y(npts), bsx(ncoef), bsy(ncoef)
        dimension xknot(nknot), yknot(nknot)
c
c.... local storage
c
        parameter (mpts = 1024)
        dimension xguess(mpts), weight(mpts)
c
c.... B-spline constraints
c
        dimension xval(2), nhard(2), ider(2), itype(2), bl(2), bu(2)

        if (npts .gt. mpts .or. ncoef .gt. mpts) then
          write(*,*) 'ERROR: increase mpts in splcrv'
          call exit(1)
        end if
c---------------------------------------------------------------------
c     Regular B-spline
c---------------------------------------------------------------------
        if (ncoef .eq. npts) then

c.... compute the knot locations
c
        call BSNAK(npts, s, korder, xknot)
        call BSNAK(npts, s, korder, yknot)
c
c       call BSOPK(npts, s, korder, xknot)
c       call BSOPK(npts, s, korder, yknot)
c
c.... compute the B-spline coefficients
c
        call BSINT(npts, s, x, korder, xknot, bsx)
        call BSINT(npts, s, y, korder, yknot, bsy)
c
c.... now that I have my knot points, recompute with contraints
c
c.... set the weights
c
c       do i = 1, npts
c         weight(i) = 1.0
c       end do
c        weight(1) = 1.0e5
c       weight(npts) = 1.0e5
c
c       nxval    = 1
c       xval(i)  = s(1)
c       nhard(1) = 0
c       ider(1)  = 1
c       itype(1) = 1
c       bl(1)    = 0.0
c       bu(1)    = 0.0
c       call CONFT(npts, s, x, weight, nxval, xval, nhard, ider, itype, 
c     &             bl, bu, korder, xknot, ncoef, bsx)
c
c       nxval    = 1
c       xval(i)  = s(1)
c       nhard(1) = 0
c       ider(1)  = 1
c       itype(1) = 1
c       bl(1)    = 1.0
c       bu(1)    = 0.0
c       call CONFT(npts, s, y, weight, nxval, xval, nhard, ider, itype, 
c     &             bl, bu, korder, yknot, ncoef, bsy)

        else

        write(*,*) 'Least-squares B-spline not implemented...'
        call exit(1)
c---------------------------------------------------------------------
c     Least-squares B-spline
c---------------------------------------------------------------------
c
c.... set the weights
c
        do i = 1, npts
          weight(i) = 1.0
        end do
        weight(1) = 1.0e5
        weight(npts) = 1.0e5
c
c.... make a guess at the note locations
c
        do i = 1, ncoef - korder + 2
           xguess(i+korder-1) = (s(npts)+0.0001)*float(i-1)/
     &                          float(ncoef-korder+1)
        end do
        do i = 1, korder - 1
           xguess(i) = xguess(korder)
           xguess(i+ncoef+1) = xguess(ncoef+1)
        end do
c
c       call BSVLS(npts, s, x, weight, korder, ncoef, xguess, xknot, bsx, ssq)
        write(*,"('Root mean square error in x = ',1pe13.6)") ssq
c       call BSVLS(npts, s, y, weight, korder, ncoef, xguess, yknot, bsy, ssq)
        write(*,"('Root mean square error in y = ',1pe13.6)") ssq
c
c.... now that I have my knot points, recompute with contraints
c
        nxval    = 1
        xval(i)  = s(1)
        nhard(1) = 0
        ider(1)  = 1
        itype(1) = 1
        bl(1)    = 0.0
        bu(1)    = 0.0
c       call CONFT(npts, s, x, weight, nxval, xval, nhard, ider, itype, 
c     &             bl, bu, korder, xknot, ncoef, bsx)

        nxval    = 1
        xval(i)  = s(1)
        nhard(1) = 0
        ider(1)  = 1
        itype(1) = 1
        bl(1)    = 1.0
        bu(1)    = 0.0
c       call CONFT(npts, s, y, weight, nxval, xval, nhard, ider, itype, 
c     &             bl, bu, korder, yknot, ncoef, bsy)

        end if
c
c.... check the spline
c
        if (.false.) then
        itime = 5
        do i = 1, npts-1
           i0 = (i-1)*itime+1
           ds = (s(i+1)-s(i))/float(itime)
           do j = 0, itime-1
              sint = s(i) + float(j) * ds
              xint   = BSDER( 0, sint, korder, xknot, ncoef, bsx )
              yint   = BSDER( 0, sint, korder, yknot, ncoef, bsy )
              dxds   = BSDER( 1, sint, korder, xknot, ncoef, bsx )
              d2xds2 = BSDER( 2, sint, korder, xknot, ncoef, bsx )
              d3xds3 = BSDER( 3, sint, korder, xknot, ncoef, bsx )
              write(16,10) sint, xint, yint, dxds, d2xds2, d3xds3
 10           format(8(1pe20.13,1x))
           end do
        end do
        sint = s(npts)
        xint   = BSDER( 0, sint, korder, xknot, ncoef, bsx )
        yint   = BSDER( 0, sint, korder, yknot, ncoef, bsy )
        dxds   = BSDER( 1, sint, korder, xknot, ncoef, bsx )
        d2xds2 = BSDER( 2, sint, korder, xknot, ncoef, bsx )
        d3xds3 = BSDER( 3, sint, korder, xknot, ncoef, bsx )
        write(16,10) sint, xint, yint, dxds, d2xds2, d3xds3
        close(32)
        stop
        end if
c
        return
        end
C***********************************************************************
      SUBROUTINE SPLINE(N,X,Y,FDP)      
C***********************************************************************
C-----THIS SUBROUTINE COMPUTES THE SECOND DERIVATIVES NEEDED 
C-----IN CUBIC SPLINE INTERPOLATION.  THE INPUT DATA ARE:    
C-----N = NUMBER OF DATA POINTS          
C-----X = ARRAY CONTAINING THE VALUES OF THE INDEPENDENT VARIABLE      
C-----    (ASSUMED TO BE IN ASCENDING ORDER)       
C-----Y = ARRAY CONTAINING THE VALUES OF THE FUNCTION AT THE 
C-----    DATA POINTS GIVEN IN THE X ARRAY         
C-----THE OUTPUT IS THE ARRAY FDP WHICH CONTAINS THE SECOND  
C-----DERIVATIVES OF THE INTERPOLATING CUBIC SPLINE.         
      IMPLICIT REAL (A-H,O-Z)
      PARAMETER (NMAX=1000)
      DIMENSION X(N),Y(N),FDP(N)
      DIMENSION A(NMAX),B(NMAX),C(NMAX),R(NMAX)
C-----COMPUTE THE COEFFICIENTS AND THE RHS OF THE EQUATIONS. 
C-----THIS ROUTINE USES THE CANTILEVER CONDITION.  THE PARAMETER       
C-----ALAMDA (LAMBDA) IS SET TO 1. BUT THIS CAN BE USER-MODIFIED.      
C-----A,B,C ARE THE THREE DIAGONALS OF THE TRIDIAGONAL SYSTEM;         
C-----R IS THE RIGHT HAND SIDE.  THESE ARE NOW ASSEMBLED.    
      IF (N .GT. NMAX) THEN
         WRITE (*,*) 'INCREASE NMAX IN SPLINE...'
         CALL EXIT(1)
      END IF
      ALAMDA = 0.    
      NM2 = N - 2    
      NM1 = N - 1    
      C(1) = X(2) - X(1)       
      DO 1 I=2,NM1   
      C(I) = X(I+1) - X(I)     
      A(I) = C(I-1)  
      B(I) = 2.*(A(I) + C(I))  
      R(I) = 6.*((Y(I+1) - Y(I))/C(I) - (Y(I) - Y(I-1))/C(I-1))        
    1 CONTINUE       
      B(2) = B(2) + ALAMDA * C(1)        
      B(NM1) = B(NM1) + ALAMDA * C(NM1)  
C-----AT THIS POINT WE COULD CALL A TRIDIAGONAL SOLVER SUBROUTINE      
C-----BUT THE NOTATION IS CLUMSY SO WE WILL SOLVE DIRECTLY.  THE       
C-----NEXT SECTION SOLVES THE SYSTEM WE HAVE JUST SET UP.    
      DO 2 I=3,NM1   
      T = A(I)/B(I-1)          
      B(I) = B(I) - T * C(I-1) 
      R(I) = R(I) - T * R(I-1) 
    2 CONTINUE       
      FDP(NM1) = R(NM1)/B(NM1) 
      DO 3 I=2,NM2   
      NMI = N - I    
      FDP(NMI) = (R(NMI) - C(NMI)*FDP(NMI+1))/B(NMI)         
    3 CONTINUE       
      FDP(1) = ALAMDA * FDP(2) 
      FDP(N) = ALAMDA * FDP(NM1)         
C-----WE NOW HAVE THE DESIRED DERIVATIVES SO WE RETURN TO THE          
C-----MAIN PROGRAM.  
      RETURN         
      END  
C***********************************************************************
      SUBROUTINE SPEVAL (N,X,Y,FDP,XX,F) 
C***********************************************************************
C-----THIS SUBROUTINE EVALUATES THE CUBIC SPLINE GIVEN       
C-----THE DERIVATIVE COMPUTED BY SUBROUTINE SPLINE.          
C-----THE INPUT PARAMETERS N,X,Y,FDP HAVE THE SAME 
C-----MEANING AS IN SPLINE.    
C-----XX = VALUE OF INDEPENDENT VARIABLE FOR WHICH 
C-----     AN INTERPOLATED VALUE IS REQUESTED      
C-----F =  THE INTERPOLATED RESULT       
      IMPLICIT REAL (A-H,O-Z)
      DIMENSION X(N),Y(N),FDP(N)      
C-----THE FIRST JOB IS TO FIND THE PROPER INTERVAL.          
      NM1 = N - 1    
      DO 1 I=1,NM1   
      IF (XX.LE.X(I+1)) GO TO 10         
    1 CONTINUE       
C-----NOW EVALUATE THE CUBIC   
   10 DXM = XX - X(I)          
      DXP = X(I+1) - XX        
      DEL = X(I+1) - X(I)      
      F = FDP(I)*DXP*(DXP*DXP/DEL - DEL)/6.        
     1   +FDP(I+1)*DXM*(DXM*DXM/DEL - DEL)/6.     
     2   +Y(I)*DXP/DEL + Y(I+1)*DXM/DEL 
      RETURN        
      END 
