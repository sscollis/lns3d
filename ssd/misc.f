c------------------------------------------------------------------------
c       subroutine scfft( )
c
c       return
c       end
c------------------------------------------------------------------------
        subroutine ccfft( )
c 
        write(*,*) 'Call to ccfft'
        stop
        end
c------------------------------------------------------------------------
        subroutine exit( i )
c
c.... A simple hack to mimic the tremain routine on the Cray
c
        stop
c       
        return
        end
c------------------------------------------------------------------------
        subroutine tremain( rtime )
c
c.... A simple hack to mimic the tremain routine on the Cray
c
        rtime = 1.0e10
c       
        return
        end
c------------------------------------------------------------------------
        function ishell( command )
c
c.... A simple hack to mimic the ishell routine on the Cray
c
        integer ishell
        character*256 command
        
        ishell = 1
        
        return
        end
#ifndef __GFORTRAN__
c------------------------------------------------------------------------
      REAL  FUNCTION SECOND( )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992 
*
*  Purpose
*  =======
*
*  SECOND returns the user time for a process in seconds.
*  This version gets the time from the system function ETIME.
*
*     .. Local Scalars ..
      REAL*4             T1
*     ..
*     .. Local Arrays ..
      REAL*4             TARRAY( 2 )
*     ..
*     .. External Functions ..
      REAL*4             ETIME
      EXTERNAL           ETIME
*     ..
*     .. Executable Statements ..
*
      T1 = ETIME( TARRAY )
      SECOND = TARRAY( 1 ) + TARRAY( 2 )
      RETURN
*
*     End of SECOND
*
      END
#endif
c------------------------------------------------------------------------
c       subroutine flush(iunit)
c
c.... A simple hack to mimic the flush routine on the Cray
c
c       integer iunit
c       
c       return
c       end
c------------------------------------------------------------------------
        subroutine scfft(i1,i2,r1,r2,r3,r4,r5,i3)
c
c.... A simple hack to mimic the FFT routine on the Cray
c
        write (*,*) 'SCFFT not supported on SGI machines...'
        call exit(1)
        ierr = 0
        
        return
        end
c------------------------------------------------------------------------
c
c.... SSD rooutines on the Iris
c
c------------------------------------------------------------------------
        subroutine opendr( iunit, index, length, it, ierr )
        
        write (*,*) 'SSD not supported on SGI machines...'
        call exit(1)
        
        return
        end
c------------------------------------------------------------------------
        subroutine readdr( iunit, ubuff, n, irec, irflag, is, ierr )
        
        ierr = 0
        
        return
        end
c------------------------------------------------------------------------
        subroutine writdr( iunit, ubuff, n, irec, irflag, is, ierr )
        
        ierr = 0
        
        return
        end
c------------------------------------------------------------------------
        subroutine waitdr( iunit, is, ierr )
        
        ierr = 0
        
        return
        end
c------------------------------------------------------------------------
        subroutine asyncdr( iunit, ierr )
        
        ierr = 0
        
        return
        end
c------------------------------------------------------------------------
        subroutine closdr( iunit, ierr )
        
        close(99)
        
        return
        end
c------------------------------------------------------------------------
