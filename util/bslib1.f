c
c.... B-Spline Routines from PPPack by c. de boor
c
c.... These are avalable on NetLib.org as public domain software
c
c.... Assembled by S. Collis
c
c
c.... Be sure to compile with -r8 otherwise you'll get roundoff errors for
c.... large nx
c
      subroutine splint ( tau, gtau, t, n, k, q, bcoef, iflag )
c  from  * a practical guide to splines *  by c. de boor    
c  calls bsplvb, banfac/slv
c
c   splint  produces the b-spline coeff.s  bcoef  of the spline of order
c   k  with knots  t(i), i=1,..., n + k , which takes on the value
c   gtau(i) at  tau(i), i=1,..., n .
c
c******  i n p u t  ******
c  tau.....array of length  n , containing data point abscissae.
c    a s s u m p t i o n . . .  tau  is strictly increasing
c  gtau.....corresponding array of length  n , containing data point or-
c        dinates
c  t.....knot sequence, of length  n+k
c  n.....number of data points and dimension of spline space  s(k,t)
c  k.....order of spline
c
c******  o u t p u t  ******
c  q.....array of size  (2*k-1)*n , containing the triangular factoriz-
c        ation of the coefficient matrix of the linear system for the b-
c        coefficients of the spline interpolant.
c           the b-coeffs for the interpolant of an additional data set
c        (tau(i),htau(i)), i=1,...,n  with the same data abscissae can
c        be obtained without going through all the calculations in this
c        routine, simply by loading  htau  into  bcoef  and then execut-
c        ing the    call banslv ( q, 2*k-1, n, k-1, k-1, bcoef )
c  bcoef.....the b-coefficients of the interpolant, of length  n
c  iflag.....an integer indicating success (= 1)  or failure (= 2)
c        the linear system to be solved is (theoretically) invertible if
c        and only if
c              t(i) .lt. tau(i) .lt. t(i+k),    all i.
c        violation of this condition is certain to lead to  iflag = 2 .
c
c******  m e t h o d  ******
c     the i-th equation of the linear system  a*bcoef = b  for the b-co-
c  effs of the interpolant enforces interpolation at  tau(i), i=1,...,n.
c  hence,  b(i) = gtau(i), all i, and  a  is a band matrix with  2k-1
c   bands (if it is invertible).
c     the matrix  a  is generated row by row and stored, diagonal by di-
c  agonal, in the  r o w s  of the array  q , with the main diagonal go-
c  ing into row  k .  see comments in the program below.
c     the banded system is then solved by a call to  banfac (which con-
c  structs the triangular factorization for  a  and stores it again in
c   q ), followed by a call to  banslv (which then obtains the solution
c   bcoef  by substitution).
c     banfac  does no pivoting, since the total positivity of the matrix
c  a  makes this unnecessary.
c
      integer iflag,k,n,   i,ilp1mx,j,jj,km1,kpkm2,left,lenq,np1
      real bcoef(n),gtau(n),q(1),t(1),tau(n),   taui
c     dimension q(2*k-1,n), t(n+k)
c  current fortran standard makes it impossible to specify precisely the
c  dimension of  q  and  t  without the introduction of otherwise super-
c  fluous additional arguments.
      np1 = n + 1
      km1 = k - 1
      kpkm2 = 2*km1
      left = k
c                zero out all entries of q
      lenq = n*(k+km1)
      do 5 i=1,lenq
    5    q(i) = 0.
c
c  ***   loop over i to construct the  n  interpolation equations
      do 30 i=1,n
         taui = tau(i)
         ilp1mx = min0(i+k,np1)
c        *** find  left  in the closed interval (i,i+k-1) such that
c                t(left) .le. tau(i) .lt. t(left+1)
c        matrix is singular if this is not possible
         left = max0(left,i)
         if (taui .lt. t(left))         go to 998
   15       if (taui .lt. t(left+1))    go to 16
            left = left + 1
            if (left .lt. ilp1mx)       go to 15
         left = left - 1
         if (taui .gt. t(left+1))       go to 998
c        *** the i-th equation enforces interpolation at taui, hence
c        a(i,j) = b(j,k,t)(taui), all j. only the  k  entries with  j =
c        left-k+1,...,left actually might be nonzero. these  k  numbers
c        are returned, in  bcoef (used for temp.storage here), by the
c        following
   16    call bsplvb ( t, k, 1, taui, left, bcoef )
c        we therefore want  bcoef(j) = b(left-k+j)(taui) to go into
c        a(i,left-k+j), i.e., into  q(i-(left+j)+2*k,(left+j)-k) since
c        a(i+j,j)  is to go into  q(i+k,j), all i,j,  if we consider  q
c        as a two-dim. array , with  2*k-1  rows (see comments in
c        banfac). in the present program, we treat  q  as an equivalent
c        one-dimensional array (because of fortran restrictions on
c        dimension statements) . we therefore want  bcoef(j) to go into
c        entry
c            i -(left+j) + 2*k + ((left+j) - k-1)*(2*k-1)
c                   =  i-left+1 + (left -k)*(2*k-1) + (2*k-2)*j
c        of  q .
         jj = i-left+1 + (left-k)*(k+km1)
         do 30 j=1,k
            jj = jj+kpkm2
   30       q(jj) = bcoef(j)
c
c     ***obtain factorization of  a  , stored again in  q.
      call banfac ( q, k+km1, n, km1, km1, iflag )
                                        go to (40,999), iflag
c     *** solve  a*bcoef = gtau  by backsubstitution
   40 do 41 i=1,n
   41    bcoef(i) = gtau(i)
      call banslv ( q, k+km1, n, km1, km1, bcoef )
                                        return
  998 iflag = 2
  999 print 699
  699 format(41h linear system in  splint  not invertible)
                                        return
      end

      subroutine bsplvb ( t, jhigh, index, x, left, biatx )
c  from  * a practical guide to splines *  by c. de boor    
c  calculates the value of all possibly nonzero b-splines at  x  of order
c
c               jout  =  max( jhigh , (j+1)*(index-1) )
c
c  with knot sequence  t .
c
c******  i n p u t  ******
c  t.....knot sequence, of length  left + jout  , assumed to be nonde-
c        creasing.  a s s u m p t i o n . . . .
c                       t(left)  .lt.  t(left + 1)   .
c   d i v i s i o n  b y  z e r o  will result if  t(left) = t(left+1)
c  jhigh,
c  index.....integers which determine the order  jout = max(jhigh,
c        (j+1)*(index-1))  of the b-splines whose values at  x  are to
c        be returned.  index  is used to avoid recalculations when seve-
c        ral columns of the triangular array of b-spline values are nee-
c        ded (e.g., in  bsplpp  or in  bsplvd ). precisely,
c                     if  index = 1 ,
c        the calculation starts from scratch and the entire triangular
c        array of b-spline values of orders 1,2,...,jhigh  is generated
c        order by order , i.e., column by column .
c                     if  index = 2 ,
c        only the b-spline values of order  j+1, j+2, ..., jout  are ge-
c        nerated, the assumption being that  biatx , j , deltal , deltar
c        are, on entry, as they were on exit at the previous call.
c           in particular, if  jhigh = 0, then  jout = j+1, i.e., just
c        the next column of b-spline values is generated.
c
c  w a r n i n g . . .  the restriction   jout .le. jmax (= 20)  is im-
c        posed arbitrarily by the dimension statement for  deltal  and
c        deltar  below, but is  n o w h e r e  c h e c k e d  for .
c
c  x.....the point at which the b-splines are to be evaluated.
c  left.....an integer chosen (usually) so that
c                  t(left) .le. x .le. t(left+1)  .
c
c******  o u t p u t  ******
c  biatx.....array of length  jout , with  biatx(i)  containing the val-
c        ue at  x  of the polynomial of order  jout  which agrees with
c        the b-spline  b(left-jout+i,jout,t)  on the interval (t(left),
c        t(left+1)) .
c
c******  m e t h o d  ******
c  the recurrence relation
c
c                       x - t(i)              t(i+j+1) - x
c     b(i,j+1)(x)  =  -----------b(i,j)(x) + ---------------b(i+1,j)(x)
c                     t(i+j)-t(i)            t(i+j+1)-t(i+1)
c
c  is used (repeatedly) to generate the (j+1)-vector  b(left-j,j+1)(x),
c  ...,b(left,j+1)(x)  from the j-vector  b(left-j+1,j)(x),...,
c  b(left,j)(x), storing the new values in  biatx  over the old. the
c  facts that
c            b(i,1) = 1  if  t(i) .le. x .lt. t(i+1)
c  and that
c            b(i,j)(x) = 0  unless  t(i) .le. x .lt. t(i+j)
c  are used. the particular organization of the calculations follows al-
c  gorithm  (8)  in chapter x of the text.
c
      integer index,jhigh,left,   i,j,jmax,jp1
      parameter (jmax = 20)
      real biatx(jhigh),t(1),x,   deltal(jmax),deltar(jmax),saved,term
C     real biatx(jhigh),t(1),x,   deltal(20),deltar(20),saved,term
c     dimension biatx(jout), t(left+jout)
c  current fortran standard makes it impossible to specify the length of
c  t  and of  biatx  precisely without the introduction of otherwise
c  superfluous additional arguments.
      data j/1/
      save j,deltal,deltar 
c
                                        go to (10,20), index
   10 j = 1
      biatx(1) = 1.
      if (j .ge. jhigh)                 go to 99
c
   20    jp1 = j + 1
         deltar(j) = t(left+j) - x
         deltal(j) = x - t(left+1-j)
         saved = 0.
         do 26 i=1,j
            term = biatx(i)/(deltar(i) + deltal(jp1-i))
            biatx(i) = saved + deltar(i)*term
   26       saved = deltal(jp1-i)*term
         biatx(jp1) = saved
         j = jp1
         if (j .lt. jhigh)              go to 20
c
   99                                   return
      end

      subroutine banfac ( w, nroww, nrow, nbandl, nbandu, iflag )
c  from  * a practical guide to splines *  by c. de boor    
c  returns in  w  the lu-factorization (without pivoting) of the banded
c  matrix  a  of order  nrow  with  (nbandl + 1 + nbandu) bands or diag-
c  onals in the work array  w .
c
c******  i n p u t  ******
c  w.....work array of size  (nroww,nrow)  containing the interesting
c        part of a banded matrix  a , with the diagonals or bands of  a
c        stored in the rows of  w , while columns of  a  correspond to
c        columns of  w . this is the storage mode used in  linpack  and
c        results in efficient innermost loops.
c           explicitly,  a  has  nbandl  bands below the diagonal
c                            +     1     (main) diagonal
c                            +   nbandu  bands above the diagonal
c        and thus, with    middle = nbandu + 1,
c          a(i+j,j)  is in  w(i+middle,j)  for i=-nbandu,...,nbandl
c                                              j=1,...,nrow .
c        for example, the interesting entries of a (1,2)-banded matrix
c        of order  9  would appear in the first  1+1+2 = 4  rows of  w
c        as follows.
c
c
c
c
c
c        all other entries of  w  not identified in this way with an en-
c        try of  a  are never referenced .
c  nroww.....row dimension of the work array  w .
c        must be  .ge.  nbandl + 1 + nbandu  .
c  nbandl.....number of bands of  a  below the main diagonal
c  nbandu.....number of bands of  a  above the main diagonal .
c
c******  o u t p u t  ******
c  iflag.....integer indicating success( = 1) or failure ( = 2) .
c     if  iflag = 1, then
c  w.....contains the lu-factorization of  a  into a unit lower triangu-
c        lar matrix  l  and an upper triangular matrix  u (both banded)
c        and stored in customary fashion over the corresponding entries
c        of  a . this makes it possible to solve any particular linear
c        system  a*x = b  for  x  by a
c              call banslv ( w, nroww, nrow, nbandl, nbandu, b )
c        with the solution x  contained in  b  on return .
c     if  iflag = 2, then
c        one of  nrow-1, nbandl,nbandu failed to be nonnegative, or else
c        one of the potential pivots was found to be zero indicating
c        that  a  does not have an lu-factorization. this implies that
c        a  is singular in case it is totally positive .
c
c******  m e t h o d  ******
c     gauss elimination  w i t h o u t  pivoting is used. the routine is
c  intended for use with matrices  a  which do not require row inter-
c  changes during factorization, especially for the  t o t a l l y
c  p o s i t i v e  matrices which occur in spline calculations.
c     the routine should not be used for an arbitrary banded matrix.
c
      integer iflag,nbandl,nbandu,nrow,nroww,   i,ipk,j,jmax,k,kmax
     *                                        ,middle,midmk,nrowm1
      real w(nroww,nrow),   factor,pivot
c
      iflag = 1
      middle = nbandu + 1
c                         w(middle,.) contains the main diagonal of  a .
      nrowm1 = nrow - 1
      if (nrowm1)                       999,900,1
    1 if (nbandl .gt. 0)                go to 10
c                a is upper triangular. check that diagonal is nonzero .
      do 5 i=1,nrowm1
         if (w(middle,i) .eq. 0.)       go to 999
    5    continue
                                        go to 900
   10 if (nbandu .gt. 0)                go to 20
c              a is lower triangular. check that diagonal is nonzero and
c                 divide each column by its diagonal .
      do 15 i=1,nrowm1
         pivot = w(middle,i)
         if(pivot .eq. 0.)              go to 999
         jmax = min0(nbandl, nrow - i)
         do 15 j=1,jmax
   15       w(middle+j,i) = w(middle+j,i)/pivot
                                        return
c
c        a  is not just a triangular matrix. construct lu factorization
   20 do 50 i=1,nrowm1
c                                  w(middle,i)  is pivot for i-th step .
         pivot = w(middle,i)
         if (pivot .eq. 0.)             go to 999
c                 jmax  is the number of (nonzero) entries in column  i
c                     below the diagonal .
         jmax = min0(nbandl,nrow - i)
c              divide each entry in column  i  below diagonal by pivot .
         do 32 j=1,jmax
   32       w(middle+j,i) = w(middle+j,i)/pivot
c                 kmax  is the number of (nonzero) entries in row  i  to
c                     the right of the diagonal .
         kmax = min0(nbandu,nrow - i)
c                  subtract  a(i,i+k)*(i-th column) from (i+k)-th column
c                  (below row  i ) .
         do 40 k=1,kmax
            ipk = i + k
            midmk = middle - k
            factor = w(midmk,ipk)
            do 40 j=1,jmax
   40          w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)*factor
   50    continue
c                                       check the last diagonal entry .
  900 if (w(middle,nrow) .ne. 0.)       return
  999 iflag = 2
                                        return
      end

      subroutine banslv ( w, nroww, nrow, nbandl, nbandu, b )
c  from  * a practical guide to splines *  by c. de boor    
c  companion routine to  banfac . it returns the solution  x  of the
c  linear system  a*x = b  in place of  b , given the lu-factorization
c  for  a  in the workarray  w .
c
c******  i n p u t  ******
c  w, nroww,nrow,nbandl,nbandu.....describe the lu-factorization of a
c        banded matrix  a  of order  nrow  as constructed in  banfac .
c        for details, see  banfac .
c  b.....right side of the system to be solved .
c
c******  o u t p u t  ******
c  b.....contains the solution  x , of order  nrow .
c
c******  m e t h o d  ******
c     (with  a = l*u, as stored in  w,) the unit lower triangular system
c  l(u*x) = b  is solved for  y = u*x, and  y  stored in  b . then the
c  upper triangular system  u*x = y  is solved for  x  . the calcul-
c  ations are so arranged that the innermost loops stay within columns.
c
      integer nbandl,nbandu,nrow,nroww,   i,j,jmax,middle,nrowm1
      real w(nroww,nrow),b(nrow)
      middle = nbandu + 1
      if (nrow .eq. 1)                  go to 49
      nrowm1 = nrow - 1
      if (nbandl .eq. 0)                go to 30
c                                 forward pass
c            for i=1,2,...,nrow-1, subtract  right side(i)*(i-th column
c            of  l )  from right side  (below i-th row) .
      do 21 i=1,nrowm1
         jmax = min0(nbandl, nrow-i)
         do 21 j=1,jmax
   21       b(i+j) = b(i+j) - b(i)*w(middle+j,i)
c                                 backward pass
c            for i=nrow,nrow-1,...,1, divide right side(i) by i-th diag-
c            onal entry of  u, then subtract  right side(i)*(i-th column
c            of  u)  from right side  (above i-th row).
   30 if (nbandu .gt. 0)                go to 40
c                                a  is lower triangular .
      do 31 i=1,nrow
   31    b(i) = b(i)/w(1,i)
                                        return
   40 i = nrow    
   41    b(i) = b(i)/w(middle,i)
         jmax = min0(nbandu,i-1)
         do 45 j=1,jmax
   45       b(i-j) = b(i-j) - b(i)*w(middle-j,i)
         i = i - 1
         if (i .gt. 1)                  go to 41
   49 b(1) = b(1)/w(middle,1)
                                        return
      end

      real function bvalue ( t, bcoef, n, k, x, jderiv )
c  from  * a practical guide to splines *  by c. de boor    
c  calls  interv
c
c  calculates value at  x  of  jderiv-th derivative of spline from b-repr.
c  the spline is taken to be continuous from the right, EXCEPT at the
c  rightmost knot, where it is taken to be continuous from the left.
c
c******  i n p u t ******
c  t, bcoef, n, k......forms the b-representation of the spline  f  to
c        be evaluated. specifically,
c  t.....knot sequence, of length  n+k, assumed nondecreasing.
c  bcoef.....b-coefficient sequence, of length  n .
c  n.....length of  bcoef  and dimension of spline(k,t),
c        a s s u m e d  positive .
c  k.....order of the spline .
c
c  w a r n i n g . . .   the restriction  k .le. kmax (=20)  is imposed
c        arbitrarily by the dimension statement for  aj, dl, dr  below,
c        but is  n o w h e r e  c h e c k e d  for.
c
c  x.....the point at which to evaluate .
c  jderiv.....integer giving the order of the derivative to be evaluated
c        a s s u m e d  to be zero or positive.
c
c******  o u t p u t  ******
c  bvalue.....the value of the (jderiv)-th derivative of  f  at  x .
c
c******  m e t h o d  ******
c     The nontrivial knot interval  (t(i),t(i+1))  containing  x  is lo-
c  cated with the aid of  interv . The  k  b-coeffs of  f  relevant for
c  this interval are then obtained from  bcoef (or taken to be zero if
c  not explicitly available) and are then differenced  jderiv  times to
c  obtain the b-coeffs of  (d**jderiv)f  relevant for that interval.
c  Precisely, with  j = jderiv, we have from x.(12) of the text that
c
c     (d**j)f  =  sum ( bcoef(.,j)*b(.,k-j,t) )
c
c  where
c                   / bcoef(.),                     ,  j .eq. 0
c                   /
c    bcoef(.,j)  =  / bcoef(.,j-1) - bcoef(.-1,j-1)
c                   / ----------------------------- ,  j .gt. 0
c                   /    (t(.+k-j) - t(.))/(k-j)
c
c     Then, we use repeatedly the fact that
c
c    sum ( a(.)*b(.,m,t)(x) )  =  sum ( a(.,x)*b(.,m-1,t)(x) )
c  with
c                 (x - t(.))*a(.) + (t(.+m-1) - x)*a(.-1)
c    a(.,x)  =    ---------------------------------------
c                 (x - t(.))      + (t(.+m-1) - x)
c
c  to write  (d**j)f(x)  eventually as a linear combination of b-splines
c  of order  1 , and the coefficient for  b(i,1,t)(x)  must then be the
c  desired number  (d**j)f(x). (see x.(17)-(19) of text).
c
      integer jderiv,k,n,   i,ilo,imk,j,jc,jcmin,jcmax,jj,kmax,kmj,km1
     *                     ,mflag,nmi,jdrvp1
      parameter (kmax = 20)
C     real bcoef(n),t(1),x,   aj(20),dl(20),dr(20),fkmj
      real bcoef(n),t(1),x,   aj(kmax),dl(kmax),dr(kmax),fkmj
c     dimension t(n+k)
c  former fortran standard made it impossible to specify the length of  t
c  precisely without the introduction of otherwise superfluous addition-
c  al arguments.
      bvalue = 0.
      if (jderiv .ge. k)                go to 99
c
c  *** Find  i   s.t.   1 .le. i .lt. n+k   and   t(i) .lt. t(i+1)   and
c      t(i) .le. x .lt. t(i+1) . If no such i can be found,  x  lies
c      outside the support of  the spline  f , hence  bvalue = 0.
c      (The asymmetry in this choice of  i  makes  f  rightcontinuous, except
c      at  t(n+k) where it is leftcontinuous.)
      call interv ( t, n+k, x, i, mflag )
      if (mflag .ne. 0)                 go to 99
c  *** if k = 1 (and jderiv = 0), bvalue = bcoef(i).
      km1 = k - 1
      if (km1 .gt. 0)                   go to 1
      bvalue = bcoef(i)
                                        go to 99
c
c  *** store the k b-spline coefficients relevant for the knot interval
c     (t(i),t(i+1)) in aj(1),...,aj(k) and compute dl(j) = x - t(i+1-j),
c     dr(j) = t(i+j) - x, j=1,...,k-1 . set any of the aj not obtainable
c     from input to zero. set any t.s not obtainable equal to t(1) or
c     to t(n+k) appropriately.
    1 jcmin = 1
      imk = i - k
      if (imk .ge. 0)                   go to 8
      jcmin = 1 - imk
      do 5 j=1,i
    5    dl(j) = x - t(i+1-j)
      do 6 j=i,km1
         aj(k-j) = 0.
    6    dl(j) = dl(i)
                                        go to 10
    8 do 9 j=1,km1
    9    dl(j) = x - t(i+1-j)
c
   10 jcmax = k
      nmi = n - i
      if (nmi .ge. 0)                   go to 18
      jcmax = k + nmi
      do 15 j=1,jcmax
   15    dr(j) = t(i+j) - x
      do 16 j=jcmax,km1
         aj(j+1) = 0.
   16    dr(j) = dr(jcmax)
                                        go to 20
   18 do 19 j=1,km1
   19    dr(j) = t(i+j) - x
c
   20 do 21 jc=jcmin,jcmax
   21    aj(jc) = bcoef(imk + jc)
c
c               *** difference the coefficients  jderiv  times.
      if (jderiv .eq. 0)                go to 30
      do 23 j=1,jderiv
         kmj = k-j
         fkmj = float(kmj)
         ilo = kmj
         do 23 jj=1,kmj
            aj(jj) = ((aj(jj+1) - aj(jj))/(dl(ilo) + dr(jj)))*fkmj
   23       ilo = ilo - 1
c
c  *** compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative,
c     given its relevant b-spline coeffs in aj(1),...,aj(k-jderiv).
   30 if (jderiv .eq. km1)              go to 39
      jdrvp1 = jderiv + 1     
      do 33 j=jdrvp1,km1
         kmj = k-j
         ilo = kmj
         do 33 jj=1,kmj
            aj(jj) = (aj(jj+1)*dl(ilo) + aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
   33       ilo = ilo - 1
   39 bvalue = aj(1)
c
   99                                   return
      end

      subroutine interv ( xt, lxt, x, left, mflag )
c  from  * a practical guide to splines *  by C. de Boor    
c  computes  left = max( i :  xt(i) .lt. xt(lxt) .and.  xt(i) .le. x )  .
c
c******  i n p u t  ******
c  xt.....a real sequence, of length  lxt , assumed to be nondecreasing
c  lxt.....number of terms in the sequence  xt .
c  x.....the point whose location with respect to the sequence  xt  is
c        to be determined.
c
c******  o u t p u t  ******
c  left, mflag.....both integers, whose value is
c
c   1     -1      if               x .lt.  xt(1)
c   i      0      if   xt(i)  .le. x .lt. xt(i+1)
c   i      0      if   xt(i)  .lt. x .eq. xt(i+1) .eq. xt(lxt)
c   i      1      if   xt(i)  .lt.        xt(i+1) .eq. xt(lxt) .lt. x
c
c        In particular,  mflag = 0  is the 'usual' case.  mflag .ne. 0
c        indicates that  x  lies outside the CLOSED interval
c        xt(1) .le. y .le. xt(lxt) . The asymmetric treatment of the
c        intervals is due to the decision to make all pp functions cont-
c        inuous from the right, but, by returning  mflag = 0  even if
C        x = xt(lxt), there is the option of having the computed pp function
c        continuous from the left at  xt(lxt) .
c
c******  m e t h o d  ******
c  The program is designed to be efficient in the common situation that
c  it is called repeatedly, with  x  taken from an increasing or decrea-
c  sing sequence. This will happen, e.g., when a pp function is to be
c  graphed. The first guess for  left  is therefore taken to be the val-
c  ue returned at the previous call and stored in the  l o c a l  varia-
c  ble  ilo . A first check ascertains that  ilo .lt. lxt (this is nec-
c  essary since the present call may have nothing to do with the previ-
c  ous call). Then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
c  ilo  and are done after just three comparisons.
c     Otherwise, we repeatedly double the difference  istep = ihi - ilo
c  while also moving  ilo  and  ihi  in the direction of  x , until
c                      xt(ilo) .le. x .lt. xt(ihi) ,
c  after which we use bisection to get, in addition, ilo+1 = ihi .
c  left = ilo  is then returned.
c
      integer left,lxt,mflag,   ihi,ilo,istep,middle
      real x,xt(lxt)
      data ilo /1/
      save ilo  
      ihi = ilo + 1
      if (ihi .lt. lxt)                 go to 20
         if (x .ge. xt(lxt))            go to 110
         if (lxt .le. 1)                go to 90
         ilo = lxt - 1
         ihi = lxt
c
   20 if (x .ge. xt(ihi))               go to 40
      if (x .ge. xt(ilo))               go to 100
c
c              **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
      istep = 1
   31    ihi = ilo
         ilo = ihi - istep
         if (ilo .le. 1)                go to 35
         if (x .ge. xt(ilo))            go to 50
         istep = istep*2
                                        go to 31
   35 ilo = 1
      if (x .lt. xt(1))                 go to 90
                                        go to 50
c              **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
   40 istep = 1
   41    ilo = ihi
         ihi = ilo + istep
         if (ihi .ge. lxt)              go to 45
         if (x .lt. xt(ihi))            go to 50
         istep = istep*2
                                        go to 41
   45 if (x .ge. xt(lxt))               go to 110
      ihi = lxt
c
c           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
   50 middle = (ilo + ihi)/2
      if (middle .eq. ilo)              go to 100
c     note. it is assumed that middle = ilo in case ihi = ilo+1 .
      if (x .lt. xt(middle))            go to 53
         ilo = middle
                                        go to 50
   53    ihi = middle
                                        go to 50
c**** set output and return.
   90 mflag = -1
      left = 1
                                        return
  100 mflag = 0
      left = ilo
                                        return
  110 mflag = 1
          if (x .eq. xt(lxt)) mflag = 0
      left = lxt
  111 if (left .eq. 1)                  return
          left = left - 1
          if (xt(left) .lt. xt(lxt))    return
                                        go to 111
      end

*DECK BINTK
      SUBROUTINE BINTK (X, Y, T, N, K, BCOEF, Q, WORK)
C***BEGIN PROLOGUE  BINTK
C***PURPOSE  Compute the B-representation of a spline which interpolates
C            given data.
C***LIBRARY   SLATEC
C***CATEGORY  E1A
C***TYPE      SINGLE PRECISION (BINTK-S, DBINTK-D)
C***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Abstract
C
C         BINTK is the SPLINT routine of the reference.
C
C         BINTK produces the B-spline coefficients, BCOEF, of the
C         B-spline of order K with knots T(I), I=1,...,N+K, which
C         takes on the value Y(I) at X(I), I=1,...,N.  The spline or
C         any of its derivatives can be evaluated by calls to BVALU.
C         The I-th equation of the linear system A*BCOEF = B for the
C         coefficients of the interpolant enforces interpolation at
C         X(I)), I=1,...,N.  Hence, B(I) = Y(I), all I, and A is
C         a band matrix with 2K-1 bands if A is invertible. The matrix
C         A is generated row by row and stored, diagonal by diagonal,
C         in the rows of Q, with the main diagonal going into row K.
C         The banded system is then solved by a call to BNFAC (which
C         constructs the triangular factorization for A and stores it
C         again in Q), followed by a call to BNSLV (which then
C         obtains the solution BCOEF by substitution). BNFAC does no
C         pivoting, since the total positivity of the matrix A makes
C         this unnecessary.  The linear system to be solved is
C         (theoretically) invertible if and only if
C                 T(I) .LT. X(I)) .LT. T(I+K),        all I.
C         Equality is permitted on the left for I=1 and on the right
C         for I=N when K knots are used at X(1) or X(N).  Otherwise,
C         violation of this condition is certain to lead to an error.
C
C     Description of Arguments
C         Input
C           X       - vector of length N containing data point abscissa
C                     in strictly increasing order.
C           Y       - corresponding vector of length N containing data
C                     point ordinates.
C           T       - knot vector of length N+K
C                     since T(1),..,T(K) .LE. X(1) and T(N+1),..,T(N+K)
C                     .GE. X(N), this leaves only N-K knots (not nec-
C                     essarily X(I)) values) interior to (X(1),X(N))
C           N       - number of data points, N .GE. K
C           K       - order of the spline, K .GE. 1
C
C         Output
C           BCOEF   - a vector of length N containing the B-spline
C                     coefficients
C           Q       - a work vector of length (2*K-1)*N, containing
C                     the triangular factorization of the coefficient
C                     matrix of the linear system being solved.  The
C                     coefficients for the interpolant of an
C                     additional data set (X(I)),YY(I)), I=1,...,N
C                     with the same abscissa can be obtained by loading
C                     YY into BCOEF and then executing
C                         CALL BNSLV (Q,2K-1,N,K-1,K-1,BCOEF)
C           WORK    - work vector of length 2*K
C
C     Error Conditions
C         Improper  input is a fatal error
C         Singular system of equations is a fatal error
C
C***REFERENCES  D. E. Amos, Computation with splines and B-splines,
C                 Report SAND78-1968, Sandia Laboratories, March 1979.
C               Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C               Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C***ROUTINES CALLED  BNFAC, BNSLV, BSPVN, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  BINTK
C
      INTEGER IFLAG, IWORK, K, N, I, ILP1MX, J, JJ, KM1, KPKM2, LEFT,
     1 LENQ, NP1
      REAL BCOEF(*), Y(*), Q(*), T(*), X(*), XI, WORK(*)
C     DIMENSION Q(2*K-1,N), T(N+K)
C***FIRST EXECUTABLE STATEMENT  BINTK
      IF(K.LT.1) GO TO 100
      IF(N.LT.K) GO TO 105
      JJ = N - 1
      IF(JJ.EQ.0) GO TO 6
      DO 5 I=1,JJ
      IF(X(I).GE.X(I+1)) GO TO 110
    5 CONTINUE
    6 CONTINUE
      NP1 = N + 1
      KM1 = K - 1
      KPKM2 = 2*KM1
      LEFT = K
C                ZERO OUT ALL ENTRIES OF Q
      LENQ = N*(K+KM1)
      DO 10 I=1,LENQ
        Q(I) = 0.0E0
   10 CONTINUE
C
C  ***   LOOP OVER I TO CONSTRUCT THE  N  INTERPOLATION EQUATIONS
      DO 50 I=1,N
        XI = X(I)
        ILP1MX = MIN(I+K,NP1)
C        *** FIND  LEFT  IN THE CLOSED INTERVAL (I,I+K-1) SUCH THAT
C                T(LEFT) .LE. X(I) .LT. T(LEFT+1)
C        MATRIX IS SINGULAR IF THIS IS NOT POSSIBLE
        LEFT = MAX(LEFT,I)
        IF (XI.LT.T(LEFT)) GO TO 80
   20   IF (XI.LT.T(LEFT+1)) GO TO 30
        LEFT = LEFT + 1
        IF (LEFT.LT.ILP1MX) GO TO 20
        LEFT = LEFT - 1
        IF (XI.GT.T(LEFT+1)) GO TO 80
C        *** THE I-TH EQUATION ENFORCES INTERPOLATION AT XI, HENCE
C        A(I,J) = B(J,K,T)(XI), ALL J. ONLY THE  K  ENTRIES WITH  J =
C        LEFT-K+1,...,LEFT ACTUALLY MIGHT BE NONZERO. THESE  K  NUMBERS
C        ARE RETURNED, IN  BCOEF (USED FOR TEMP. STORAGE HERE), BY THE
C        FOLLOWING
   30   CALL BSPVN(T, K, K, 1, XI, LEFT, BCOEF, WORK, IWORK)
C        WE THEREFORE WANT  BCOEF(J) = B(LEFT-K+J)(XI) TO GO INTO
C        A(I,LEFT-K+J), I.E., INTO  Q(I-(LEFT+J)+2*K,(LEFT+J)-K) SINCE
C        A(I+J,J)  IS TO GO INTO  Q(I+K,J), ALL I,J,  IF WE CONSIDER  Q
C        AS A TWO-DIM. ARRAY , WITH  2*K-1  ROWS (SEE COMMENTS IN
C        BNFAC). IN THE PRESENT PROGRAM, WE TREAT  Q  AS AN EQUIVALENT
C        ONE-DIMENSIONAL ARRAY (BECAUSE OF FORTRAN RESTRICTIONS ON
C        DIMENSION STATEMENTS) . WE THEREFORE WANT  BCOEF(J) TO GO INTO
C        ENTRY
C            I -(LEFT+J) + 2*K + ((LEFT+J) - K-1)*(2*K-1)
C                   =  I-LEFT+1 + (LEFT -K)*(2*K-1) + (2*K-2)*J
C        OF  Q .
        JJ = I - LEFT + 1 + (LEFT-K)*(K+KM1)
        DO 40 J=1,K
          JJ = JJ + KPKM2
          Q(JJ) = BCOEF(J)
   40   CONTINUE
   50 CONTINUE
C
C     ***OBTAIN FACTORIZATION OF  A  , STORED AGAIN IN  Q.
      CALL BNFAC(Q, K+KM1, N, KM1, KM1, IFLAG)
      GO TO (60, 90), IFLAG
C     *** SOLVE  A*BCOEF = Y  BY BACKSUBSTITUTION
   60 DO 70 I=1,N
        BCOEF(I) = Y(I)
   70 CONTINUE
      CALL BNSLV(Q, K+KM1, N, KM1, KM1, BCOEF)
      RETURN
C
C
   80 CONTINUE
      CALL XERMSG ('SLATEC', 'BINTK',
     +   'SOME ABSCISSA WAS NOT IN THE SUPPORT OF THE CORRESPONDING ' //
     +   'BASIS FUNCTION AND THE SYSTEM IS SINGULAR.', 2, 1)
      RETURN
   90 CONTINUE
      CALL XERMSG ('SLATEC', 'BINTK',
     +   'THE SYSTEM OF SOLVER DETECTS A SINGULAR SYSTEM ALTHOUGH ' //
     +   'THE THEORETICAL CONDITIONS FOR A SOLUTION WERE SATISFIED.',
     +   8, 1)
      RETURN
  100 CONTINUE
      CALL XERMSG ('SLATEC', 'BINTK', 'K DOES NOT SATISFY K.GE.1', 2,
     +   1)
      RETURN
  105 CONTINUE
      CALL XERMSG ('SLATEC', 'BINTK', 'N DOES NOT SATISFY N.GE.K', 2,
     +   1)
      RETURN
  110 CONTINUE
      CALL XERMSG ('SLATEC', 'BINTK',
     +   'X(I) DOES NOT SATISFY X(I).LT.X(I+1) FOR SOME I', 2, 1)
      RETURN
      END
*DECK BNFAC
      SUBROUTINE BNFAC (W, NROWW, NROW, NBANDL, NBANDU, IFLAG)
C***BEGIN PROLOGUE  BNFAC
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BINT4 and BINTK
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (BNFAC-S, DBNFAC-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C  BNFAC is the BANFAC routine from
C        * A Practical Guide to Splines *  by C. de Boor
C
C  Returns in  W  the lu-factorization (without pivoting) of the banded
C  matrix  A  of order  NROW  with  (NBANDL + 1 + NBANDU) bands or diag-
C  onals in the work array  W .
C
C *****  I N P U T  ******
C  W.....Work array of size  (NROWW,NROW)  containing the interesting
C        part of a banded matrix  A , with the diagonals or bands of  A
C        stored in the rows of  W , while columns of  A  correspond to
C        columns of  W . This is the storage mode used in  LINPACK  and
C        results in efficient innermost loops.
C           Explicitly,  A  has  NBANDL  bands below the diagonal
C                            +     1     (main) diagonal
C                            +   NBANDU  bands above the diagonal
C        and thus, with    MIDDLE = NBANDU + 1,
C          A(I+J,J)  is in  W(I+MIDDLE,J)  for I=-NBANDU,...,NBANDL
C                                              J=1,...,NROW .
C        For example, the interesting entries of A (1,2)-banded matrix
C        of order  9  would appear in the first  1+1+2 = 4  rows of  W
C        as follows.
C                          13 24 35 46 57 68 79
C                       12 23 34 45 56 67 78 89
C                    11 22 33 44 55 66 77 88 99
C                    21 32 43 54 65 76 87 98
C
C        All other entries of  W  not identified in this way with an en-
C        try of  A  are never referenced .
C  NROWW.....Row dimension of the work array  W .
C        must be  .GE.  NBANDL + 1 + NBANDU  .
C  NBANDL.....Number of bands of  A  below the main diagonal
C  NBANDU.....Number of bands of  A  above the main diagonal .
C
C *****  O U T P U T  ******
C  IFLAG.....Integer indicating success( = 1) or failure ( = 2) .
C     If  IFLAG = 1, then
C  W.....contains the LU-factorization of  A  into a unit lower triangu-
C        lar matrix  L  and an upper triangular matrix  U (both banded)
C        and stored in customary fashion over the corresponding entries
C        of  A . This makes it possible to solve any particular linear
C        system  A*X = B  for  X  by A
C              CALL BNSLV ( W, NROWW, NROW, NBANDL, NBANDU, B )
C        with the solution X  contained in  B  on return .
C     If  IFLAG = 2, then
C        one of  NROW-1, NBANDL,NBANDU failed to be nonnegative, or else
C        one of the potential pivots was found to be zero indicating
C        that  A  does not have an LU-factorization. This implies that
C        A  is singular in case it is totally positive .
C
C *****  M E T H O D  ******
C     Gauss elimination  W I T H O U T  pivoting is used. The routine is
C  intended for use with matrices  A  which do not require row inter-
C  changes during factorization, especially for the  T O T A L L Y
C  P O S I T I V E  matrices which occur in spline calculations.
C     The routine should not be used for an arbitrary banded matrix.
C
C***SEE ALSO  BINT4, BINTK
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  BNFAC
C
      INTEGER IFLAG, NBANDL, NBANDU, NROW, NROWW, I, IPK, J, JMAX, K,
     1 KMAX, MIDDLE, MIDMK, NROWM1
      REAL W(NROWW,*), FACTOR, PIVOT
C
C***FIRST EXECUTABLE STATEMENT  BNFAC
      IFLAG = 1
      MIDDLE = NBANDU + 1
C                         W(MIDDLE,.) CONTAINS THE MAIN DIAGONAL OF  A .
      NROWM1 = NROW - 1
      IF (NROWM1) 120, 110, 10
   10 IF (NBANDL.GT.0) GO TO 30
C                A IS UPPER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO .
      DO 20 I=1,NROWM1
        IF (W(MIDDLE,I).EQ.0.0E0) GO TO 120
   20 CONTINUE
      GO TO 110
   30 IF (NBANDU.GT.0) GO TO 60
C              A IS LOWER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO AND
C                 DIVIDE EACH COLUMN BY ITS DIAGONAL .
      DO 50 I=1,NROWM1
        PIVOT = W(MIDDLE,I)
        IF (PIVOT.EQ.0.0E0) GO TO 120
        JMAX = MIN(NBANDL,NROW-I)
        DO 40 J=1,JMAX
          W(MIDDLE+J,I) = W(MIDDLE+J,I)/PIVOT
   40   CONTINUE
   50 CONTINUE
      RETURN
C
C        A  IS NOT JUST A TRIANGULAR MATRIX. CONSTRUCT LU FACTORIZATION
   60 DO 100 I=1,NROWM1
C                                  W(MIDDLE,I)  IS PIVOT FOR I-TH STEP .
        PIVOT = W(MIDDLE,I)
        IF (PIVOT.EQ.0.0E0) GO TO 120
C                 JMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN COLUMN  I
C                     BELOW THE DIAGONAL .
        JMAX = MIN(NBANDL,NROW-I)
C              DIVIDE EACH ENTRY IN COLUMN  I  BELOW DIAGONAL BY PIVOT .
        DO 70 J=1,JMAX
          W(MIDDLE+J,I) = W(MIDDLE+J,I)/PIVOT
   70   CONTINUE
C                 KMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN ROW  I  TO
C                     THE RIGHT OF THE DIAGONAL .
        KMAX = MIN(NBANDU,NROW-I)
C                  SUBTRACT  A(I,I+K)*(I-TH COLUMN) FROM (I+K)-TH COLUMN
C                  (BELOW ROW  I ) .
        DO 90 K=1,KMAX
          IPK = I + K
          MIDMK = MIDDLE - K
          FACTOR = W(MIDMK,IPK)
          DO 80 J=1,JMAX
            W(MIDMK+J,IPK) = W(MIDMK+J,IPK) - W(MIDDLE+J,I)*FACTOR
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE
C                                       CHECK THE LAST DIAGONAL ENTRY .
  110 IF (W(MIDDLE,NROW).NE.0.0E0) RETURN
  120 IFLAG = 2
      RETURN
      END
*DECK BNSLV
      SUBROUTINE BNSLV (W, NROWW, NROW, NBANDL, NBANDU, B)
C***BEGIN PROLOGUE  BNSLV
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BINT4 and BINTK
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (BNSLV-S, DBNSLV-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C  BNSLV is the BANSLV routine from
C        * A Practical Guide to Splines *  by C. de Boor
C
C  Companion routine to  BNFAC . It returns the solution  X  of the
C  linear system  A*X = B  in place of  B , given the LU-factorization
C  for  A  in the work array  W from BNFAC.
C
C *****  I N P U T  ******
C  W, NROWW,NROW,NBANDL,NBANDU.....Describe the LU-factorization of a
C        banded matrix  A  of order  NROW  as constructed in  BNFAC .
C        For details, see  BNFAC .
C  B.....Right side of the system to be solved .
C
C *****  O U T P U T  ******
C  B.....Contains the solution  X , of order  NROW .
C
C *****  M E T H O D  ******
C     (With  A = L*U, as stored in  W,) the unit lower triangular system
C  L(U*X) = B  is solved for  Y = U*X, and  Y  stored in  B . Then the
C  upper triangular system  U*X = Y  is solved for  X  . The calcul-
C  ations are so arranged that the innermost loops stay within columns.
C
C***SEE ALSO  BINT4, BINTK
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  BNSLV
C
      INTEGER NBANDL, NBANDU, NROW, NROWW, I, J, JMAX, MIDDLE, NROWM1
      REAL W(NROWW,*), B(*)
C***FIRST EXECUTABLE STATEMENT  BNSLV
      MIDDLE = NBANDU + 1
      IF (NROW.EQ.1) GO TO 80
      NROWM1 = NROW - 1
      IF (NBANDL.EQ.0) GO TO 30
C                                 FORWARD PASS
C            FOR I=1,2,...,NROW-1, SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN
C            OF  L )  FROM RIGHT SIDE  (BELOW I-TH ROW) .
      DO 20 I=1,NROWM1
        JMAX = MIN(NBANDL,NROW-I)
        DO 10 J=1,JMAX
          B(I+J) = B(I+J) - B(I)*W(MIDDLE+J,I)
   10   CONTINUE
   20 CONTINUE
C                                 BACKWARD PASS
C            FOR I=NROW,NROW-1,...,1, DIVIDE RIGHT SIDE(I) BY I-TH DIAG-
C            ONAL ENTRY OF  U, THEN SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN
C            OF  U)  FROM RIGHT SIDE  (ABOVE I-TH ROW).
   30 IF (NBANDU.GT.0) GO TO 50
C                                A  IS LOWER TRIANGULAR .
      DO 40 I=1,NROW
        B(I) = B(I)/W(1,I)
   40 CONTINUE
      RETURN
   50 I = NROW
   60 B(I) = B(I)/W(MIDDLE,I)
      JMAX = MIN(NBANDU,I-1)
      DO 70 J=1,JMAX
        B(I-J) = B(I-J) - B(I)*W(MIDDLE-J,I)
   70 CONTINUE
      I = I - 1
      IF (I.GT.1) GO TO 60
   80 B(1) = B(1)/W(MIDDLE,1)
      RETURN
      END
*DECK BSPVN
      SUBROUTINE BSPVN (T, JHIGH, K, INDEX, X, ILEFT, VNIKX, WORK,
     +   IWORK)
C***BEGIN PROLOGUE  BSPVN
C***PURPOSE  Calculate the value of all (possibly) nonzero basis
C            functions at X.
C***LIBRARY   SLATEC
C***CATEGORY  E3, K6
C***TYPE      SINGLE PRECISION (BSPVN-S, DBSPVN-D)
C***KEYWORDS  EVALUATION OF B-SPLINE
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Abstract
C         BSPVN is the BSPLVN routine of the reference.
C
C         BSPVN calculates the value of all (possibly) nonzero basis
C         functions at X of order MAX(JHIGH,(J+1)*(INDEX-1)), where
C         T(K) .LE. X .LE. T(N+1) and J=IWORK is set inside the routine
C         on the first call when INDEX=1.  ILEFT is such that T(ILEFT)
C         .LE. X .LT. T(ILEFT+1).  A call to INTRV(T,N+1,X,ILO,ILEFT,
C         MFLAG) produces the proper ILEFT.  BSPVN calculates using the
C         basic algorithm needed in BSPVD.  If only basis functions are
C         desired, setting JHIGH=K and INDEX=1 can be faster than
C         calling BSPVD, but extra coding is required for derivatives
C         (INDEX=2) and BSPVD is set up for this purpose.
C
C         Left limiting values are set up as described in BSPVD.
C
C     Description of Arguments
C         Input
C          T       - knot vector of length N+K, where
C                    N = number of B-spline basis functions
C                    N = sum of knot multiplicities-K
C          JHIGH   - order of B-spline, 1 .LE. JHIGH .LE. K
C          K       - highest possible order
C          INDEX   - INDEX = 1 gives basis functions of order JHIGH
C                          = 2 denotes previous entry with WORK, IWORK
C                              values saved for subsequent calls to
C                              BSPVN.
C          X       - argument of basis functions,
C                    T(K) .LE. X .LE. T(N+1)
C          ILEFT   - largest integer such that
C                    T(ILEFT) .LE. X .LT. T(ILEFT+1)
C
C         Output
C          VNIKX   - vector of length K for spline values.
C          WORK    - a work vector of length 2*K
C          IWORK   - a work parameter.  Both WORK and IWORK contain
C                    information necessary to continue for INDEX = 2.
C                    When INDEX = 1 exclusively, these are scratch
C                    variables and can be used for other purposes.
C
C     Error Conditions
C         Improper input is a fatal error.
C
C***REFERENCES  Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  BSPVN
C
      INTEGER ILEFT, IMJP1, INDEX, IPJ, IWORK, JHIGH, JP1, JP1ML, K, L
      REAL T, VM, VMPREV, VNIKX, WORK, X
C     DIMENSION T(ILEFT+JHIGH)
      DIMENSION T(*), VNIKX(*), WORK(*)
C     CONTENT OF J, DELTAM, DELTAP IS EXPECTED UNCHANGED BETWEEN CALLS.
C     WORK(I) = DELTAP(I), WORK(K+I) = DELTAM(I), I = 1,K
C***FIRST EXECUTABLE STATEMENT  BSPVN
      IF(K.LT.1) GO TO 90
      IF(JHIGH.GT.K .OR. JHIGH.LT.1) GO TO 100
      IF(INDEX.LT.1 .OR. INDEX.GT.2) GO TO 105
      IF(X.LT.T(ILEFT) .OR. X.GT.T(ILEFT+1)) GO TO 110
      GO TO (10, 20), INDEX
   10 IWORK = 1
      VNIKX(1) = 1.0E0
      IF (IWORK.GE.JHIGH) GO TO 40
C
   20 IPJ = ILEFT + IWORK
      WORK(IWORK) = T(IPJ) - X
      IMJP1 = ILEFT - IWORK + 1
      WORK(K+IWORK) = X - T(IMJP1)
      VMPREV = 0.0E0
      JP1 = IWORK + 1
      DO 30 L=1,IWORK
        JP1ML = JP1 - L
        VM = VNIKX(L)/(WORK(L)+WORK(K+JP1ML))
        VNIKX(L) = VM*WORK(L) + VMPREV
        VMPREV = VM*WORK(K+JP1ML)
   30 CONTINUE
      VNIKX(JP1) = VMPREV
      IWORK = JP1
      IF (IWORK.LT.JHIGH) GO TO 20
C
   40 RETURN
C
C
   90 CONTINUE
      CALL XERMSG ('SLATEC', 'BSPVN', 'K DOES NOT SATISFY K.GE.1', 2,
     +   1)
      RETURN
  100 CONTINUE
      CALL XERMSG ('SLATEC', 'BSPVN',
     +   'JHIGH DOES NOT SATISFY 1.LE.JHIGH.LE.K', 2, 1)
      RETURN
  105 CONTINUE
      CALL XERMSG ('SLATEC', 'BSPVN', 'INDEX IS NOT 1 OR 2', 2, 1)
      RETURN
  110 CONTINUE
      CALL XERMSG ('SLATEC', 'BSPVN',
     +   'X DOES NOT SATISFY T(ILEFT).LE.X.LE.T(ILEFT+1)', 2, 1)
      RETURN
      END
*DECK FDUMP
      SUBROUTINE FDUMP
C***BEGIN PROLOGUE  FDUMP
C***PURPOSE  Symbolic dump (should be locally written).
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3
C***TYPE      ALL (FDUMP-A)
C***KEYWORDS  ERROR, XERMSG
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C        ***Note*** Machine Dependent Routine
C        FDUMP is intended to be replaced by a locally written
C        version which produces a symbolic dump.  Failing this,
C        it should be replaced by a version which prints the
C        subprogram nesting list.  Note that this dump must be
C        printed on each of up to five files, as indicated by the
C        XGETUA routine.  See XSETUA and XGETUA for details.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
*DECK I1MACH
      INTEGER FUNCTION I1MACH (I)
C***BEGIN PROLOGUE  I1MACH
C***PURPOSE  Return integer machine dependent constants.
C***LIBRARY   SLATEC
C***CATEGORY  R1
C***TYPE      INTEGER (I1MACH-I)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Fox, P. A., (Bell Labs)
C           Hall, A. D., (Bell Labs)
C           Schryer, N. L., (Bell Labs)
C***DESCRIPTION
C
C   I1MACH can be used to obtain machine-dependent parameters for the
C   local machine environment.  It is a function subprogram with one
C   (input) argument and can be referenced as follows:
C
C        K = I1MACH(I)
C
C   where I=1,...,16.  The (output) value of K above is determined by
C   the (input) value of I.  The results for various values of I are
C   discussed below.
C
C   I/O unit numbers:
C     I1MACH( 1) = the standard input unit.
C     I1MACH( 2) = the standard output unit.
C     I1MACH( 3) = the standard punch unit.
C     I1MACH( 4) = the standard error message unit.
C
C   Words:
C     I1MACH( 5) = the number of bits per integer storage unit.
C     I1MACH( 6) = the number of characters per integer storage unit.
C
C   Integers:
C     assume integers are represented in the S-digit, base-A form
C
C                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C                where 0 .LE. X(I) .LT. A for I=0,...,S-1.
C     I1MACH( 7) = A, the base.
C     I1MACH( 8) = S, the number of base-A digits.
C     I1MACH( 9) = A**S - 1, the largest magnitude.
C
C   Floating-Point Numbers:
C     Assume floating-point numbers are represented in the T-digit,
C     base-B form
C                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C                where 0 .LE. X(I) .LT. B for I=1,...,T,
C                0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
C     I1MACH(10) = B, the base.
C
C   Single-Precision:
C     I1MACH(11) = T, the number of base-B digits.
C     I1MACH(12) = EMIN, the smallest exponent E.
C     I1MACH(13) = EMAX, the largest exponent E.
C
C   Double-Precision:
C     I1MACH(14) = T, the number of base-B digits.
C     I1MACH(15) = EMIN, the smallest exponent E.
C     I1MACH(16) = EMAX, the largest exponent E.
C
C   To alter this function for a particular environment, the desired
C   set of DATA statements should be activated by removing the C from
C   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
C   checked for consistency with the local operating system.
C
C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
C                 a portable library, ACM Transactions on Mathematical
C                 Software 4, 2 (June 1978), pp. 177-188.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   891012  Added VAX G-floating constants.  (WRB)
C   891012  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900618  Added DEC RISC constants.  (WRB)
C   900723  Added IBM RS 6000 constants.  (WRB)
C   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16.
C           (RWC)
C   910710  Added HP 730 constants.  (SMR)
C   911114  Added Convex IEEE constants.  (WRB)
C   920121  Added SUN -r8 compiler option constants.  (WRB)
C   920229  Added Touchstone Delta i860 constants.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   920625  Added Convex -p8 and -pd8 compiler option constants.
C           (BKS, WRB)
C   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
C   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler
C           options.  (DWL, RWC and WRB).
C***END PROLOGUE  I1MACH
C
      INTEGER IMACH(16),OUTPUT
      SAVE IMACH
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT COMPILER
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1022 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        129 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1025 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA IMACH( 1) /          7 /
C     DATA IMACH( 2) /          2 /
C     DATA IMACH( 3) /          2 /
C     DATA IMACH( 4) /          2 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         33 /
C     DATA IMACH( 9) / Z1FFFFFFFF /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -256 /
C     DATA IMACH(13) /        255 /
C     DATA IMACH(14) /         60 /
C     DATA IMACH(15) /       -256 /
C     DATA IMACH(16) /        255 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         48 /
C     DATA IMACH( 6) /          6 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /          8 /
C     DATA IMACH(11) /         13 /
C     DATA IMACH(12) /        -50 /
C     DATA IMACH(13) /         76 /
C     DATA IMACH(14) /         26 /
C     DATA IMACH(15) /        -50 /
C     DATA IMACH(16) /         76 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         48 /
C     DATA IMACH( 6) /          6 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /          8 /
C     DATA IMACH(11) /         13 /
C     DATA IMACH(12) /        -50 /
C     DATA IMACH(13) /         76 /
C     DATA IMACH(14) /         26 /
C     DATA IMACH(15) /     -32754 /
C     DATA IMACH(16) /      32780 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          8 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /      -4095 /
C     DATA IMACH(13) /       4094 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /      -4095 /
C     DATA IMACH(16) /       4094 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /    6LOUTPUT/
C     DATA IMACH( 5) /         60 /
C     DATA IMACH( 6) /         10 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         48 /
C     DATA IMACH( 9) / 00007777777777777777B /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /       -929 /
C     DATA IMACH(13) /       1070 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /       -929 /
C     DATA IMACH(16) /       1069 /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          0 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / Z'7FFFFFFF' /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1022 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fn COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fi COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -p8 COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         53 /
C     DATA IMACH(12) /      -1023 /
C     DATA IMACH(13) /       1023 /
C     DATA IMACH(14) /        113 /
C     DATA IMACH(15) /     -16383 /
C     DATA IMACH(16) /      16383 /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -pd8 COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 9223372036854775807 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         53 /
C     DATA IMACH(12) /      -1023 /
C     DATA IMACH(13) /       1023 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE CRAY
C     USING THE 46 BIT INTEGER COMPILER OPTION
C
C     DATA IMACH( 1) /        100 /
C     DATA IMACH( 2) /        101 /
C     DATA IMACH( 3) /        102 /
C     DATA IMACH( 4) /        101 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          8 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         46 /
C     DATA IMACH( 9) / 1777777777777777B /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /      -8189 /
C     DATA IMACH(13) /       8190 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /      -8099 /
C     DATA IMACH(16) /       8190 /
C
C     MACHINE CONSTANTS FOR THE CRAY
C     USING THE 64 BIT INTEGER COMPILER OPTION
C
C     DATA IMACH( 1) /        100 /
C     DATA IMACH( 2) /        101 /
C     DATA IMACH( 3) /        102 /
C     DATA IMACH( 4) /        101 /
C     DATA IMACH( 5) /         64 /
C     DATA IMACH( 6) /          8 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         63 /
C     DATA IMACH( 9) / 777777777777777777777B /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         47 /
C     DATA IMACH(12) /      -8189 /
C     DATA IMACH(13) /       8190 /
C     DATA IMACH(14) /         94 /
C     DATA IMACH(15) /      -8099 /
C     DATA IMACH(16) /       8190 /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     DATA IMACH( 1) /         11 /
C     DATA IMACH( 2) /         12 /
C     DATA IMACH( 3) /          8 /
C     DATA IMACH( 4) /         10 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /         16 /
C     DATA IMACH(11) /          6 /
C     DATA IMACH(12) /        -64 /
C     DATA IMACH(13) /         63 /
C     DATA IMACH(14) /         14 /
C     DATA IMACH(15) /        -64 /
C     DATA IMACH(16) /         63 /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING G_FLOAT
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING IEEE_FLOAT
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE DEC RISC
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING D_FLOATING
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING G_FLOATING
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1023 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         32 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1022 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          0 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         24 /
C     DATA IMACH( 6) /          3 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         23 /
C     DATA IMACH( 9) /    8388607 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         23 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         38 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /         43 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          6 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         63 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HP 730
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          4 /
C     DATA IMACH( 4) /          1 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         23 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         39 /
C     DATA IMACH(15) /       -128 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          4 /
C     DATA IMACH( 4) /          1 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         23 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         55 /
C     DATA IMACH(15) /       -128 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          7 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         32 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -126 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1015 /
C     DATA IMACH(16) /       1017 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          7 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) /  Z7FFFFFFF /
C     DATA IMACH(10) /         16 /
C     DATA IMACH(11) /          6 /
C     DATA IMACH(12) /        -64 /
C     DATA IMACH(13) /         63 /
C     DATA IMACH(14) /         14 /
C     DATA IMACH(15) /        -64 /
C     DATA IMACH(16) /         63 /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          0 /
C     DATA IMACH( 4) /          0 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          0 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE INTEL i860
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          5 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         54 /
C     DATA IMACH(15) /       -101 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          5 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         62 /
C     DATA IMACH(15) /       -128 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          5 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -125 /
C     DATA IMACH(13) /        128 /
C     DATA IMACH(14) /         53 /
C     DATA IMACH(15) /      -1021 /
C     DATA IMACH(16) /       1024 /
C
C     MACHINE CONSTANTS FOR THE SUN
C     USING THE -r8 COMPILER OPTION
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          6 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         32 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         53 /
C     DATA IMACH(12) /      -1021 /
C     DATA IMACH(13) /       1024 /
C     DATA IMACH(14) /        113 /
C     DATA IMACH(15) /     -16381 /
C     DATA IMACH(16) /      16384 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
C
C     DATA IMACH( 1) /          5 /
C     DATA IMACH( 2) /          6 /
C     DATA IMACH( 3) /          1 /
C     DATA IMACH( 4) /          6 /
C     DATA IMACH( 5) /         36 /
C     DATA IMACH( 6) /          4 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         27 /
C     DATA IMACH(12) /       -128 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         60 /
C     DATA IMACH(15) /      -1024 /
C     DATA IMACH(16) /       1023 /
C
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
C
C     DATA IMACH( 1) /          1 /
C     DATA IMACH( 2) /          1 /
C     DATA IMACH( 3) /          0 /
C     DATA IMACH( 4) /          1 /
C     DATA IMACH( 5) /         16 /
C     DATA IMACH( 6) /          2 /
C     DATA IMACH( 7) /          2 /
C     DATA IMACH( 8) /         15 /
C     DATA IMACH( 9) /      32767 /
C     DATA IMACH(10) /          2 /
C     DATA IMACH(11) /         24 /
C     DATA IMACH(12) /       -127 /
C     DATA IMACH(13) /        127 /
C     DATA IMACH(14) /         56 /
C     DATA IMACH(15) /       -127 /
C     DATA IMACH(16) /        127 /
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10
C
      I1MACH = IMACH(I)
      RETURN
C
   10 CONTINUE
      WRITE (UNIT = OUTPUT, FMT = 9000)
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
C
C     CALL FDUMP
C
      STOP
      END
*DECK J4SAVE
      FUNCTION J4SAVE (IWHICH, IVALUE, ISET)
C***BEGIN PROLOGUE  J4SAVE
C***SUBSIDIARY
C***PURPOSE  Save or recall global variables needed by error
C            handling routines.
C***LIBRARY   SLATEC (XERROR)
C***TYPE      INTEGER (J4SAVE-I)
C***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        J4SAVE saves and recalls several global variables needed
C        by the library error handling routines.
C
C     Description of Parameters
C      --Input--
C        IWHICH - Index of item desired.
C                = 1 Refers to current error number.
C                = 2 Refers to current error control flag.
C                = 3 Refers to current unit number to which error
C                    messages are to be sent.  (0 means use standard.)
C                = 4 Refers to the maximum number of times any
C                     message is to be printed (as set by XERMAX).
C                = 5 Refers to the total number of units to which
C                     each error message is to be written.
C                = 6 Refers to the 2nd unit for error messages
C                = 7 Refers to the 3rd unit for error messages
C                = 8 Refers to the 4th unit for error messages
C                = 9 Refers to the 5th unit for error messages
C        IVALUE - The value to be set for the IWHICH-th parameter,
C                 if ISET is .TRUE. .
C        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
C                 given the value, IVALUE.  If ISET=.FALSE., the
C                 IWHICH-th parameter will be unchanged, and IVALUE
C                 is a dummy parameter.
C      --Output--
C        The (old) value of the IWHICH-th parameter will be returned
C        in the function value, J4SAVE.
C
C***SEE ALSO  XERMSG
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900205  Minor modifications to prologue.  (WRB)
C   900402  Added TYPE section.  (WRB)
C   910411  Added KEYWORDS section.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
C***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
*DECK XERCNT
      SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
C***BEGIN PROLOGUE  XERCNT
C***SUBSIDIARY
C***PURPOSE  Allow user control over handling of errors.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERCNT-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        Allows user control over handling of individual errors.
C        Just after each message is recorded, but before it is
C        processed any further (i.e., before it is printed or
C        a decision to abort is made), a call is made to XERCNT.
C        If the user has provided his own version of XERCNT, he
C        can then override the value of KONTROL used in processing
C        this message by redefining its value.
C        KONTRL may be set to any value from -2 to 2.
C        The meanings for KONTRL are the same as in XSETF, except
C        that the value of KONTRL changes only for this message.
C        If KONTRL is set to a value outside the range from -2 to 2,
C        it will be moved back into that range.
C
C     Description of Parameters
C
C      --Input--
C        LIBRAR - the library that the routine is in.
C        SUBROU - the subroutine that XERMSG is being called from
C        MESSG  - the first 20 characters of the error message.
C        NERR   - same as in the call to XERMSG.
C        LEVEL  - same as in the call to XERMSG.
C        KONTRL - the current value of the control flag as set
C                 by a call to XSETF.
C
C      --Output--
C        KONTRL - the new value of KONTRL.  If KONTRL is not
C                 defined, it will remain at its original value.
C                 This changed value of control affects only
C                 the current occurrence of the current message.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900206  Routine changed from user-callable to subsidiary.  (WRB)
C   900510  Changed calling sequence to include LIBRARY and SUBROUTINE
C           names, changed routine name from XERCTL to XERCNT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERCNT
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
C***FIRST EXECUTABLE STATEMENT  XERCNT
      RETURN
      END
*DECK XERHLT
      SUBROUTINE XERHLT (MESSG)
C***BEGIN PROLOGUE  XERHLT
C***SUBSIDIARY
C***PURPOSE  Abort program execution and print error message.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERHLT-A)
C***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        ***Note*** machine dependent routine
C        XERHLT aborts the execution of the program.
C        The error message causing the abort is given in the calling
C        sequence, in case one needs it for printing on a dayfile,
C        for example.
C
C     Description of Parameters
C        MESSG is as in XERMSG.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900206  Routine changed from user-callable to subsidiary.  (WRB)
C   900510  Changed calling sequence to delete length of character
C           and changed routine name from XERABT to XERHLT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERHLT
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERHLT
      STOP
      END
*DECK XERMSG
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
C***BEGIN PROLOGUE  XERMSG
C***PURPOSE  Process error messages for SLATEC and other libraries.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERMSG-A)
C***KEYWORDS  ERROR MESSAGE, XERROR
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
C***DESCRIPTION
C
C   XERMSG processes a diagnostic message in a manner determined by the
C   value of LEVEL and the current value of the library error control
C   flag, KONTRL.  See subroutine XSETF for details.
C
C    LIBRAR   A character constant (or character variable) with the name
C             of the library.  This will be 'SLATEC' for the SLATEC
C             Common Math Library.  The error handling package is
C             general enough to be used by many libraries
C             simultaneously, so it is desirable for the routine that
C             detects and reports an error to identify the library name
C             as well as the routine name.
C
C    SUBROU   A character constant (or character variable) with the name
C             of the routine that detected the error.  Usually it is the
C             name of the routine that is calling XERMSG.  There are
C             some instances where a user callable library routine calls
C             lower level subsidiary routines where the error is
C             detected.  In such cases it may be more informative to
C             supply the name of the routine the user called rather than
C             the name of the subsidiary routine that detected the
C             error.
C
C    MESSG    A character constant (or character variable) with the text
C             of the error or warning message.  In the example below,
C             the message is a character constant that contains a
C             generic message.
C
C                   CALL XERMSG ('SLATEC', 'MMPY',
C                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
C                  *3, 1)
C
C             It is possible (and is sometimes desirable) to generate a
C             specific message--e.g., one that contains actual numeric
C             values.  Specific numeric values can be converted into
C             character strings using formatted WRITE statements into
C             character variables.  This is called standard Fortran
C             internal file I/O and is exemplified in the first three
C             lines of the following example.  You can also catenate
C             substrings of characters to construct the error message.
C             Here is an example showing the use of both writing to
C             an internal file and catenating character strings.
C
C                   CHARACTER*5 CHARN, CHARL
C                   WRITE (CHARN,10) N
C                   WRITE (CHARL,10) LDA
C                10 FORMAT(I5)
C                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
C                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
C                  *   CHARL, 3, 1)
C
C             There are two subtleties worth mentioning.  One is that
C             the // for character catenation is used to construct the
C             error message so that no single character constant is
C             continued to the next line.  This avoids confusion as to
C             whether there are trailing blanks at the end of the line.
C             The second is that by catenating the parts of the message
C             as an actual argument rather than encoding the entire
C             message into one large character variable, we avoid
C             having to know how long the message will be in order to
C             declare an adequate length for that large character
C             variable.  XERMSG calls XERPRN to print the message using
C             multiple lines if necessary.  If the message is very long,
C             XERPRN will break it into pieces of 72 characters (as
C             requested by XERMSG) for printing on multiple lines.
C             Also, XERMSG asks XERPRN to prefix each line with ' *  '
C             so that the total line length could be 76 characters.
C             Note also that XERPRN scans the error message backwards
C             to ignore trailing blanks.  Another feature is that
C             the substring '$$' is treated as a new line sentinel
C             by XERPRN.  If you want to construct a multiline
C             message without having to count out multiples of 72
C             characters, just use '$$' as a separator.  '$$'
C             obviously must occur within 72 characters of the
C             start of each line to have its intended effect since
C             XERPRN is asked to wrap around at 72 characters in
C             addition to looking for '$$'.
C
C    NERR     An integer value that is chosen by the library routine's
C             author.  It must be in the range -99 to 999 (three
C             printable digits).  Each distinct error should have its
C             own error number.  These error numbers should be described
C             in the machine readable documentation for the routine.
C             The error numbers need be unique only within each routine,
C             so it is reasonable for each routine to start enumerating
C             errors from 1 and proceeding to the next integer.
C
C    LEVEL    An integer value in the range 0 to 2 that indicates the
C             level (severity) of the error.  Their meanings are
C
C            -1  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.  An attempt is made to only print this
C                message once.
C
C             0  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.
C
C             1  A recoverable error.  This is used even if the error is
C                so serious that the routine cannot return any useful
C                answer.  If the user has told the error package to
C                return after recoverable errors, then XERMSG will
C                return to the Library routine which can then return to
C                the user's routine.  The user may also permit the error
C                package to terminate the program upon encountering a
C                recoverable error.
C
C             2  A fatal error.  XERMSG will not return to its caller
C                after it receives a fatal error.  This level should
C                hardly ever be used; it is much better to allow the
C                user a chance to recover.  An example of one of the few
C                cases in which it is permissible to declare a level 2
C                error is a reverse communication Library routine that
C                is likely to be called repeatedly until it integrates
C                across some interval.  If there is a serious error in
C                the input such that another step cannot be taken and
C                the Library routine is called again without the input
C                error having been corrected by the caller, the Library
C                routine will probably be called forever with improper
C                input.  In this case, it is reasonable to declare the
C                error to be fatal.
C
C    Each of the arguments to XERMSG is input; none will be modified by
C    XERMSG.  A routine may make multiple calls to XERMSG with warning
C    level messages; however, after a call to XERMSG with a recoverable
C    error, the routine should return to the user.  Do not try to call
C    XERMSG with a second recoverable error after the first recoverable
C    error because the error package saves the error number.  The user
C    can retrieve this error number by calling another entry point in
C    the error handling package and then clear the error number when
C    recovering from the error.  Calling XERMSG in succession causes the
C    old error number to be overwritten by the latest error number.
C    This is considered harmless for error numbers associated with
C    warning messages but must not be done for error numbers of serious
C    errors.  After a call to XERMSG with a recoverable error, the user
C    must be given a chance to call NUMXER or XERCLR to retrieve or
C    clear the error number.
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
C***REVISION HISTORY  (YYMMDD)
C   880101  DATE WRITTEN
C   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
C           THERE ARE TWO BASIC CHANGES.
C           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
C               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
C               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
C               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
C               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
C               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
C               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
C               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
C           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
C               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
C               OF LOWER CASE.
C   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
C           THE PRINCIPAL CHANGES ARE
C           1.  CLARIFY COMMENTS IN THE PROLOGUES
C           2.  RENAME XRPRNT TO XERPRN
C           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
C               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
C               CHARACTER FOR NEW RECORDS.
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
C           CLEAN UP THE CODING.
C   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
C           PREFIX.
C   891013  REVISED TO CORRECT COMMENTS.
C   891214  Prologue converted to Version 4.0 format.  (WRB)
C   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
C           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
C           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
C           XERCTL to XERCNT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERMSG
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8 XLIBR, XSUBR
      CHARACTER*72  TEMP
      CHARACTER*20  LFIRST
C***FIRST EXECUTABLE STATEMENT  XERMSG
      LKNTRL = J4SAVE (2, 0, .FALSE.)
      MAXMES = J4SAVE (4, 0, .FALSE.)
C
C       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
C       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
C          SHOULD BE PRINTED.
C
C       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
C          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
C          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
C
      IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR.
     *   LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' //
     *      'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '//
     *      'JOB ABORT DUE TO FATAL ERROR.', 72)
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
         CALL XERHLT (' ***XERMSG -- INVALID INPUT')
         RETURN
      ENDIF
C
C       RECORD THE MESSAGE.
C
      I = J4SAVE (1, NERR, .TRUE.)
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
C
C       HANDLE PRINT-ONCE WARNING MESSAGES.
C
      IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN
C
C       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
C
      XLIBR  = LIBRAR
      XSUBR  = SUBROU
      LFIRST = MESSG
      LERR   = NERR
      LLEVEL = LEVEL
      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
C
      LKNTRL = MAX(-2, MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
C
C       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
C       ZERO AND THE ERROR IS NOT FATAL.
C
      IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30
      IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30
      IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30
      IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
C
C       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
C       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
C       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
C       IS NOT ZERO.
C
      IF (LKNTRL .NE. 0) THEN
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '
         I = MIN(LEN(SUBROU), 16)
         TEMP(22:21+I) = SUBROU(1:I)
         TEMP(22+I:33+I) = ' IN LIBRARY '
         LTEMP = 33 + I
         I = MIN(LEN(LIBRAR), 16)
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
         LTEMP = LTEMP + I + 1
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
C       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
C       FROM EACH OF THE FOLLOWING THREE OPTIONS.
C       1.  LEVEL OF THE MESSAGE
C              'INFORMATIVE MESSAGE'
C              'POTENTIALLY RECOVERABLE ERROR'
C              'FATAL ERROR'
C       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
C              'PROG CONTINUES'
C              'PROG ABORTED'
C       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
C           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
C           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
C              'TRACEBACK REQUESTED'
C              'TRACEBACK NOT REQUESTED'
C       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
C       EXCEED 74 CHARACTERS.
C       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
C
      IF (LKNTRL .GT. 0) THEN
C
C       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
C
         IF (LEVEL .LE. 0) THEN
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'
            LTEMP = 20
         ELSEIF (LEVEL .EQ. 1) THEN
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
            LTEMP = 30
         ELSE
            TEMP(1:12) = 'FATAL ERROR,'
            LTEMP = 12
         ENDIF
C
C       THEN WHETHER THE PROGRAM WILL CONTINUE.
C
         IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR.
     *       (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
            LTEMP = LTEMP + 14
         ELSE
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
            LTEMP = LTEMP + 16
         ENDIF
C
C       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
C
         IF (LKNTRL .GT. 0) THEN
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
            LTEMP = LTEMP + 20
         ELSE
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
            LTEMP = LTEMP + 24
         ENDIF
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       NOW SEND OUT THE MESSAGE.
C
      CALL XERPRN (' *  ', -1, MESSG, 72)
C
C       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
C          TRACEBACK.
C
      IF (LKNTRL .GT. 0) THEN
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
         DO 10 I=16,22
            IF (TEMP(I:I) .NE. ' ') GO TO 20
   10    CONTINUE
C
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
         CALL FDUMP
      ENDIF
C
C       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
C
      IF (LKNTRL .NE. 0) THEN
         CALL XERPRN (' *  ', -1, ' ', 72)
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
         CALL XERPRN ('    ',  0, ' ', 72)
      ENDIF
C
C       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
C       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
C
   30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN
C
C       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
C       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
C       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
C
      IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
         IF (LEVEL .EQ. 1) THEN
            CALL XERPRN
     *         (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
         ELSE
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
         ENDIF
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
         CALL XERHLT (' ')
      ELSE
         CALL XERHLT (MESSG)
      ENDIF
      RETURN
      END
*DECK XERPRN
      SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)
C***BEGIN PROLOGUE  XERPRN
C***SUBSIDIARY
C***PURPOSE  Print error messages processed by XERMSG.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERPRN-A)
C***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
C***DESCRIPTION
C
C This routine sends one or more lines to each of the (up to five)
C logical units to which error messages are to be sent.  This routine
C is called several times by XERMSG, sometimes with a single line to
C print and sometimes with a (potentially very long) message that may
C wrap around into multiple lines.
C
C PREFIX  Input argument of type CHARACTER.  This argument contains
C         characters to be put at the beginning of each line before
C         the body of the message.  No more than 16 characters of
C         PREFIX will be used.
C
C NPREF   Input argument of type INTEGER.  This argument is the number
C         of characters to use from PREFIX.  If it is negative, the
C         intrinsic function LEN is used to determine its length.  If
C         it is zero, PREFIX is not used.  If it exceeds 16 or if
C         LEN(PREFIX) exceeds 16, only the first 16 characters will be
C         used.  If NPREF is positive and the length of PREFIX is less
C         than NPREF, a copy of PREFIX extended with blanks to length
C         NPREF will be used.
C
C MESSG   Input argument of type CHARACTER.  This is the text of a
C         message to be printed.  If it is a long message, it will be
C         broken into pieces for printing on multiple lines.  Each line
C         will start with the appropriate prefix and be followed by a
C         piece of the message.  NWRAP is the number of characters per
C         piece; that is, after each NWRAP characters, we break and
C         start a new line.  In addition the characters '$$' embedded
C         in MESSG are a sentinel for a new line.  The counting of
C         characters up to NWRAP starts over for each new line.  The
C         value of NWRAP typically used by XERMSG is 72 since many
C         older error messages in the SLATEC Library are laid out to
C         rely on wrap-around every 72 characters.
C
C NWRAP   Input argument of type INTEGER.  This gives the maximum size
C         piece into which to break MESSG for printing on multiple
C         lines.  An embedded '$$' ends a line, and the count restarts
C         at the following character.  If a line break does not occur
C         on a blank (it would split a word) that word is moved to the
C         next line.  Values of NWRAP less than 16 will be treated as
C         16.  Values of NWRAP greater than 132 will be treated as 132.
C         The actual line length will be NPREF + NWRAP after NPREF has
C         been adjusted to fall between 0 and 16 and NWRAP has been
C         adjusted to fall between 16 and 132.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  I1MACH, XGETUA
C***REVISION HISTORY  (YYMMDD)
C   880621  DATE WRITTEN
C   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
C           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
C           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
C           SLASH CHARACTER IN FORMAT STATEMENTS.
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
C           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
C           LINES TO BE PRINTED.
C   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
C           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
C   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
C   891214  Prologue converted to Version 4.0 format.  (WRB)
C   900510  Added code to break messages between words.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERPRN
      CHARACTER*(*) PREFIX, MESSG
      INTEGER NPREF, NWRAP
      CHARACTER*148 CBUFF
      INTEGER IU(5), NUNIT
      CHARACTER*2 NEWLIN
      PARAMETER (NEWLIN = '$$')
C***FIRST EXECUTABLE STATEMENT  XERPRN
      CALL XGETUA(IU,NUNIT)
C
C       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
C       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
C       ERROR MESSAGE UNIT.
C
      N = I1MACH(4)
      DO 10 I=1,NUNIT
         IF (IU(I) .EQ. 0) IU(I) = N
   10 CONTINUE
C
C       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
C       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
C       THE REST OF THIS ROUTINE.
C
      IF ( NPREF .LT. 0 ) THEN
         LPREF = LEN(PREFIX)
      ELSE
         LPREF = NPREF
      ENDIF
      LPREF = MIN(16, LPREF)
      IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX
C
C       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
C       TIME FROM MESSG TO PRINT ON ONE LINE.
C
      LWRAP = MAX(16, MIN(132, NWRAP))
C
C       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
C
      LENMSG = LEN(MESSG)
      N = LENMSG
      DO 20 I=1,N
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30
         LENMSG = LENMSG - 1
   20 CONTINUE
   30 CONTINUE
C
C       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
C
      IF (LENMSG .EQ. 0) THEN
         CBUFF(LPREF+1:LPREF+1) = ' '
         DO 40 I=1,NUNIT
            WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
   40    CONTINUE
         RETURN
      ENDIF
C
C       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
C       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
C       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
C       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
C
C       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
C       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
C       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
C       OF THE SECOND ARGUMENT.
C
C       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
C       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
C       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
C       POSITION NEXTC.
C
C       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
C                       REMAINDER OF THE CHARACTER STRING.  LPIECE
C                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
C                       WHICHEVER IS LESS.
C
C       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
C                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
C                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
C                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
C                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
C                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
C                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
C                       SHOULD BE INCREMENTED BY 2.
C
C       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
C
C       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
C                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
C                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
C                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
C                       AT THE END OF A LINE.
C
      NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
      IF (LPIECE .EQ. 0) THEN
C
C       THERE WAS NO NEW LINE SENTINEL FOUND.
C
         IDELTA = 0
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
         IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
            DO 52 I=LPIECE+1,2,-1
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                  LPIECE = I-1
                  IDELTA = 1
                  GOTO 54
               ENDIF
   52       CONTINUE
         ENDIF
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSEIF (LPIECE .EQ. 1) THEN
C
C       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
C       DON'T PRINT A BLANK LINE.
C
         NEXTC = NEXTC + 2
         GO TO 50
      ELSEIF (LPIECE .GT. LWRAP+1) THEN
C
C       LPIECE SHOULD BE SET DOWN TO LWRAP.
C
         IDELTA = 0
         LPIECE = LWRAP
         DO 56 I=LPIECE+1,2,-1
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
               LPIECE = I-1
               IDELTA = 1
               GOTO 58
            ENDIF
   56    CONTINUE
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSE
C
C       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
C       WE SHOULD DECREMENT LPIECE BY ONE.
C
         LPIECE = LPIECE - 1
         CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC  = NEXTC + LPIECE + 2
      ENDIF
C
C       PRINT
C
      DO 60 I=1,NUNIT
         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
   60 CONTINUE
C
      IF (NEXTC .LE. LENMSG) GO TO 50
      RETURN
      END
*DECK XERSVE
      SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL,
     +   ICOUNT)
C***BEGIN PROLOGUE  XERSVE
C***SUBSIDIARY
C***PURPOSE  Record that an error has occurred.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3
C***TYPE      ALL (XERSVE-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C *Usage:
C
C        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
C        CHARACTER * (len) LIBRAR, SUBROU, MESSG
C
C        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
C
C *Arguments:
C
C        LIBRAR :IN    is the library that the message is from.
C        SUBROU :IN    is the subroutine that the message is from.
C        MESSG  :IN    is the message to be saved.
C        KFLAG  :IN    indicates the action to be performed.
C                      when KFLAG > 0, the message in MESSG is saved.
C                      when KFLAG=0 the tables will be dumped and
C                      cleared.
C                      when KFLAG < 0, the tables will be dumped and
C                      not cleared.
C        NERR   :IN    is the error number.
C        LEVEL  :IN    is the error severity.
C        ICOUNT :OUT   the number of times this message has been seen,
C                      or zero if the table has overflowed and does not
C                      contain this message specifically.  When KFLAG=0,
C                      ICOUNT will not be altered.
C
C *Description:
C
C   Record that this error occurred and possibly dump and clear the
C   tables.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  I1MACH, XGETUA
C***REVISION HISTORY  (YYMMDD)
C   800319  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900413  Routine modified to remove reference to KFLAG.  (WRB)
C   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
C           sequence, use IF-THEN-ELSE, make number of saved entries
C           easily changeable, changed routine name from XERSAV to
C           XERSVE.  (RWC)
C   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERSVE
      PARAMETER (LENTAB=10)
      INTEGER LUN(5)
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
      CHARACTER*20 MESTAB(LENTAB), MES
      DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
      SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
      DATA KOUNTX/0/, NMSG/0/
C***FIRST EXECUTABLE STATEMENT  XERSVE
C
      IF (KFLAG.LE.0) THEN
C
C        Dump the table.
C
         IF (NMSG.EQ.0) RETURN
C
C        Print to each unit.
C
         CALL XGETUA (LUN, NUNIT)
         DO 20 KUNIT = 1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C
C           Print the table header.
C
            WRITE (IUNIT,9000)
C
C           Print body of table.
C
            DO 10 I = 1,NMSG
               WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I),
     *            NERTAB(I),LEVTAB(I),KOUNT(I)
   10       CONTINUE
C
C           Print number of other errors.
C
            IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX
            WRITE (IUNIT,9030)
   20    CONTINUE
C
C        Clear the error tables.
C
         IF (KFLAG.EQ.0) THEN
            NMSG = 0
            KOUNTX = 0
         ENDIF
      ELSE
C
C        PROCESS A MESSAGE...
C        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
C
         LIB = LIBRAR
         SUB = SUBROU
         MES = MESSG
         DO 30 I = 1,NMSG
            IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND.
     *         MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND.
     *         LEVEL.EQ.LEVTAB(I)) THEN
                  KOUNT(I) = KOUNT(I) + 1
                  ICOUNT = KOUNT(I)
                  RETURN
            ENDIF
   30    CONTINUE
C
         IF (NMSG.LT.LENTAB) THEN
C
C           Empty slot found for new message.
C
            NMSG = NMSG + 1
            LIBTAB(I) = LIB
            SUBTAB(I) = SUB
            MESTAB(I) = MES
            NERTAB(I) = NERR
            LEVTAB(I) = LEVEL
            KOUNT (I) = 1
            ICOUNT    = 1
         ELSE
C
C           Table is full.
C
            KOUNTX = KOUNTX+1
            ICOUNT = 0
         ENDIF
      ENDIF
      RETURN
C
C     Formats.
C
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' /
     +   ' LIBRARY    SUBROUTINE MESSAGE START             NERR',
     +   '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
      END
*DECK XGETUA
      SUBROUTINE XGETUA (IUNITA, N)
C***BEGIN PROLOGUE  XGETUA
C***PURPOSE  Return unit number(s) to which error messages are being
C            sent.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XGETUA-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        XGETUA may be called to determine the unit number or numbers
C        to which error messages are being sent.
C        These unit numbers may have been set by a call to XSETUN,
C        or a call to XSETUA, or may be a default value.
C
C     Description of Parameters
C      --Output--
C        IUNIT - an array of one to five unit numbers, depending
C                on the value of N.  A value of zero refers to the
C                default unit, as defined by the I1MACH machine
C                constant routine.  Only IUNIT(1),...,IUNIT(N) are
C                defined by XGETUA.  The values of IUNIT(N+1),...,
C                IUNIT(5) are not defined (for N .LT. 5) or altered
C                in any way by XGETUA.
C        N     - the number of units to which copies of the
C                error messages are being sent.  N will be in the
C                range from 1 to 5.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  J4SAVE
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
C***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
*DECK BVALU
      FUNCTION BVALU (T, A, N, K, IDERIV, X, INBV, WORK)
C***BEGIN PROLOGUE  BVALU
C***PURPOSE  Evaluate the B-representation of a B-spline at X for the
C            function value or any of its derivatives.
C***LIBRARY   SLATEC
C***CATEGORY  E3, K6
C***TYPE      SINGLE PRECISION (BVALU-S, DBVALU-D)
C***KEYWORDS  DIFFERENTIATION OF B-SPLINE, EVALUATION OF B-SPLINE
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Abstract
C         BVALU is the BVALUE function of the reference.
C
C         BVALU evaluates the B-representation (T,A,N,K) of a B-spline
C         at X for the function value on IDERIV = 0 or any of its
C         derivatives on IDERIV = 1,2,...,K-1.  Right limiting values
C         (right derivatives) are returned except at the right end
C         point X=T(N+1) where left limiting values are computed.  The
C         spline is defined on T(K) .LE. X .LE. T(N+1).  BVALU returns
C         a fatal error message when X is outside of this interval.
C
C         To compute left derivatives or left limiting values at a
C         knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
C
C         BVALU calls INTRV
C
C     Description of Arguments
C         Input
C          T       - knot vector of length N+K
C          A       - B-spline coefficient vector of length N
C          N       - number of B-spline coefficients
C                    N = sum of knot multiplicities-K
C          K       - order of the B-spline, K .GE. 1
C          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1
C                    IDERIV=0 returns the B-spline value
C          X       - argument, T(K) .LE. X .LE. T(N+1)
C          INBV    - an initialization parameter which must be set
C                    to 1 the first time BVALU is called.
C
C         Output
C          INBV    - INBV contains information for efficient process-
C                    ing after the initial call and INBV must not
C                    be changed by the user.  Distinct splines require
C                    distinct INBV parameters.
C          WORK    - work vector of length 3*K.
C          BVALU   - value of the IDERIV-th derivative at X
C
C     Error Conditions
C         An improper input is a fatal error
C
C***REFERENCES  Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C***ROUTINES CALLED  INTRV, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  BVALU
C
      INTEGER I,IDERIV,IDERP1,IHI,IHMKMJ,ILO,IMK,IMKPJ, INBV, IPJ,
     1 IP1, IP1MJ, J, JJ, J1, J2, K, KMIDER, KMJ, KM1, KPK, MFLAG, N
      REAL A, FKMJ, T, WORK, X
C     DIMENSION T(N+K), WORK(3*K)
      DIMENSION T(*), A(*), WORK(*)
C***FIRST EXECUTABLE STATEMENT  BVALU
      BVALU = 0.0E0
      IF(K.LT.1) GO TO 102
      IF(N.LT.K) GO TO 101
      IF(IDERIV.LT.0 .OR. IDERIV.GE.K) GO TO 110
      KMIDER = K - IDERIV
C
C *** FIND *I* IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)
C     (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).
      KM1 = K - 1
      CALL INTRV(T, N+1, X, INBV, I, MFLAG)
      IF (X.LT.T(K)) GO TO 120
      IF (MFLAG.EQ.0) GO TO 20
      IF (X.GT.T(I)) GO TO 130
   10 IF (I.EQ.K) GO TO 140
      I = I - 1
      IF (X.EQ.T(I)) GO TO 10
C
C *** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES
C     WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K
C
   20 IMK = I - K
      DO 30 J=1,K
        IMKPJ = IMK + J
        WORK(J) = A(IMKPJ)
   30 CONTINUE
      IF (IDERIV.EQ.0) GO TO 60
      DO 50 J=1,IDERIV
        KMJ = K - J
        FKMJ = KMJ
        DO 40 JJ=1,KMJ
          IHI = I + JJ
          IHMKMJ = IHI - KMJ
          WORK(JJ) = (WORK(JJ+1)-WORK(JJ))/(T(IHI)-T(IHMKMJ))*FKMJ
   40   CONTINUE
   50 CONTINUE
C
C *** COMPUTE VALUE AT *X* IN (T(I),(T(I+1)) OF IDERIV-TH DERIVATIVE,
C     GIVEN ITS RELEVANT B-SPLINE COEFF. IN AJ(1),...,AJ(K-IDERIV).
   60 IF (IDERIV.EQ.KM1) GO TO 100
      IP1 = I + 1
      KPK = K + K
      J1 = K + 1
      J2 = KPK + 1
      DO 70 J=1,KMIDER
        IPJ = I + J
        WORK(J1) = T(IPJ) - X
        IP1MJ = IP1 - J
        WORK(J2) = X - T(IP1MJ)
        J1 = J1 + 1
        J2 = J2 + 1
   70 CONTINUE
      IDERP1 = IDERIV + 1
      DO 90 J=IDERP1,KM1
        KMJ = K - J
        ILO = KMJ
        DO 80 JJ=1,KMJ
          WORK(JJ) = (WORK(JJ+1)*WORK(KPK+ILO)+WORK(JJ)
     1              *WORK(K+JJ))/(WORK(KPK+ILO)+WORK(K+JJ))
          ILO = ILO - 1
   80   CONTINUE
   90 CONTINUE
  100 BVALU = WORK(1)
      RETURN
C
C
  101 CONTINUE
      CALL XERMSG ('SLATEC', 'BVALU', 'N DOES NOT SATISFY N.GE.K', 2,
     +   1)
      RETURN
  102 CONTINUE
      CALL XERMSG ('SLATEC', 'BVALU', 'K DOES NOT SATISFY K.GE.1', 2,
     +   1)
      RETURN
  110 CONTINUE
      CALL XERMSG ('SLATEC', 'BVALU',
     +   'IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K', 2, 1)
      RETURN
  120 CONTINUE
      CALL XERMSG ('SLATEC', 'BVALU',
     +   'X IS N0T GREATER THAN OR EQUAL TO T(K)', 2, 1)
      RETURN
  130 CONTINUE
      CALL XERMSG ('SLATEC', 'BVALU',
     +   'X IS NOT LESS THAN OR EQUAL TO T(N+1)', 2, 1)
      RETURN
  140 CONTINUE
      CALL XERMSG ('SLATEC', 'BVALU',
     +   'A LEFT LIMITING VALUE CANNOT BE OBTAINED AT T(K)', 2, 1)
      RETURN
      END
*DECK BSQAD
      SUBROUTINE BSQAD (T, BCOEF, N, K, X1, X2, BQUAD, WORK)
C***BEGIN PROLOGUE  BSQAD
C***PURPOSE  Compute the integral of a K-th order B-spline using the
C            B-representation.
C***LIBRARY   SLATEC
C***CATEGORY  H2A2A1, E3, K6
C***TYPE      SINGLE PRECISION (BSQAD-S, DBSQAD-D)
C***KEYWORDS  INTEGRAL OF B-SPLINES, QUADRATURE
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C         BSQAD computes the integral on (X1,X2) of a K-th order
C         B-spline using the B-representation (T,BCOEF,N,K).  Orders
C         K as high as 20 are permitted by applying a 2, 6, or 10
C         point Gauss formula on subintervals of (X1,X2) which are
C         formed by included (distinct) knots.
C
C         If orders K greater than 20 are needed, use BFQAD with
C         F(X) = 1.
C
C     Description of Arguments
C         Input
C           T      - knot array of length N+K
C           BCOEF  - B-spline coefficient array of length N
C           N      - length of coefficient array
C           K      - order of B-spline, 1 .LE. K .LE. 20
C           X1,X2  - end points of quadrature interval in
C                    T(K) .LE. X .LE. T(N+1)
C
C         Output
C           BQUAD  - integral of the B-spline over (X1,X2)
C           WORK   - work vector of length 3*K
C
C     Error Conditions
C         Improper input is a fatal error
C
C***REFERENCES  D. E. Amos, Quadrature subroutines for splines and
C                 B-splines, Report SAND79-1825, Sandia Laboratories,
C                 December 1979.
C***ROUTINES CALLED  BVALU, INTRV, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  BSQAD
C
      INTEGER I,IL1,IL2,ILO,INBV, JF,K,LEFT,M,MF,MFLAG,N, NPK, NP1
      REAL A, AA, B, BB, BCOEF, BMA, BPA, BQUAD, C1, GPTS, GWTS, GX, Q,
     1 SUM, T, TA, TB, WORK, X1, X2, Y1, Y2
      REAL BVALU
      DIMENSION T(*), BCOEF(*), GPTS(9), GWTS(9), SUM(5), WORK(*)
C
      SAVE GPTS, GWTS
      DATA GPTS(1), GPTS(2), GPTS(3), GPTS(4), GPTS(5), GPTS(6),
     1     GPTS(7), GPTS(8), GPTS(9)/
     2     5.77350269189625764E-01,     2.38619186083196909E-01,
     3     6.61209386466264514E-01,     9.32469514203152028E-01,
     4     1.48874338981631211E-01,     4.33395394129247191E-01,
     5     6.79409568299024406E-01,     8.65063366688984511E-01,
     6     9.73906528517171720E-01/
      DATA GWTS(1), GWTS(2), GWTS(3), GWTS(4), GWTS(5), GWTS(6),
     1     GWTS(7), GWTS(8), GWTS(9)/
     2     1.00000000000000000E+00,     4.67913934572691047E-01,
     3     3.60761573048138608E-01,     1.71324492379170345E-01,
     4     2.95524224714752870E-01,     2.69266719309996355E-01,
     5     2.19086362515982044E-01,     1.49451349150580593E-01,
     6     6.66713443086881376E-02/
C
C***FIRST EXECUTABLE STATEMENT  BSQAD
      BQUAD = 0.0E0
      IF(K.LT.1 .OR. K.GT.20) GO TO 65
      IF(N.LT.K) GO TO 70
      AA = MIN(X1,X2)
      BB = MAX(X1,X2)
      IF (AA.LT.T(K)) GO TO 60
      NP1 = N + 1
      IF (BB.GT.T(NP1)) GO TO 60
      IF (AA.EQ.BB) RETURN
      NPK = N + K
C     SELECTION OF 2, 6, OR 10 POINT GAUSS FORMULA
      JF = 0
      MF = 1
      IF (K.LE.4) GO TO 10
      JF = 1
      MF = 3
      IF (K.LE.12) GO TO 10
      JF = 4
      MF = 5
   10 CONTINUE
C
      DO 20 I=1,MF
        SUM(I) = 0.0E0
   20 CONTINUE
      ILO = 1
      INBV = 1
      CALL INTRV(T, NPK, AA, ILO, IL1, MFLAG)
      CALL INTRV(T, NPK, BB, ILO, IL2, MFLAG)
      IF (IL2.GE.NP1) IL2 = N
      DO 40 LEFT=IL1,IL2
        TA = T(LEFT)
        TB = T(LEFT+1)
        IF (TA.EQ.TB) GO TO 40
        A = MAX(AA,TA)
        B = MIN(BB,TB)
        BMA = 0.5E0*(B-A)
        BPA = 0.5E0*(B+A)
        DO 30 M=1,MF
          C1 = BMA*GPTS(JF+M)
          GX = -C1 + BPA
          Y2 = BVALU(T,BCOEF,N,K,0,GX,INBV,WORK)
          GX = C1 + BPA
          Y1 = BVALU(T,BCOEF,N,K,0,GX,INBV,WORK)
          SUM(M) = SUM(M) + (Y1+Y2)*BMA
   30   CONTINUE
   40 CONTINUE
      Q = 0.0E0
      DO 50 M=1,MF
        Q = Q + GWTS(JF+M)*SUM(M)
   50 CONTINUE
      IF (X1.GT.X2) Q = -Q
      BQUAD = Q
      RETURN
C
C
   60 CONTINUE
      CALL XERMSG ('SLATEC', 'BSQAD',
     +   'X1 OR X2 OR BOTH DO NOT SATISFY T(K).LE.X.LE.T(N+1)', 2, 1)
      RETURN
   65 CONTINUE
      CALL XERMSG ('SLATEC', 'BSQAD', 'K DOES NOT SATISFY 1.LE.K.LE.20'
     +   , 2, 1)
      RETURN
   70 CONTINUE
      CALL XERMSG ('SLATEC', 'BSQAD', 'N DOES NOT SATISFY N.GE.K', 2,
     +   1)
      RETURN
      END
*DECK INTRV
      SUBROUTINE INTRV (XT, LXT, X, ILO, ILEFT, MFLAG)
C***BEGIN PROLOGUE  INTRV
C***PURPOSE  Compute the largest integer ILEFT in 1 .LE. ILEFT .LE. LXT
C            such that XT(ILEFT) .LE. X where XT(*) is a subdivision
C            of the X interval.
C***LIBRARY   SLATEC
C***CATEGORY  E3, K6
C***TYPE      SINGLE PRECISION (INTRV-S, DINTRV-D)
C***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Abstract
C         INTRV is the INTERV routine of the reference.
C
C         INTRV computes the largest integer ILEFT in 1 .LE. ILEFT .LE.
C         LXT such that XT(ILEFT) .LE. X where XT(*) is a subdivision of
C         the X interval.  Precisely,
C
C                      X .LT. XT(1)                1         -1
C         if  XT(I) .LE. X .LT. XT(I+1)  then  ILEFT=I  , MFLAG=0
C           XT(LXT) .LE. X                         LXT        1,
C
C         That is, when multiplicities are present in the break point
C         to the left of X, the largest index is taken for ILEFT.
C
C     Description of Arguments
C         Input
C          XT      - XT is a knot or break point vector of length LXT
C          LXT     - length of the XT vector
C          X       - argument
C          ILO     - an initialization parameter which must be set
C                    to 1 the first time the spline array XT is
C                    processed by INTRV.
C
C         Output
C          ILO     - ILO contains information for efficient process-
C                    ing after the initial call, and ILO must not be
C                    changed by the user.  Distinct splines require
C                    distinct ILO parameters.
C          ILEFT   - largest integer satisfying XT(ILEFT) .LE. X
C          MFLAG   - signals when X lies out of bounds
C
C     Error Conditions
C         None
C
C***REFERENCES  Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  INTRV
C
      INTEGER IHI, ILEFT, ILO, ISTEP, LXT, MFLAG, MIDDLE
      REAL X, XT
      DIMENSION XT(*)
C***FIRST EXECUTABLE STATEMENT  INTRV
      IHI = ILO + 1
      IF (IHI.LT.LXT) GO TO 10
      IF (X.GE.XT(LXT)) GO TO 110
      IF (LXT.LE.1) GO TO 90
      ILO = LXT - 1
      IHI = LXT
C
   10 IF (X.GE.XT(IHI)) GO TO 40
      IF (X.GE.XT(ILO)) GO TO 100
C
C *** NOW X .LT. XT(IHI) . FIND LOWER BOUND
      ISTEP = 1
   20 IHI = ILO
      ILO = IHI - ISTEP
      IF (ILO.LE.1) GO TO 30
      IF (X.GE.XT(ILO)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 20
   30 ILO = 1
      IF (X.LT.XT(1)) GO TO 90
      GO TO 70
C *** NOW X .GE. XT(ILO) . FIND UPPER BOUND
   40 ISTEP = 1
   50 ILO = IHI
      IHI = ILO + ISTEP
      IF (IHI.GE.LXT) GO TO 60
      IF (X.LT.XT(IHI)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 50
   60 IF (X.GE.XT(LXT)) GO TO 110
      IHI = LXT
C
C *** NOW XT(ILO) .LE. X .LT. XT(IHI) . NARROW THE INTERVAL
   70 MIDDLE = (ILO+IHI)/2
      IF (MIDDLE.EQ.ILO) GO TO 100
C     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1
      IF (X.LT.XT(MIDDLE)) GO TO 80
      ILO = MIDDLE
      GO TO 70
   80 IHI = MIDDLE
      GO TO 70
C *** SET OUTPUT AND RETURN
   90 MFLAG = -1
      ILEFT = 1
      RETURN
  100 MFLAG = 0
      ILEFT = ILO
      RETURN
  110 MFLAG = 1
      ILEFT = LXT
      RETURN
      END

*DECK GAUS8
      SUBROUTINE GAUS8 (FUN, A, B, ERR, ANS, IERR)
C***BEGIN PROLOGUE  GAUS8
C***PURPOSE  Integrate a real function of one variable over a finite
C            interval using an adaptive 8-point Legendre-Gauss
C            algorithm.  Intended primarily for high accuracy
C            integration or integration of smooth functions.
C***LIBRARY   SLATEC
C***CATEGORY  H2A1A1
C***TYPE      SINGLE PRECISION (GAUS8-S, DGAUS8-D)
C***KEYWORDS  ADAPTIVE QUADRATURE, AUTOMATIC INTEGRATOR,
C             GAUSS QUADRATURE, NUMERICAL INTEGRATION
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        GAUS8 integrates real functions of one variable over finite
C        intervals using an adaptive 8-point Legendre-Gauss algorithm.
C        GAUS8 is intended primarily for high accuracy integration
C        or integration of smooth functions.
C
C     Description of Arguments
C
C        Input--
C        FUN - name of external function to be integrated.  This name
C              must be in an EXTERNAL statement in the calling program.
C              FUN must be a REAL function of one REAL argument.  The
C              value of the argument to FUN is the variable of
C              integration which ranges from A to B.
C        A   - lower limit of integration
C        B   - upper limit of integration (may be less than A)
C        ERR - is a requested pseudorelative error tolerance.  Normally
C              pick a value of ABS(ERR) so that STOL .LT. ABS(ERR) .LE.
C              1.0E-3 where STOL is the single precision unit roundoff
C              R1MACH(4).  ANS will normally have no more error than
C              ABS(ERR) times the integral of the absolute value of
C              FUN(X).  Usually, smaller values for ERR yield more
C              accuracy and require more function evaluations.
C
C              A negative value for ERR causes an estimate of the
C              absolute error in ANS to be returned in ERR.  Note that
C              ERR must be a variable (not a constant) in this case.
C              Note also that the user must reset the value of ERR
C              before making any more calls that use the variable ERR.
C
C        Output--
C        ERR - will be an estimate of the absolute error in ANS if the
C              input value of ERR was negative.  (ERR is unchanged if
C              the input value of ERR was non-negative.)  The estimated
C              error is solely for information to the user and should
C              not be used as a correction to the computed integral.
C        ANS - computed value of integral
C        IERR- a status code
C            --Normal codes
C               1 ANS most likely meets requested error tolerance,
C                 or A=B.
C              -1 A and B are too nearly equal to allow normal
C                 integration.  ANS is set to zero.
C            --Abnormal code
C               2 ANS probably does not meet requested error tolerance.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  I1MACH, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   810223  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  GAUS8
      INTEGER IERR, K, KML, KMX, L, LMN, LMX, LR, MXL, NBITS,
     1 NIB, NLMN, NLMX
      INTEGER I1MACH
      REAL A, AA, AE, ANIB, ANS, AREA, B, C, CE, EE, EF, EPS, ERR, EST,
     1 GL, GLR, GR, HH, SQ2, TOL, VL, VR, W1, W2, W3, W4, X1, X2, X3,
     2 X4, X, H
      REAL R1MACH, G8, FUN
      DIMENSION AA(30), HH(30), LR(30), VL(30), GR(30)
      SAVE X1, X2, X3, X4, W1, W2, W3, W4, SQ2,
     1 NLMN, KMX, KML
      DATA X1, X2, X3, X4/
     1     1.83434642495649805E-01,     5.25532409916328986E-01,
     2     7.96666477413626740E-01,     9.60289856497536232E-01/
      DATA W1, W2, W3, W4/
     1     3.62683783378361983E-01,     3.13706645877887287E-01,
     2     2.22381034453374471E-01,     1.01228536290376259E-01/
      DATA SQ2/1.41421356E0/
      DATA NLMN/1/,KMX/5000/,KML/6/
      G8(X,H)=H*((W1*(FUN(X-X1*H) + FUN(X+X1*H))
     1           +W2*(FUN(X-X2*H) + FUN(X+X2*H)))
     2          +(W3*(FUN(X-X3*H) + FUN(X+X3*H))
     3           +W4*(FUN(X-X4*H) + FUN(X+X4*H))))
C***FIRST EXECUTABLE STATEMENT  GAUS8
C
C     Initialize
C
      K = I1MACH(11)
      ANIB = R1MACH(5)*K/0.30102000E0
      NBITS = ANIB
      NLMX = MIN(30,(NBITS*5)/8)
      ANS = 0.0E0
      IERR = 1
      CE = 0.0E0
      IF (A .EQ. B) GO TO 140
      LMX = NLMX
      LMN = NLMN
      IF (B .EQ. 0.0E0) GO TO 10
      IF (SIGN(1.0E0,B)*A .LE. 0.0E0) GO TO 10
      C = ABS(1.0E0-A/B)
      IF (C .GT. 0.1E0) GO TO 10
      IF (C .LE. 0.0E0) GO TO 140
      ANIB = 0.5E0 - LOG(C)/0.69314718E0
      NIB = ANIB
      LMX = MIN(NLMX,NBITS-NIB-7)
      IF (LMX .LT. 1) GO TO 130
      LMN = MIN(LMN,LMX)
   10 TOL = MAX(ABS(ERR),2.0E0**(5-NBITS))/2.0E0
      IF (ERR .EQ. 0.0E0) TOL = SQRT(R1MACH(4))
      EPS = TOL
      HH(1) = (B-A)/4.0E0
      AA(1) = A
      LR(1) = 1
      L = 1
      EST = G8(AA(L)+2.0E0*HH(L),2.0E0*HH(L))
      K = 8
      AREA = ABS(EST)
      EF = 0.5E0
      MXL = 0
C
C     Compute refined estimates, estimate the error, etc.
C
   20 GL = G8(AA(L)+HH(L),HH(L))
      GR(L) = G8(AA(L)+3.0E0*HH(L),HH(L))
      K = K + 16
      AREA = AREA + (ABS(GL)+ABS(GR(L))-ABS(EST))
C     IF (L .LT. LMN) GO TO 11
      GLR = GL + GR(L)
      EE = ABS(EST-GLR)*EF
      AE = MAX(EPS*AREA,TOL*ABS(GLR))
      IF (EE-AE) 40, 40, 50
   30 MXL = 1
   40 CE = CE + (EST-GLR)
      IF (LR(L)) 60, 60, 80
C
C     Consider the left half of this level
C
   50 IF (K .GT. KMX) LMX = KML
      IF (L .GE. LMX) GO TO 30
      L = L + 1
      EPS = EPS*0.5E0
      EF = EF/SQ2
      HH(L) = HH(L-1)*0.5E0
      LR(L) = -1
      AA(L) = AA(L-1)
      EST = GL
      GO TO 20
C
C     Proceed to right half at this level
C
   60 VL(L) = GLR
   70 EST = GR(L-1)
      LR(L) = 1
      AA(L) = AA(L) + 4.0E0*HH(L)
      GO TO 20
C
C     Return one level
C
   80 VR = GLR
   90 IF (L .LE. 1) GO TO 120
      L = L - 1
      EPS = EPS*2.0E0
      EF = EF*SQ2
      IF (LR(L)) 100, 100, 110
  100 VL(L) = VL(L+1) + VR
      GO TO 70
  110 VR = VL(L+1) + VR
      GO TO 90
C
C     Exit
C
  120 ANS = VR
      IF ((MXL.EQ.0) .OR. (ABS(CE).LE.2.0E0*TOL*AREA)) GO TO 140
      IERR = 2
      CALL XERMSG ('SLATEC', 'GAUS8',
     +   'ANS is probably insufficiently accurate.', 3, 1)
      GO TO 140
  130 IERR = -1
      CALL XERMSG ('SLATEC', 'GAUS8',
     +   'A and B are too nearly equal to allow normal integration. $$'
     +   // 'ANS is set to zero and IERR to -1.', 1, -1)
  140 IF (ERR .LT. 0.0E0) ERR = CE
      RETURN
      END
*DECK R1MACH
      REAL FUNCTION R1MACH (I)
C***BEGIN PROLOGUE  R1MACH
C***PURPOSE  Return floating point machine dependent constants.
C***LIBRARY   SLATEC
C***CATEGORY  R1
C***TYPE      SINGLE PRECISION (R1MACH-S, D1MACH-D)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Fox, P. A., (Bell Labs)
C           Hall, A. D., (Bell Labs)
C           Schryer, N. L., (Bell Labs)
C***DESCRIPTION
C
C   R1MACH can be used to obtain machine-dependent parameters for the
C   local machine environment.  It is a function subprogram with one
C   (input) argument, and can be referenced as follows:
C
C        A = R1MACH(I)
C
C   where I=1,...,5.  The (output) value of A above is determined by
C   the (input) value of I.  The results for various values of I are
C   discussed below.
C
C   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
C   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C   R1MACH(3) = B**(-T), the smallest relative spacing.
C   R1MACH(4) = B**(1-T), the largest relative spacing.
C   R1MACH(5) = LOG10(B)
C
C   Assume single precision numbers are represented in the T-digit,
C   base-B form
C
C              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
C   EMIN .LE. E .LE. EMAX.
C
C   The values of B, T, EMIN and EMAX are provided in I1MACH as
C   follows:
C   I1MACH(10) = B, the base.
C   I1MACH(11) = T, the number of base-B digits.
C   I1MACH(12) = EMIN, the smallest exponent E.
C   I1MACH(13) = EMAX, the largest exponent E.
C
C   To alter this function for a particular environment, the desired
C   set of DATA statements should be activated by removing the C from
C   column 1.  Also, the values of R1MACH(1) - R1MACH(4) should be
C   checked for consistency with the local operating system.
C
C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
C                 a portable library, ACM Transactions on Mathematical
C                 Software 4, 2 (June 1978), pp. 177-188.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790101  DATE WRITTEN
C   890213  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900618  Added DEC RISC constants.  (WRB)
C   900723  Added IBM RS 6000 constants.  (WRB)
C   910710  Added HP 730 constants.  (SMR)
C   911114  Added Convex IEEE constants.  (WRB)
C   920121  Added SUN -r8 compiler option constants.  (WRB)
C   920229  Added Touchstone Delta i860 constants.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   920625  Added CONVEX -p8 and -pd8 compiler option constants.
C           (BKS, WRB)
C   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
C***END PROLOGUE  R1MACH
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
C
      REAL RMACH(5)
      SAVE RMACH
C
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7F7FFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7EFFFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA SMALL(1) / 16#00800000 /
C     DATA LARGE(1) / 16#7FFFFFFF /
C     DATA RIGHT(1) / 16#33800000 /
C     DATA DIVER(1) / 16#34000000 /
C     DATA LOG10(1) / 16#3E9A209B /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA RMACH(1) / Z400800000 /
C     DATA RMACH(2) / Z5FFFFFFFF /
C     DATA RMACH(3) / Z4E9800000 /
C     DATA RMACH(4) / Z4EA800000 /
C     DATA RMACH(5) / Z500E730E8 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS
C
C     DATA RMACH(1) / O1771000000000000 /
C     DATA RMACH(2) / O0777777777777777 /
C     DATA RMACH(3) / O1311000000000000 /
C     DATA RMACH(4) / O1301000000000000 /
C     DATA RMACH(5) / O1157163034761675 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA RMACH(1) / Z"3001800000000000" /
C     DATA RMACH(2) / Z"4FFEFFFFFFFFFFFE" /
C     DATA RMACH(3) / Z"3FD2800000000000" /
C     DATA RMACH(4) / Z"3FD3800000000000" /
C     DATA RMACH(5) / Z"3FFF9A209A84FBCF" /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA RMACH(1) / 00564000000000000000B /
C     DATA RMACH(2) / 37767777777777777776B /
C     DATA RMACH(3) / 16414000000000000000B /
C     DATA RMACH(4) / 16424000000000000000B /
C     DATA RMACH(5) / 17164642023241175720B /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7F7FFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fn COMPILER OPTION
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7FFFFFFF' /
C     DATA RMACH(3) / Z'34800000' /
C     DATA RMACH(4) / Z'35000000' /
C     DATA RMACH(5) / Z'3F9A209B' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fi COMPILER OPTION
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -p8 OR -pd8 COMPILER OPTION
C
C     DATA RMACH(1) / Z'0010000000000000' /
C     DATA RMACH(2) / Z'7FFFFFFFFFFFFFFF' /
C     DATA RMACH(3) / Z'3CC0000000000000' /
C     DATA RMACH(4) / Z'3CD0000000000000' /
C     DATA RMACH(5) / Z'3FF34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE CRAY
C
C     DATA RMACH(1) / 200034000000000000000B /
C     DATA RMACH(2) / 577767777777777777776B /
C     DATA RMACH(3) / 377224000000000000000B /
C     DATA RMACH(4) / 377234000000000000000B /
C     DATA RMACH(5) / 377774642023241175720B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
C     STATIC RMACH(5)
C
C     DATA SMALL /    20K,       0 /
C     DATA LARGE / 77777K, 177777K /
C     DATA RIGHT / 35420K,       0 /
C     DATA DIVER / 36020K,       0 /
C     DATA LOG10 / 40423K,  42023K /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING G_FLOAT
C
C     DATA RMACH(1) / '00000080'X /
C     DATA RMACH(2) / 'FFFF7FFF'X /
C     DATA RMACH(3) / '00003480'X /
C     DATA RMACH(4) / '00003500'X /
C     DATA RMACH(5) / '209B3F9A'X /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING IEEE_FLOAT
C
C     DATA RMACH(1) / '00800000'X /
C     DATA RMACH(2) / '7F7FFFFF'X /
C     DATA RMACH(3) / '33800000'X /
C     DATA RMACH(4) / '34000000'X /
C     DATA RMACH(5) / '3E9A209B'X /
C
C     MACHINE CONSTANTS FOR THE DEC RISC
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1) /       128 /
C     DATA LARGE(1) /    -32769 /
C     DATA RIGHT(1) /     13440 /
C     DATA DIVER(1) /     13568 /
C     DATA LOG10(1) / 547045274 /
C
C     DATA SMALL(1) / Z00000080 /
C     DATA LARGE(1) / ZFFFF7FFF /
C     DATA RIGHT(1) / Z00003480 /
C     DATA DIVER(1) / Z00003500 /
C     DATA LOG10(1) / Z209B3F9A /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C     (ASSUMING REAL*4 IS THE DEFAULT REAL)
C
C     DATA SMALL(1) / '00800000'X /
C     DATA LARGE(1) / '7F7FFFFF'X /
C     DATA RIGHT(1) / '33800000'X /
C     DATA DIVER(1) / '34000000'X /
C     DATA LOG10(1) / '3E9A209B'X /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
C     DATA LARGE(1), LARGE(2) / '37777777, '00000177 /
C     DATA RIGHT(1), RIGHT(2) / '20000000, '00000352 /
C     DATA DIVER(1), DIVER(2) / '20000000, '00000353 /
C     DATA LOG10(1), LOG10(2) / '23210115, '00000377 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA RMACH(1) / O402400000000 /
C     DATA RMACH(2) / O376777777777 /
C     DATA RMACH(3) / O714400000000 /
C     DATA RMACH(4) / O716400000000 /
C     DATA RMACH(5) / O776464202324 /
C
C     MACHINE CONSTANTS FOR THE HP 730
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) / 40000B,       1 /
C     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
C     DATA DIVER(1), DIVER(2) / 40000B,    327B /
C     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) / 40000B,       1 /
C     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
C     DATA DIVER(1), DIVER(2) / 40000B,    327B /
C     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA SMALL(1) / 00004000000B /
C     DATA LARGE(1) / 17677777777B /
C     DATA RIGHT(1) / 06340000000B /
C     DATA DIVER(1) / 06400000000B /
C     DATA LOG10(1) / 07646420233B /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA RMACH(1) / Z00100000 /
C     DATA RMACH(2) / Z7FFFFFFF /
C     DATA RMACH(3) / Z3B100000 /
C     DATA RMACH(4) / Z3C100000 /
C     DATA RMACH(5) / Z41134413 /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C
C     DATA SMALL(1) / 1.18E-38      /
C     DATA LARGE(1) / 3.40E+38      /
C     DATA RIGHT(1) / 0.595E-07     /
C     DATA DIVER(1) / 1.19E-07      /
C     DATA LOG10(1) / 0.30102999566 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE INTEL i860
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR)
C
C     DATA RMACH(1) / "000400000000 /
C     DATA RMACH(2) / "377777777777 /
C     DATA RMACH(3) / "146400000000 /
C     DATA RMACH(4) / "147400000000 /
C     DATA RMACH(5) / "177464202324 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1) /    8388608 /
C     DATA LARGE(1) / 2147483647 /
C     DATA RIGHT(1) /  880803840 /
C     DATA DIVER(1) /  889192448 /
C     DATA LOG10(1) / 1067065499 /
C
C     DATA RMACH(1) / O00040000000 /
C     DATA RMACH(2) / O17777777777 /
C     DATA RMACH(3) / O06440000000 /
C     DATA RMACH(4) / O06500000000 /
C     DATA RMACH(5) / O07746420233 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /   128,     0 /
C     DATA LARGE(1), LARGE(2) / 32767,    -1 /
C     DATA RIGHT(1), RIGHT(2) / 13440,     0 /
C     DATA DIVER(1), DIVER(2) / 13568,     0 /
C     DATA LOG10(1), LOG10(2) / 16282,  8347 /
C
C     DATA SMALL(1), SMALL(2) / O000200, O000000 /
C     DATA LARGE(1), LARGE(2) / O077777, O177777 /
C     DATA RIGHT(1), RIGHT(2) / O032200, O000000 /
C     DATA DIVER(1), DIVER(2) / O032400, O000000 /
C     DATA LOG10(1), LOG10(2) / O037632, O020233 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE SUN
C     USING THE -r8 COMPILER OPTION
C
C     DATA RMACH(1) / Z'0010000000000000' /
C     DATA RMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA RMACH(3) / Z'3CA0000000000000' /
C     DATA RMACH(4) / Z'3CB0000000000000' /
C     DATA RMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES
C
C     DATA RMACH(1) / O000400000000 /
C     DATA RMACH(2) / O377777777777 /
C     DATA RMACH(3) / O146400000000 /
C     DATA RMACH(4) / O147400000000 /
C     DATA RMACH(5) / O177464202324 /
C
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
C
C     DATA SMALL(1), SMALL(2) /     0,    256/
C     DATA LARGE(1), LARGE(2) /    -1,   -129/
C     DATA RIGHT(1), RIGHT(2) /     0,  26880/
C     DATA DIVER(1), DIVER(2) /     0,  27136/
C     DATA LOG10(1), LOG10(2) /  8347,  32538/
C
C***FIRST EXECUTABLE STATEMENT  R1MACH
      IF (I .LT. 1 .OR. I .GT. 5) CALL XERMSG ('SLATEC', 'R1MACH',
     +   'I OUT OF BOUNDS', 1, 2)
C
      R1MACH = RMACH(I)
      RETURN
C
      END
