C***********************************************************************
      SUBROUTINE SPLINE (N,X,Y,FDP)      
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
      DIMENSION X(N),Y(N),A(N),B(N),C(N),R(N),FDP(N)  
C-----COMPUTE THE COEFFICIENTS AND THE RHS OF THE EQUATIONS. 
C-----THIS ROUTINE USES THE CANTILEVER CONDITION.  THE PARAMETER       
C-----ALAMDA (LAMBDA) IS SET TO 1. BUT THIS CAN BE USER-MODIFIED.      
C-----A,B,C ARE THE THREE DIAGONALS OF THE TRIDIAGONAL SYSTEM;         
C-----R IS THE RIGHT HAND SIDE.  THESE ARE NOW ASSEMBLED.    
      ALAMDA = 1.    
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
