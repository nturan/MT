C$Modified: Thu Jan 11 15:07:56 2007 by puwer $
C
C ----------------------------------------------------------------------
C
      SUBROUTINE MOM2CX(ESUM,MASS1,MASS2,COSTH1,PHI1 , P1,P2)
      implicit none
C
C This subroutine sets up two four-momenta in the two particle rest     
C frame.                                                                
C                                                                       
C INPUT:                                                                
C       double precision    ESUM           : energy sum of particle 1 and 2         
C       double precision    MASS1          : mass            of particle 1          
C       double precision    MASS2          : mass            of particle 2          
C       double precision    COSTH1         : cos(theta)      of particle 1          
C       double precision    PHI1           : azimuthal angle of particle 1          
C                                                                       
C OUTPUT:                                                               
C       double precision    P1(0:3)        : four-momentum of particle 1            
C       double precision    P2(0:3)        : four-momentum of particle 2            
C
      DOUBLE PRECISION   P1(0:3),P2(0:3),
     &        ESUM,MASS1,MASS2,COSTH1,PHI1,MD2,ED,PP,SINTH1
      double precision r_half
      parameter( r_half=1.d0/2.d0)
C
      MD2=(MASS1-MASS2)*(MASS1+MASS2)
      ED=MD2/ESUM
      IF (MASS1*MASS2.EQ.0.) THEN
      PP=(ESUM-ABS(ED))*r_half
C
      ELSE
      PP=SQRT((MD2/ESUM)**2-2.0d0*(MASS1**2+MASS2**2)+ESUM**2)*r_half
      ENDIF
      SINTH1=SQRT((1.0d0-COSTH1)*(1.0d0+COSTH1))
C
      P1(0) = MAX((ESUM+ED)*r_half,0.d0)
      P1(1) = PP*SINTH1*COS(PHI1)
      P1(2) = PP*SINTH1*SIN(PHI1)
      P1(3) = PP*COSTH1
C
      P2(0) = MAX((ESUM-ED)*r_half,0.d0)
      P2(1) = -P1(1)
      P2(2) = -P1(2)
      P2(3) = -P1(3)
C
      RETURN
      END
