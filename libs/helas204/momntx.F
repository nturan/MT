C$Modified: Thu Jan 11 15:37:44 2007 by puwer $
C **********************************************************************
C
      SUBROUTINE MOMNTX(ENERGY,MASS,COSTH,PHI , P)
C
C This subroutine sets up a four-momentum from the four inputs.         
C                                                                       
C INPUT:                                                                
C       double precision    ENERGY         : energy                                 
C       double precision    MASS           : mass                                   
C       double precision    COSTH          : cos(theta)                             
C       double precision    PHI            : azimuthal angle                        
C                                                                       
C OUTPUT:                                                               
C       double precision    P(0:3)         : four-momentum                          
C
      implicit none
      DOUBLE PRECISION    P(0:3),ENERGY,MASS,COSTH,PHI,PP,SINTH
C
      P(0) = ENERGY
      IF (ENERGY.EQ.MASS) THEN
         P(1) = 0.d0
         P(2) = 0.d0
         P(3) = 0.d0
      ELSE
         PP=SQRT((ENERGY-MASS)*(ENERGY+MASS))
         SINTH=SQRT((1.d0-COSTH)*(1.d0+COSTH))
         P(3) = PP*COSTH
         IF (PHI.EQ.0.) THEN
            P(1) = PP*SINTH
            P(2) = 0.d0
         ELSE
            P(1) = PP*SINTH*COS(PHI)
            P(2) = PP*SINTH*SIN(PHI)
         ENDIF
      ENDIF
      RETURN
      END
C
