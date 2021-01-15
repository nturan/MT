C$Modified: Thu Jan 11 11:17:28 2007 by puwer $
C ----------------------------------------------------------------------
C
      SUBROUTINE HSSXXX(S1,S2,G,SMASS,SWIDTH , HSS)
      implicit none
C
C This subroutine computes an off-shell scalar current from the three-  
C scalar coupling.                                                      
C                                                                       
C INPUT:                                                                
C       double complex S1(3)          : first  scalar                        S1
C       double complex S2(3)          : second scalar                        S2
C       double precision    G              : coupling constant                  GHHH
C       double precision    SMASS          : mass  of OUTPUT scalar S'              
C       double precision    SWIDTH         : width of OUTPUT scalar S'              
C                                                                       
C OUTPUT:                                                               
C       double complex HSS(3)         : scalar current              J(S':S1,S2)
C
      DOUBLE COMPLEX S1(3),S2(3),HSS(3),DG
      DOUBLE PRECISION  Q(0:3),G,SMASS,SWIDTH,Q2
C
      HSS(2) = S1(2)+S2(2)
      HSS(3) = S1(3)+S2(3)
C
      Q(0)=dble( HSS(2))
      Q(1)=dble( HSS(3))
      Q(2)=dIMAG(HSS(3))
      Q(3)=dIMAG(HSS(2))
      Q2=Q(0)**2-(Q(1)**2+Q(2)**2+Q(3)**2)
C
      DG=-G/dCMPLX( Q2-SMASS**2, MAX(SIGN(SMASS*SWIDTH ,Q2),0.d0))
C
      HSS(1) = DG*S1(1)*S2(1)
C
      RETURN
      END
C
C ======================================================================
