C$Modified: Thu Jan 11 15:00:50 2007 by puwer $
C
C ----------------------------------------------------------------------
C
      SUBROUTINE COUP4X(SW2,ZMASS,FMASS , GCHF)
      implicit none
C
C This subroutine sets up the coupling constant for the fermion-fermion-
C Higgs vertex in the STANDARD MODEL.  The coupling is DOUBLE COMPLEX and the  
C array of the coupling specifies the chirality of the flowing-IN       
C fermion.  GCHF(1) denotes a left-handed coupling, and GCHF(2) a right-
C handed coupling.                                                      
C                                                                       
C INPUT:                                                                
C       double precision    SW2            : square of sine of the weak angle       
C       double precision    ZMASS          : Z       mass                           
C       double precision    FMASS          : fermion mass                           
C                                                                       
C OUTPUT:                                                               
C       double complex GCHF(2)        : coupling of fermion and Higgs          
C
      DOUBLE COMPLEX GCHF(2)
      DOUBLE PRECISION    SW2,ZMASS,FMASS,ALPHA,FOURPI,EZ,G
C
      ALPHA=1.d0/128.d0
C      ALPHA=1./DOUBLE PRECISION(137.0359895)
      FOURPI=4.D0*3.141592653589793238462643383279502884197d0
      EZ=SQRT(ALPHA*FOURPI)/SQRT(SW2*(1.d0-SW2))
      G=EZ*FMASS/ZMASS/2.d0
C
      GCHF(1) = DCMPLX( -G )
      GCHF(2) = DCMPLX( -G )
C
      RETURN
      END
