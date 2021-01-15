C$Modified: Thu Jan 11 15:38:26 2007 by puwer $
C
C ----------------------------------------------------------------------
C
      SUBROUTINE SXXXXX(P,NSS , SC)
      implicit none
C
C This subroutine computes a double complex SCALAR wavefunction.               
C                                                                       
C INPUT:                                                                
C       double precision    P(0:3)         : four-momentum of scalar boson          
C       integer NSS  = -1 or 1 : +1 for final, -1 for initial           
C                                                                       
C OUTPUT:                                                               
C       double complex SC(3)          : scalar wavefunction                   S
C
      DOUBLE COMPLEX SC(3)
      DOUBLE PRECISION    P(0:3)
      INTEGER NSS
C
      SC(1) = DCMPLX( 1.d0 )
      SC(2) = DCMPLX(P(0),P(3))*NSS
      SC(3) = DCMPLX(P(1),P(2))*NSS
C
      RETURN
      END
