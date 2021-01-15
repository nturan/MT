C$Modified: Thu Jan 11 11:10:06 2007 by puwer $
C ======================================================================
C
      SUBROUTINE SSSSXX(S1,S2,S3,S4,G , VERTEX)
C
C This subroutine computes an amplitude of the four-scalar coupling.    
C                                                                       
C INPUT:                                                                
C       double complex S1(3)          : first  scalar                        S1
C       double complex S2(3)          : second scalar                        S2
C       double complex S3(3)          : third  scalar                        S3
C       double complex S4(3)          : fourth scalar                        S4
C       double precision    G              : coupling constant                 GHHHH
C                                                                       
C OUTPUT:                                                               
C       double complex VERTEX         : amplitude            Gamma(S1,S2,S3,S4)
C
      implicit none
      DOUBLE COMPLEX S1(3),S2(3),S3(3),S4(3),VERTEX
      DOUBLE PRECISION     G
C
      VERTEX = G*S1(1)*S2(1)*S3(1)*S4(1)
C
      RETURN
      END
C
