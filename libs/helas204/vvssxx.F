C$Modified: Thu Jan 11 11:08:00 2007 by puwer $
C
      SUBROUTINE VVSSXX(V1,V2,S1,S2,G , VERTEX)
      implicit none
C
C This subroutine computes an amplitude of the vector-vector-scalar-    
C scalar coupling.                                                      
C                                                                       
C INPUT:                                                                
C       double complex V1(6)          : first  vector                        V1
C       double complex V2(6)          : second vector                        V2
C       double complex S1(3)          : first  scalar                        S1
C       double complex S2(3)          : second scalar                        S2
C       double precision    G              : coupling constant                 GVVHH
C                                                                       
C OUTPUT:                                                               
C       double complex VERTEX         : amplitude            Gamma(V1,V2,S1,S2)
C
      DOUBLE COMPLEX V1(6),V2(6),S1(3),S2(3),VERTEX
      DOUBLE PRECISION    G
C
      VERTEX = G*S1(1)*S2(1)
     &        *(V1(1)*V2(1)-V1(2)*V2(2)-V1(3)*V2(3)-V1(4)*V2(4))
C
      RETURN
      END
C
