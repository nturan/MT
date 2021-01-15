C$Modified: Thu Jan 11 11:33:01 2007 by puwer $
C
C ======================================================================
C
      SUBROUTINE IOSXXX(FI,FO,SC,GC , VERTEX)
      implicit none
C
C This subroutine computes an amplitude of the fermion-fermion-scalar   
C coupling.                                                             
C                                                                       
C INPUT:                                                                
C       double complex FI(6)          : flow-in  fermion                   |FI>
C       double complex FO(6)          : flow-out fermion                   <FO|
C       double complex SC(3)          : input    scalar                      S 
C       double complex GC(2)          : coupling constants                 GCHF
C                                                                       
C OUTPUT:                                                               
C       double complex VERTEX         : amplitude                     <FO|S|FI>
C
      DOUBLE COMPLEX FI(6),FO(6),SC(3),GC(2),VERTEX
C
      VERTEX = SC(1)*( GC(1)*(FI(1)*FO(1)+FI(2)*FO(2))
     &                +GC(2)*(FI(3)*FO(3)+FI(4)*FO(4)) )
C
      RETURN          
      END
