C$Modified: Thu Jan 11 11:00:59 2007 by puwer $
c
c ======================================================================
c
      subroutine vvsxxx(v1,v2,sc,g , vertex)
      implicit none
c
c this subroutine computes an amplitude of the vector-vector-scalar     
c coupling.                                                             
c                                                                       
c input:                                                                
c       double complex v1(6)          : first  vector                        v1
c       double complex v2(6)          : second vector                        v2
c       double complex sc(3)          : input  scalar                        s 
c       double precision    g              : coupling constant                  gvvh
c                                                                       
c output:                                                               
c       double complex vertex         : amplitude                gamma(v1,v2,s)
c
      double complex v1(6),v2(6),sc(3),vertex
      double precision    g
c
      vertex = g*sc(1)*(v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)-v1(4)*v2(4))
c
      return
      end
