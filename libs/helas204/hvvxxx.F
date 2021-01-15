C$Modified: Thu Jan 11 11:32:33 2007 by puwer $
c
c ----------------------------------------------------------------------
c
      subroutine hvvxxx(v1,v2,g,smass,swidth , hvv)
      implicit none
c
c this subroutine computes an off-shell scalar current from the vector- 
c vector-scalar coupling.                                               
c                                                                       
c input:                                                                
c       double complex v1(6)          : first  vector                        v1
c       double complex v2(6)          : second vector                        v2
c       double precision    g              : coupling constant                  gvvh
c       double precision    smass          : mass  of output scalar s               
c       double precision    swidth         : width of output scalar s               
c                                                                       
c output:                                                               
c       double complex hvv(3)         : off-shell scalar current     j(s:v1,v2)
c
      double complex v1(6),v2(6),hvv(3),dg
      double precision    q(0:3),g,smass,swidth,q2
c
      double precision r_zero
      parameter( r_zero=0.0d0 )
c
      hvv(2) = v1(5)+v2(5)
      hvv(3) = v1(6)+v2(6)
c
      q(0)=dble( hvv(2))
      q(1)=dble( hvv(3))
      q(2)=dimag(hvv(3))
      q(3)=dimag(hvv(2))
      q2=q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
c
      dg=-g/dcmplx( q2-smass**2 , max(sign( smass*swidth ,q2),r_zero) )
c
      hvv(1) = dg*(v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)-v1(4)*v2(4))
c
      return
      end
