C$Modified: Thu Jan 11 10:55:34 2007 by puwer $
c
c ----------------------------------------------------------------------
c
      subroutine jvsxxx(vc,sc,g,vmass,vwidth , jvs)
      implicit double precision(a-h,o-z)
c
c this subroutine computes an off-shell vector current from the vector- 
c vector-scalar coupling.  the vector propagator is given in feynman    
c gauge for a massless vector and in unitary gauge for a massive vector.
c                                                                       
c input:                                                                
c       double complex vc(6)          : input vector                          v
c       double complex sc(3)          : input scalar                          s
c       double precision    g              : coupling constant                  gvvh
c       double precision    vmass          : mass  of output vector v'              
c       double precision    vwidth         : width of output vector v'              
c                                                                       
c output:                                                               
c       double complex jvs(6)         : vector current             j^mu(v':v,s)
c
      double complex vc(6),sc(3),jvs(6),dg,vk
      double precision    q(0:3),vmass,vwidth,q2,vm2,g
c
      jvs(5) = vc(5)+sc(2)
      jvs(6) = vc(6)+sc(3)
c
      q(0)=dble( jvs(5))
      q(1)=dble( jvs(6))
      q(2)=dimag(jvs(6))
      q(3)=dimag(jvs(5))
      q2=q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2=vmass**2
c
      if (vmass.eq.0.) goto 10
c
      dg=g*sc(1)/dcmplx( q2-vm2 , max(dsign( vmass*vwidth ,q2),0.d0) )
c  for the running width, use below instead of the above dg.
c      dg=g*sc(1)/dcmplx( q2-vm2 , max( vwidth*q2/vmass ,0.) )
c
      vk=(-q(0)*vc(1)+q(1)*vc(2)+q(2)*vc(3)+q(3)*vc(4))/vm2
c
      jvs(1) = dg*(q(0)*vk+vc(1))
      jvs(2) = dg*(q(1)*vk+vc(2))
      jvs(3) = dg*(q(2)*vk+vc(3))
      jvs(4) = dg*(q(3)*vk+vc(4))
c
      return
c
  10  dg=g*sc(1)/q2
c
      jvs(1) = dg*vc(1)
      jvs(2) = dg*vc(2)
      jvs(3) = dg*vc(3)
      jvs(4) = dg*vc(4)
c
      return
      end


