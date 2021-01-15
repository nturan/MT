C$Modified: Thu Jan 11 10:54:24 2007 by puwer $
c
c ----------------------------------------------------------------------
c
      subroutine jtioxx(fi,fo,g , jio)
      implicit none
c
c this subroutine computes an off-shell vector current from an external 
c fermion pair.  the vector boson propagator is not included in this
c routine.
c                                                                       
c input:                                                                
c       double complex fi(6)          : flow-in  fermion                   |fi>
c       double complex fo(6)          : flow-out fermion                   <fo|
c       double precision    g(2)           : coupling constants                  gvf
c                                                                       
c output:                                                               
c       double complex jio(6)         : vector current          j^mu(<fo|v|fi>)
c
      double complex fi(6),fo(6),jio(6)
      double precision    g(2)
c
      double precision r_zero, r_one
      parameter( r_zero=0.0d0, r_one=1.0d0 )
      double complex c_imag
**      parameter( c_imag=dcmplx( r_zero, r_one ) )
      parameter( c_imag=( 0.d0, 1.d0 ) )
c
      jio(5) = fo(5)-fi(5)
      jio(6) = fo(6)-fi(6)
c
      if ( g(2) .ne. r_zero ) then
         jio(1) = ( g(1)*( fo(3)*fi(1)+fo(4)*fi(2))
     &             +g(2)*( fo(1)*fi(3)+fo(2)*fi(4)) )
         jio(2) = (-g(1)*( fo(3)*fi(2)+fo(4)*fi(1))
     &             +g(2)*( fo(1)*fi(4)+fo(2)*fi(3)) )
         jio(3) = ( g(1)*( fo(3)*fi(2)-fo(4)*fi(1))
     &             +g(2)*(-fo(1)*fi(4)+fo(2)*fi(3)) )*c_imag
         jio(4) = ( g(1)*(-fo(3)*fi(1)+fo(4)*fi(2))
     &             +g(2)*( fo(1)*fi(3)-fo(2)*fi(4)) )
c
      else
         jio(1) =  ( fo(3)*fi(1)+fo(4)*fi(2))*g(1)
         jio(2) = -( fo(3)*fi(2)+fo(4)*fi(1))*g(1)
         jio(3) =  ( fo(3)*fi(2)-fo(4)*fi(1))*dcmplx(r_zero,g(1))
         jio(4) =  (-fo(3)*fi(1)+fo(4)*fi(2))*g(1)
      end if
c
      return
      end
