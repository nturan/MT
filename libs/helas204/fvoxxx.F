C$Modified: Thu Jan 11 11:31:57 2007 by puwer $
c
c ----------------------------------------------------------------------
c
      subroutine fvoxxx(fo,vc,g,fmass,fwidth , fvo)
      implicit none
c
c this subroutine computes an off-shell fermion wavefunction from a     
c flowing-out external fermion and a vector boson.                      
c                                                                       
c input:                                                                
c       double complex fo(6)          : flow-out fermion                   <fo|
c       double complex vc(6)          : input    vector                      v 
c       double precision    g(2)           : coupling constants                  gvf
c       double precision    fmass          : mass  of output fermion f'             
c       double precision    fwidth         : width of output fermion f'             
c                                                                       
c output:                                                               
c       double complex fvo(6)         : off-shell fermion             <fo,v,f'|
c
      double complex fo(6),vc(6),fvo(6),sl1,sl2,sr1,sr2,d
      double precision    g(2),pf(0:3),fmass,fwidth,pf2
c
      double precision r_zero, r_one
      parameter( r_zero=0.0d0, r_one=1.0d0 )
      double complex c_imag
**      parameter( c_imag=dcmplx( r_zero, r_one ) )
      parameter( c_imag=( 0.d0, 1.d0 ) )
c
      fvo(5) = fo(5)+vc(5)
      fvo(6) = fo(6)+vc(6)
c
      pf(0)=dble( fvo(5))
      pf(1)=dble( fvo(6))
      pf(2)=dimag(fvo(6))
      pf(3)=dimag(fvo(5))
      pf2=pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)
c
      d=-r_one/dcmplx( pf2-fmass**2,max(sign(fmass*fwidth,pf2),r_zero))
      sl1= (vc(1)+       vc(4))*fo(3)
     &    +(vc(2)+c_imag*vc(3))*fo(4)
      sl2= (vc(2)-c_imag*vc(3))*fo(3)
     &    +(vc(1)-       vc(4))*fo(4)
c
      if ( g(2) .ne. r_zero ) then
         sr1= (vc(1)-       vc(4))*fo(1)
     &       -(vc(2)+c_imag*vc(3))*fo(2)
         sr2=-(vc(2)-c_imag*vc(3))*fo(1)
     &       +(vc(1)+       vc(4))*fo(2)
c
         fvo(1) = ( g(2)*( (pf(0)+pf(3))*sr1        +fvo(6)*sr2)
     &             +g(1)*fmass*sl1)*d
         fvo(2) = ( g(2)*( dconjg(fvo(6))*sr1 +(pf(0)-pf(3))*sr2)
     &             +g(1)*fmass*sl2)*d
         fvo(3) = ( g(1)*( (pf(0)-pf(3))*sl1        -fvo(6)*sl2)
     &             +g(2)*fmass*sr1)*d
         fvo(4) = ( g(1)*(-dconjg(fvo(6))*sl1 +(pf(0)+pf(3))*sl2)
     &             +g(2)*fmass*sr2)*d
c
      else
         fvo(1) = g(1)*fmass*sl1*d
         fvo(2) = g(1)*fmass*sl2*d
         fvo(3) = g(1)*( (pf(0)-pf(3))*sl1        -fvo(6)*sl2)*d
         fvo(4) = g(1)*(-dconjg(fvo(6))*sl1 +(pf(0)+pf(3))*sl2)*d
      end if
c
      return          
      end
