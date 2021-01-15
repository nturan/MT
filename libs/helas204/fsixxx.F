C$Modified: Thu Jan 11 11:31:18 2007 by puwer $
c
c ----------------------------------------------------------------------
c
      subroutine fsixxx(fi,sc,gc,fmass,fwidth , fsi)
      implicit none
c
c this subroutine computes an off-shell fermion wavefunction from a     
c flowing-in external fermion and a vector boson.                       
c                                                                       
c input:                                                                
c       double complex fi(6)          : flow-in  fermion                   |fi>
c       double complex sc(3)          : input    scalar                      s 
c       double complex gc(2)          : coupling constants                 gchf
c       double precision    fmass          : mass  of output fermion f'             
c       double precision    fwidth         : width of output fermion f'             
c                                                                       
c output:                                                               
c       double complex fsi(6)         : off-shell fermion             |f',s,fi>
c
      double complex fi(6),sc(3),fsi(6),gc(2),sl1,sl2,sr1,sr2,ds
      double precision     pf(0:3),fmass,fwidth,pf2,p0p3,p0m3
c
      fsi(5) = fi(5)-sc(2)
      fsi(6) = fi(6)-sc(3)
c
      pf(0)=dble( fsi(5))
      pf(1)=dble( fsi(6))
      pf(2)=dimag(fsi(6))
      pf(3)=dimag(fsi(5))
      pf2=pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)
c
      ds=-sc(1)/dcmplx(pf2-fmass**2,max(dsign(fmass*fwidth ,pf2),0d0))
      p0p3=pf(0)+pf(3)
      p0m3=pf(0)-pf(3)
      sl1=gc(1)*(p0p3*fi(1)+dconjg(fsi(6))*fi(2))
      sl2=gc(1)*(p0m3*fi(2)      +fsi(6) *fi(1))
      sr1=gc(2)*(p0m3*fi(3)-dconjg(fsi(6))*fi(4))
      sr2=gc(2)*(p0p3*fi(4)      -fsi(6) *fi(3))
c
      fsi(1) = ( gc(1)*fmass*fi(1) + sr1 )*ds
      fsi(2) = ( gc(1)*fmass*fi(2) + sr2 )*ds
      fsi(3) = ( gc(2)*fmass*fi(3) + sl1 )*ds
      fsi(4) = ( gc(2)*fmass*fi(4) + sl2 )*ds
c
      return          
      end
