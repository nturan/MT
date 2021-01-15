C$Modified: Thu Jan 11 11:33:56 2007 by puwer $
c
c ----------------------------------------------------------------------
c
      subroutine fsoxxx(fo,sc,gc,fmass,fwidth , fso)
      implicit none
c
c this subroutine computes an off-shell fermion wavefunction from a     
c flowing-out external fermion and a vector boson.                      
c                                                                       
c input:                                                                
c       double complex fo(6)          : flow-out fermion                   <fo|
c       double complex sc(6)          : input    scalar                      s 
c       double complex gc(2)          : coupling constants                 gchf
c       double precision*8     fmass          : mass  of output fermion f'             
c       double precision*8     fwidth         : width of output fermion f'             
c                                                                       
c output:                                                               
c       double complex fso(6)         : off-shell fermion             <fo,s,f'|
c
      double complex fo(6),sc(6),fso(6),gc(2),sl1,sl2,sr1,sr2,ds
      double precision     pf(0:3),fmass,fwidth,pf2,p0p3,p0m3
c
      fso(5) = fo(5)+sc(2)
      fso(6) = fo(6)+sc(3)
c
      pf(0)=dble( fso(5))
      pf(1)=dble( fso(6))
      pf(2)=dimag(fso(6))
      pf(3)=dimag(fso(5))
      pf2=pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)
c
      ds=-sc(1)/dcmplx(pf2-fmass**2,max(dsign(fmass*fwidth ,pf2),0d0))
      p0p3=pf(0)+pf(3)
      p0m3=pf(0)-pf(3)
      sl1=gc(2)*(p0p3*fo(3)      +fso(6) *fo(4))
      sl2=gc(2)*(p0m3*fo(4)+dconjg(fso(6))*fo(3))
      sr1=gc(1)*(p0m3*fo(1)      -fso(6) *fo(2))
      sr2=gc(1)*(p0p3*fo(2)-dconjg(fso(6))*fo(1))
c
      fso(1) = ( gc(1)*fmass*fo(1) + sl1 )*ds
      fso(2) = ( gc(1)*fmass*fo(2) + sl2 )*ds
      fso(3) = ( gc(2)*fmass*fo(3) + sr1 )*ds
      fso(4) = ( gc(2)*fmass*fo(4) + sr2 )*ds
c
      return          
      end
