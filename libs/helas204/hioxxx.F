C$Modified: Thu Jan 11 11:32:20 2007 by puwer $
      subroutine hioxxx(fi,fo,gc,smass,swidth , hio)                         
      implicit none
c                                                                            
c this subroutine computes an off-shell scalar current from an external       
c fermion pair.                                                               
c                                                                              
c input:                                                                      
c       double complex fi(6)          : flow-in  fermion                   |fi>     
c       double complex fo(6)          : flow-out fermion                   <fo|       
c       double complex gc(2)          : coupling constants                 gchf       
c       double precision    smass          : mass  of output scalar s                      
c       double precision    swidth         : width of output scalar s                      
c                                                                              
c output:                                                                      
c       double complex hio(3)         : scalar current             j(<fi|s|fo>)       
c                                                                              
      double complex fi(6),fo(6),hio(3),gc(2),dn                                   
      double precision q(0:3),smass,swidth,q2                                           
c                                                                              
      hio(2) = fo(5)-fi(5)                                                     
      hio(3) = fo(6)-fi(6)                                                     
c                                                                              
      q(0)=dble( hio(2))                                                       
      q(1)=dble( hio(3))                                                       
      q(2)=dimag(hio(3))                                                       
      q(3)=dimag(hio(2))                                                       
      q2=q(0)**2-(q(1)**2+q(2)**2+q(3)**2)                                     
c                                                                              
      dn=-dcmplx(q2-smass**2,dmax1(dsign(smass*swidth,q2),0.d0))               
c                                                                              
      hio(1) = ( gc(1)*(fo(1)*fi(1)+fo(2)*fi(2))                               
     &          +gc(2)*(fo(3)*fi(3)+fo(4)*fi(4)) )/dn                          
c                                                                              
      return                                                                   
      end                                                                      

