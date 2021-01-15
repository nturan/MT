C$Modified: Thu Jan 11 15:27:32 2007 by puwer $
c
c ----------------------------------------------------------------------
c
      subroutine coup3x(sw2,zmass,hmass , 
     &                  gwwh,gzzh,ghhh,gwwhh,gzzhh,ghhhh)
      implicit none
c
c this subroutine sets up the coupling constants of the gauge bosons and
c higgs boson in the standard model.                                    
c                                                                       
c input:                                                                
c       double precision    sw2            : square of sine of the weak angle       
c       double precision    zmass          : mass of z                              
c       double precision    hmass          : mass of higgs                          
c                                                                       
c output:                                                               
c       double precision    gwwh           : dimensionful  coupling of w-,w+,h      
c       double precision    gzzh           : dimensionful  coupling of z, z, h      
c       double precision    ghhh           : dimensionful  coupling of h, h, h      
c       double precision    gwwhh          : dimensionful  coupling of w-,w+,h, h   
c       double precision    gzzhh          : dimensionful  coupling of z, z, h, h   
c       double precision    ghhhh          : dimensionless coupling of h, h, h, h   
c
      double precision    sw2,zmass,hmass,gwwh,gzzh,ghhh,gwwhh,gzzhh,
     &        ghhhh,alpha,fourpi,ee2,sc2,v
c
      double precision r_half, r_one, r_two, r_three, r_four, r_ote
      double precision r_pi, r_ialph
      parameter( r_half=1.d0/2.d0, r_one=1.0d0, r_two=2.0d0, 
     &     r_three=3.0d0 )
      parameter( r_four=4.0d0, r_ote=128.0d0 )
      parameter( r_pi=3.141592653589793238462643383279502884197d0,
     &     r_ialph=137.0359895d0 )
c
      alpha = r_one / r_ote
c      alpha = r_one / r_ialph
      fourpi = r_four * r_pi
      ee2=alpha*fourpi
      sc2=sw2*( r_one - sw2 )
      v = r_two * zmass*sqrt(sc2)/sqrt(ee2)
c
      gwwh  =   ee2/sw2*r_half*v
      gzzh  =   ee2/sc2*r_half*v
      ghhh  =  -hmass**2/v*r_three
      gwwhh =   ee2/sw2*r_half
      gzzhh =   ee2/sc2*r_half
      ghhhh = -(hmass/v)**2*r_three
c
      return
      end
