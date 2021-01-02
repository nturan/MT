! This is part of BSYpptt package by S. Badger, R. Sattler, V. Yundin
! For details see README file

      subroutine bsyppttinit(topmass,lightflavours,schemechoice)
      implicit none
      include 'bsyparams.h'
      double precision topmass
      integer lightflavours,schemechoice

!     0 = FDH, 1 = HV
      scheme=schemechoice
      nf=lightflavours
      m=topmass
      m2=m*m
      pi=4*datan(1D0)

      call qlinit()

      end

      subroutine bsyggttsq(pp1,pp2,pp3,pp4,mur2,sqamp,sqtree)
      implicit none
      include 'bsyparams.h'
      double precision pp1(4),pp2(4),pp3(4),pp4(4),mur2
      double precision p1(4),p2(4),p3(4),p4(4)
      double precision xbeta, sig23, sig14
      double complex sq0aa,sq0ab,sq0ba,sq0bb
      double complex PggR23,PggL23,PggL32,PggSL
      double complex sqA0,sqA0a,sqA0b,sqA0sl,sqAL,sqAR,sqAf,sqAH,sqAsl
      double complex ren0,ren1,sch
      double complex sqamp(-2:0)
      double complex sqtree(1:2,1:2)
      include 'bsymom.h'

      call bsyggttint()
      call bsyggtt()

      include 'ampggsums.f'

      sqtree(1,1)=sq0aa
      sqtree(1,2)=sq0ab
      sqtree(2,1)=sq0ba
      sqtree(2,2)=sq0bb
      
      xbeta=dsqrt(xbeta2)
      if (s14 .gt. 0D0) then
        sig23=1D0
      else
        sig23=0D0
      endif
      if (s14-2D0*m2 .gt. 0D0) then
        sig14=1D0
      else
        sig14=0D0
      endif
! note that this means for 1.2>0 and 1.3>0 the branch cuts for the poles
! are not correct. Since the ttbar will always be in the final state 
! it is not a problem.

      PggR23=-0.5D0+(s14-2D0*m2)/(s14*xbeta)*
     # dcmplx( dlog(abs((1D0-xbeta)/(1D0+xbeta))), pi*sig14 )

      PggL23=+0.5D0-dlog(m2*musq/(dreal(sab(1))*dreal(sab(1))))
     # -dcmplx( dlog(musq/dabs(s14)), pi*sig23 )

      PggL32=+0.5D0-dlog(m2*musq/(dreal(sab(4))*dreal(sab(4))))
     # -dcmplx( dlog(musq/dabs(s14)), pi*sig23 )

      PggSL=(
     # -dreal(sab(1))*dlog(m2*musq/(dreal(sab(1))*dreal(sab(1))))
     # -dreal(sab(4))*dlog(m2*musq/(dreal(sab(4))*dreal(sab(4))))
     # )/s14 - dcmplx( dlog(musq/abs(s14)) , pi*sig23 )
     # -(s14-2D0*m2)/(s14*xbeta)*
     # dcmplx( dlog(dabs((1D0-xbeta)/(1D0+xbeta))) , pi*sig14 )

      sqamp( 0)=2D0*(sqAL+sqAR+sqAsl+sqAH+nf*sqAf)
      sqamp(-1)=2D0*(PggR23*(sqA0a + sqA0b)
     # +9D0*(PggL23*sqA0a + PggL32*sqA0b + PggSL*sqA0sl))/3D0
      sqamp(-2)=-12D0*sqA0

      ren1 = 11D0 - 2D0/3D0*nf + 4D0
      ren0 = -1D0 + 4D0/3D0*( 3D0*dlog(musq/m2) + 5D0 )

      sch = 0d0
      if ( scheme.eq.1 ) then
        sch = 1D0
      endif

      sqamp(0) = sqamp(0)-2D0*(ren0+sch)*sqA0
      sqamp(-1) = sqamp(-1)-2D0*ren1*sqA0

      sqamp(0) = sqamp(0)/256D0 
      sqamp(-1) = sqamp(-1)/256D0 
      sqamp(-2) = sqamp(-2)/256D0 
      
      sqtree(1,1)=sqtree(1,1)/256D0
      sqtree(1,2)=sqtree(1,2)/256D0
      sqtree(2,1)=sqtree(2,1)/256D0
      sqtree(2,2)=sqtree(2,2)/256D0
 
      end

      subroutine bsyqqttsq(pp1,pp2,pp3,pp4,mur2,sqamp,sqtree)
      implicit none
      include 'bsyparams.h'
      double precision pp1(4),pp2(4),pp3(4),pp4(4),mur2
      double precision p1(4),p2(4),p3(4),p4(4)
      double precision xbeta, sig23, sig14
      double complex Pqqlc,Pqqlc32,Pqqslc
      double complex sqA0,sqAlc,sqAsl,sqAf,sqAH
      double complex ren0,ren1,sch
      double complex sqamp(-2:0)
      double complex sqtree
      include 'bsymom.h'

      call bsyqqttint()
      call bsyqqtt()

      include 'ampqqsums.f'
      sqtree=sqA0

      xbeta=dsqrt(xbeta2)
      if (s14 .gt. 0D0) then
        sig23=1D0
      else
        sig23=0D0
      endif
      if (s14-2D0*m2 .gt. 0D0) then
        sig14=1D0
      else
        sig14=0D0
      endif

      Pqqlc32=8D0/3D0-dlog(m2*musq/(dreal(sab(1))*dreal(sab(1))))
      Pqqlc=-8D0/3D0+dlog(m2*musq/(dreal(sab(4))*dreal(sab(4))))

      Pqqslc =-1D0-dcmplx(dlog(musq/abs(s14)),pi*sig23)
     # -(s14-2D0*m2)/(s14*xbeta)*
     # dcmplx(dlog(dabs((1D0-xbeta)/(1D0+xbeta))),pi*sig14)

      sqamp(0)=2D0*(sqAlc+sqAsl+(sqAH+nf*sqAf))
      sqamp(-1)=-2D0*sqA0*(
     # 9D0*Pqqlc-2D0*(Pqqlc+Pqqlc32-0.5D0*Pqqslc)+(nf+1)*2D0)/3D0
      sqamp(-2)=-16D0/3D0*sqA0

      ren1 = 11D0 - 2D0/3D0*(nf+1) + 4D0
      ren0 = -1D0 + 4D0/3D0*( 3D0*dlog(musq/m2) +5D0 )
     # - 2D0/3D0*dlog(musq/m2)

      sch = 0d0
      if( scheme.eq.1 ) then
        sch = 4D0/3D0
      endif

      sqamp(0) = sqamp(0)-2D0*(ren0+sch)*sqA0
      sqamp(-1) = sqamp(-1)-2D0*ren1*sqA0

      sqamp(0) = sqamp(0)/36D0 
      sqamp(-1) = sqamp(-1)/36D0 
      sqamp(-2) = sqamp(-2)/36D0 
      
      sqtree=sqtree/36D0

      end
