! This is part of BSYpptt package by S. Badger, R. Sattler, V. Yundin
! For details see README file

      subroutine bsyspadd(p,n)
      implicit none
      integer n
      double precision p(4)

      include 'bsyparams.h'

      double precision s0
      double complex kp,km,ep,em

      s0=1D0
      s0=dsign(s0,p(4))
      kp=dsqrt(s0*(p(4)+p(2)))
      km=dsqrt(s0*(p(4)-p(2)))
      if (dabs(p(4))-dabs(p(2)) .eq. 0D0) then
        ep=(0D0,1D0)
      else
        ep=s0*dcmplx(p(1),p(3))
        ep=ep/(km*kp)
      endif
      em=dconjg(ep)
      if (s0 .lt. 0d0) then
        kp=dcmplx(-dimag(kp),dreal(kp))
        km=dcmplx(-dimag(km),dreal(km))
      endif
!       write (*,*) kp,km,ep,em
      sp(n,1)=kp
      sp(n,2)=km
      sp(n,3)=ep
      sp(n,4)=em
      end

      double complex function bsyA(i,j)
      implicit none
      integer i,j

      include 'bsyparams.h'

!       ncomplex val=m1.km*m2.kp*m1.ep-m2.km*m1.kp*m2.ep;
      bsyA=sp(i,2)*sp(j,1)*sp(i,3)-sp(j,2)*sp(i,1)*sp(j,3)
      return
      end


      double complex function bsyB(i,j)
      implicit none
      integer i,j

      include 'bsyparams.h'

!       ncomplex val=-m1.km*m2.kp*m1.em+m2.km*m1.kp*m2.em;
      bsyB=-sp(i,2)*sp(j,1)*sp(i,4)+sp(j,2)*sp(i,1)*sp(j,4)
      return
      end

      double complex function bsyAB(i,p,j)
      implicit none
      integer i,j
      double precision p(4)

      include 'bsyparams.h'

!       ncomplex val=-m1.km*m2.kp*m1.em+m2.km*m1.kp*m2.em;
      bsyAB=sp(i,1)*sp(j,1)*(p(4)-p(2))
     &+sp(i,2)*sp(j,2)*sp(i,3)*sp(j,4)*(p(4)+p(2))
     &-sp(i,2)*sp(j,1)*sp(i,3)*dcmplx(p(1),-p(3))
     &-sp(i,1)*sp(j,2)*sp(j,4)*dcmplx(p(1),+p(3))

      return
      end
