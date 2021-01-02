! This is part of BSYpptt package by S. Badger, R. Sattler, V. Yundin
! For details see README file

      subroutine bsyqqttint()
!     0 -> t(1) qbar(2) q(3) tbar(4)
!     initialize basis scalar one-loop integrals

      implicit none
      include 'bsyparams.h'

      double precision zip
      parameter (zip=0D0)
      double complex qlI4,qlI3,qlI2

      double complex aint0,aint1,aint2

      aint0=(2*sab(3)*sab(2))/((sab(1)**2*sa(1)*sb(1)))
      aint1=(sab(1)+2*m2)/((sab(1)))
      aint2=(aint1)/((sab(1)))

      scin(1)=qlI2(m2,zip,m2,musq,0)
      scin(2)=qlI2(s12,zip,m2,musq,0)
      scin(3)=-scin(1)+scin(2)
      scin(4)=qlI2(s14,zip,zip,m2,0)
      scin(5)=qlI3(zip,zip,s14,zip,zip,zip,musq,0)
      scin(6)=qlI3(s12,zip,m2,m2,zip,zip,musq,0)
      scin(7)=qlI3(s14,m2,m2,zip,zip,m2,musq,0)
      scin(8)=qlI4(zip,m2,m2,zip,s12,s14,zip,zip,m2,zip,musq,0)
      scin(9)=qlI2(s14,zip,zip,musq,0)
      scin(10)=-2+scin(1)
      scin(11)=2-scin(1)+qlI2(s14,m2,m2,musq,0)
      scin(12)=qlI3(s14,m2,m2,m2,m2,zip,musq,0)
      scin(13)=qlI2(s13,zip,m2,musq,0)
      scin(14)=-scin(1)+scin(13)
      scin(15)=qlI3(s13,zip,m2,m2,zip,zip,musq,0)
      scin(16)=qlI4(zip,m2,m2,zip,s13,s14,zip,zip,m2,zip,musq,0)

      end
