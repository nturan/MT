! This is part of BSYpptt package by S. Badger, R. Sattler, V. Yundin
! For details see README file

      subroutine bsyspinit(p1)
      implicit none
      double precision p1(4)
      include 'bsyparams.h'
      double complex bsyA,bsyB,bsyAB
      Sa(1)=bsyA(1,2)
      Sa(2)=bsyA(1,3)
      Sa(3)=bsyA(1,4)
      Sb(1)=bsyB(1,2)
      Sb(2)=bsyB(1,3)
      Sb(3)=bsyB(1,4)
      sab(1)=bsyAB(1,p1,1)
      sab(2)=bsyAB(1,p1,2)
      sab(3)=bsyAB(2,p1,1)
      sab(4)=bsyAB(2,p1,2)

      end
