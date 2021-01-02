! This is part of BSYpptt package by S. Badger, R. Sattler, V. Yundin
! For details see README file

      integer maxsprs,maxsp2,maxsp3
      parameter(maxsprs=4)
      parameter(maxsp2=3)
      parameter(maxsp3=4)

      double precision m,m2,musq,xbeta2,s12,s13,s14,pi
      integer nf,scheme
      common/bsyproc/m,m2,musq,xbeta2,s12,s13,s14,pi,nf,scheme

      double complex sp(maxsprs,4)
      double complex sa(maxsp2)
      double complex sb(maxsp2)
      double complex sab(maxsp3)
      common/bsyspin/sp,sa,sb,sab

      integer hels
      parameter(hels=16)
      double complex ggA0(hels),ggB0(hels)
      double complex ggAf(hels),ggBf(hels)
      double complex ggAh(hels),ggBh(hels)
      double complex ggAl(hels),ggBl(hels)
      double complex ggAr(hels),ggBr(hels)
      double complex ggAsl(hels),ggBsl(hels)
      common/bsygghel/ggA0,ggB0,ggAf,ggBf,ggAh,ggBh,
     #                ggAl,ggBl,ggAr,ggBr,ggAsl,ggBsl

      double complex qqA0(hels),qqB0(hels)
      double complex qqAf(hels),qqBf(hels)
      double complex qqAh(hels),qqBh(hels)
      double complex qqAsl(hels),qqBsl(hels)
      double complex qqAlc(hels),qqBlc(hels)
      common/bsyqqhel/qqA0,qqB0,qqAf,qqBf,qqAh,qqBh,
     #                qqAsl,qqBsl,qqAlc,qqBlc

      integer intbasis
      parameter (intbasis=21)
      double complex scin(intbasis)
      common/bsyint/scin
