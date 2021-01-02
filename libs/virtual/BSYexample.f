! This is part of BSYpptt package by S. Badger, R. Sattler, V. Yundin
! For details see README file

      program main
      implicit none
      double precision p1(4),p2(4),p3(4),p4(4),Ecms
      double precision m,pi,theta,phi,mur2,gs
      double complex ggamp(-2:0),qqamp(-2:0)
      double complex ggtree(1:2,1:2),qqtree

      m=1.75D0
      Ecms=dsqrt(m**2+1d0)
      pi=4*datan(1D0)
      theta=pi/3D0
      phi=pi/4D0

      p3(4)=Ecms
      p3(1)=-1D0
      p3(2)=0D0
      p3(3)=0D0

      p4(4)=Ecms
      p4(1)=1D0
      p4(2)=0D0
      p4(3)=0D0

      p1(4)=Ecms
      p1(1)=-Ecms*dsin(theta)
      p1(2)=-Ecms*dcos(theta)*dcos(phi)
      p1(3)=-Ecms*dcos(theta)*dsin(phi)

      p2(4)=Ecms
      p2(1)=Ecms*dsin(theta)
      p2(2)=Ecms*dcos(theta)*dcos(phi)
      p2(3)=Ecms*dcos(theta)*dsin(phi)

      write(*,*) "p1",p1
      write(*,*) "p2",p2
      write(*,*) "p3",p3
      write(*,*) "p4",p4

      call bsyppttinit(m,5,1)

      mur2=1.75**2;
      write (*,*) "Mu_R^2=", mur2
      call bsyggttsq(p1,p2,p3,p4,mur2,ggamp,ggtree)
      write(*,*) "gg_tree",4D0*( 
     #  +(ggtree(1,1)+ggtree(2,2))*16D0/3D0 
     #  -(ggtree(1,2)+ggtree(1,2))*2D0/3D0
     # )*64D0*4D0
      write(*,*) "gg(-2)",ggamp(-2)*64D0*4D0
      write(*,*) "gg(-1)",ggamp(-1)*64D0*4D0
      write(*,*) "gg( 0)",ggamp( 0)*64D0*4D0
      call bsyqqttsq(p1,p2,p3,p4,mur2,qqamp,qqtree)
      write(*,*) "qq_tree",qqtree*9D0*4D0
      write(*,*) "qq(-2)",qqamp(-2)*9D0*4D0
      write(*,*) "qq(-1)",qqamp(-1)*9D0*4D0
      write(*,*) "qq( 0)",qqamp( 0)*9D0*4D0

      mur2=4*1.75**2;
      write (*,*) "Mu_R^2=", mur2
      call bsyggttsq(p1,p2,p3,p4,mur2,ggamp,ggtree)
      write(*,*) "gg_tree",4D0*( 
     #  +(ggtree(1,1)+ggtree(2,2))*16D0/3D0 
     #  -(ggtree(1,2)+ggtree(1,2))*2D0/3D0
     # )*64D0*4D0
      write(*,*) "gg(-2)",ggamp(-2)*64D0*4D0
      write(*,*) "gg(-1)",ggamp(-1)*64D0*4D0
      write(*,*) "gg( 0)",ggamp( 0)*64D0*4D0
      call bsyqqttsq(p1,p2,p3,p4,mur2,qqamp,qqtree)
      write(*,*) "qq_tree",qqtree*9D0*4D0
      write(*,*) "qq(-2)",qqamp(-2)*9D0*4D0
      write(*,*) "qq(-1)",qqamp(-1)*9D0*4D0
      write(*,*) "qq( 0)",qqamp( 0)*9D0*4D0

      write(*,*)  "-----------------------------------------------"
      write(*,*)  "-----------------------------------------------"

      write(*,*) "p1",p1
      write(*,*) "p2",p2
      write(*,*) "p3",p3
      write(*,*) "p4",p4


      end
