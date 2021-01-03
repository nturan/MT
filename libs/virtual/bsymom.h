! This is part of BSYpptt package by S. Badger, R. Sattler, V. Yundin
! For details see README file

      double precision x
      double precision q1(4),q4(4)

!     take gg   -> ttb to 0 -> tggtb
!     take qqb ->  ttb to 0 ->  tqbqtb

!     modified on 16.05.20 by turan
!     change the energy axis of momenta to 4

      p1(1) = pp3(2)
      p1(2) = pp3(3)
      p1(3) = pp3(4)
      p1(4) = pp3(1)

      p2(1) = -pp1(2)
      p2(2) = -pp1(3)
      p2(3) = -pp1(4)
      p2(4) = -pp1(1)
      
      p3(1) = -pp2(2)
      p3(2) = -pp2(3)
      p3(3) = -pp2(4)
      p3(4) = -pp2(1)
 
      p4(1) = pp4(2)
      p4(2) = pp4(3)
      p4(3) = pp4(4)
      p4(4) = pp4(1)
      
      x=(p1(4)**2-(p1(1)**2+p1(2)**2+p1(3)**2))
     &/(2D0*(p1(4)*p2(4)-(p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3))))
      q1(1)=p1(1)-x*p2(1)
      q1(2)=p1(2)-x*p2(2)
      q1(3)=p1(3)-x*p2(3)
      q1(4)=p1(4)-x*p2(4)

      x=(p4(4)**2-(p4(1)**2+p4(2)**2+p4(3)**2))
     &/(2D0*(p4(4)*p2(4)-(p4(1)*p2(1)+p4(2)*p2(2)+p4(3)*p2(3))))
      q4(1)=p4(1)-x*p2(1)
      q4(2)=p4(2)-x*p2(2)
      q4(3)=p4(3)-x*p2(3)
      q4(4)=p4(4)-x*p2(4)

      call bsyspadd(p2,1)
      call bsyspadd(p3,2)
      call bsyspadd(q1,3)
      call bsyspadd(q4,4)
      call bsyspinit(p1)
     
      musq=mur2
      s12=dreal(sab(1)+m2)
      s13=dreal(sab(4)+m2)
      s14=dreal(-sa(1)*sb(1))
      xbeta2=1D0-4D0*m2/s14
