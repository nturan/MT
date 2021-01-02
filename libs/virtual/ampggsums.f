! This is part of BSYpptt package by S. Badger, R. Sattler, V. Yundin
! For details see README file

      sq0aa=(+ggA0(10)*dconjg(ggA0(10))+ggA0(13)*dconjg(ggA0(13))+ggA0(1
     #4)*dconjg(ggA0(14))+ggA0(7)*dconjg(ggA0(7))+ggA0(4)*dconjg(ggA0(4)
     #)+ggA0(3)*dconjg(ggA0(3))+ggA0(11)*dconjg(ggA0(11))+ggA0(12)*dconj
     #g(ggA0(12))+ggA0(15)*dconjg(ggA0(15))+ggA0(6)*dconjg(ggA0(6))+ggA0
     #(5)*dconjg(ggA0(5))+ggA0(2)*dconjg(ggA0(2)))

      sq0ab=(+ggA0(10)*dconjg(ggB0(10))+ggA0(13)*dconjg(ggB0(16))+ggA0(1
     #4)*dconjg(ggB0(7))+ggA0(7)*dconjg(ggB0(14))+ggA0(4)*dconjg(ggB0(4)
     #)+ggA0(3)*dconjg(ggB0(3))+ggA0(11)*dconjg(ggB0(11))+ggA0(12)*dconj
     #g(ggB0(12))+ggA0(15)*dconjg(ggB0(6))+ggA0(6)*dconjg(ggB0(15))+ggA0
     #(5)*dconjg(ggB0(8))+ggA0(2)*dconjg(ggB0(2)))

      sq0ba=(+ggB0(10)*dconjg(ggA0(10))+ggB0(16)*dconjg(ggA0(13))+ggB0(1
     #4)*dconjg(ggA0(7))+ggB0(7)*dconjg(ggA0(14))+ggB0(4)*dconjg(ggA0(4)
     #)+ggB0(3)*dconjg(ggA0(3))+ggB0(11)*dconjg(ggA0(11))+ggB0(12)*dconj
     #g(ggA0(12))+ggB0(15)*dconjg(ggA0(6))+ggB0(6)*dconjg(ggA0(15))+ggB0
     #(8)*dconjg(ggA0(5))+ggB0(2)*dconjg(ggA0(2)))

      sq0bb=(+ggB0(10)*dconjg(ggB0(10))+ggB0(14)*dconjg(ggB0(14))+ggB0(1
     #6)*dconjg(ggB0(16))+ggB0(7)*dconjg(ggB0(7))+ggB0(4)*dconjg(ggB0(4)
     #)+ggB0(3)*dconjg(ggB0(3))+ggB0(11)*dconjg(ggB0(11))+ggB0(12)*dconj
     #g(ggB0(12))+ggB0(15)*dconjg(ggB0(15))+ggB0(8)*dconjg(ggB0(8))+ggB0
     #(6)*dconjg(ggB0(6))+ggB0(2)*dconjg(ggB0(2)))

      sqA0a=4*((-2D0/3D0)*(sq0ab)+(16D0/3D0)*(sq0aa))

      sqA0b=4*((-2D0/3D0)*(sq0ba)+(16D0/3D0)*(sq0bb))

      sqA0sl=4*((+2D0/3D0)*(sq0aa+sq0ab+sq0ba+sq0bb))

      sqA0=4*((-2D0/3D0)*(sq0ab+sq0ba)+(16D0/3D0)*(sq0aa+sq0bb))

      sqAL=4*((-2D0)*(+ggAl(10)*dconjg(ggB0(10))+ggBl(10)*dconjg(ggA0(10
     #))+ggAl(13)*dconjg(ggB0(16))+ggBl(16)*dconjg(ggA0(13))+ggAl(14)*dc
     #onjg(ggB0(7))+ggBl(14)*dconjg(ggA0(7))+ggBl(7)*dconjg(ggA0(14))+gg
     #Al(7)*dconjg(ggB0(14))+ggAl(4)*dconjg(ggB0(4))+ggBl(4)*dconjg(ggA0
     #(4))+ggAl(3)*dconjg(ggB0(3))+ggBl(3)*dconjg(ggA0(3))+ggAl(11)*dcon
     #jg(ggB0(11))+ggBl(11)*dconjg(ggA0(11))+ggAl(12)*dconjg(ggB0(12))+g
     #gBl(12)*dconjg(ggA0(12))+ggAl(15)*dconjg(ggB0(6))+ggBl(15)*dconjg(
     #ggA0(6))+ggBl(6)*dconjg(ggA0(15))+ggAl(6)*dconjg(ggB0(15))+ggBl(8)
     #*dconjg(ggA0(5))+ggAl(5)*dconjg(ggB0(8))+ggAl(2)*dconjg(ggB0(2))+g
     #gBl(2)*dconjg(ggA0(2)))+(16D0)*(+ggAl(10)*dconjg(ggA0(10))+ggBl(10
     #)*dconjg(ggB0(10))+ggAl(13)*dconjg(ggA0(13))+ggAl(14)*dconjg(ggA0(
     #14))+ggBl(14)*dconjg(ggB0(14))+ggBl(16)*dconjg(ggB0(16))+ggAl(7)*d
     #conjg(ggA0(7))+ggBl(7)*dconjg(ggB0(7))+ggAl(4)*dconjg(ggA0(4))+ggB
     #l(4)*dconjg(ggB0(4))+ggAl(3)*dconjg(ggA0(3))+ggBl(3)*dconjg(ggB0(3
     #))+ggAl(11)*dconjg(ggA0(11))+ggBl(11)*dconjg(ggB0(11))+ggAl(12)*dc
     #onjg(ggA0(12))+ggBl(12)*dconjg(ggB0(12))+ggAl(15)*dconjg(ggA0(15))
     #+ggBl(15)*dconjg(ggB0(15))+ggBl(8)*dconjg(ggB0(8))+ggAl(6)*dconjg(
     #ggA0(6))+ggBl(6)*dconjg(ggB0(6))+ggAl(5)*dconjg(ggA0(5))+ggAl(2)*d
     #conjg(ggA0(2))+ggBl(2)*dconjg(ggB0(2))))

      sqAf=4*((-2D0/3D0)*(+ggAf(10)*dconjg(ggB0(10))+ggBf(10)*dconjg(ggA
     #0(10))+ggAf(3)*dconjg(ggB0(3))+ggBf(3)*dconjg(ggA0(3))+ggAf(11)*dc
     #onjg(ggB0(11))+ggBf(11)*dconjg(ggA0(11))+ggAf(2)*dconjg(ggB0(2))+g
     #gBf(2)*dconjg(ggA0(2)))+(16D0/3D0)*(+ggAf(10)*dconjg(ggA0(10))+ggB
     #f(10)*dconjg(ggB0(10))+ggAf(3)*dconjg(ggA0(3))+ggBf(3)*dconjg(ggB0
     #(3))+ggAf(11)*dconjg(ggA0(11))+ggBf(11)*dconjg(ggB0(11))+ggAf(2)*d
     #conjg(ggA0(2))+ggBf(2)*dconjg(ggB0(2))))

      sqAh=4*((-2D0/3D0)*(+ggAh(10)*dconjg(ggB0(10))+ggBh(10)*dconjg(ggA
     #0(10))+ggAh(3)*dconjg(ggB0(3))+ggBh(3)*dconjg(ggA0(3))+ggAh(11)*dc
     #onjg(ggB0(11))+ggBh(11)*dconjg(ggA0(11))+ggAh(2)*dconjg(ggB0(2))+g
     #gBh(2)*dconjg(ggA0(2)))+(16D0/3D0)*(+ggAh(10)*dconjg(ggA0(10))+ggB
     #h(10)*dconjg(ggB0(10))+ggAh(3)*dconjg(ggA0(3))+ggBh(3)*dconjg(ggB0
     #(3))+ggAh(11)*dconjg(ggA0(11))+ggBh(11)*dconjg(ggB0(11))+ggAh(2)*d
     #conjg(ggA0(2))+ggBh(2)*dconjg(ggB0(2))))

      sqAr=4*((-2D0/9D0)*(+ggAr(10)*dconjg(ggB0(10))+ggBr(10)*dconjg(ggA
     #0(10))+ggAr(13)*dconjg(ggB0(16))+ggBr(16)*dconjg(ggA0(13))+ggAr(14
     #)*dconjg(ggB0(7))+ggBr(14)*dconjg(ggA0(7))+ggBr(7)*dconjg(ggA0(14)
     #)+ggAr(7)*dconjg(ggB0(14))+ggAr(4)*dconjg(ggB0(4))+ggBr(4)*dconjg(
     #ggA0(4))+ggAr(3)*dconjg(ggB0(3))+ggBr(3)*dconjg(ggA0(3))+ggAr(11)*
     #dconjg(ggB0(11))+ggBr(11)*dconjg(ggA0(11))+ggAr(12)*dconjg(ggB0(12
     #))+ggBr(12)*dconjg(ggA0(12))+ggAr(15)*dconjg(ggB0(6))+ggBr(15)*dco
     #njg(ggA0(6))+ggBr(6)*dconjg(ggA0(15))+ggAr(6)*dconjg(ggB0(15))+ggB
     #r(8)*dconjg(ggA0(5))+ggAr(5)*dconjg(ggB0(8))+ggAr(2)*dconjg(ggB0(2
     #))+ggBr(2)*dconjg(ggA0(2)))+(16D0/9D0)*(+ggAr(10)*dconjg(ggA0(10))
     #+ggBr(10)*dconjg(ggB0(10))+ggAr(13)*dconjg(ggA0(13))+ggAr(14)*dcon
     #jg(ggA0(14))+ggBr(14)*dconjg(ggB0(14))+ggBr(16)*dconjg(ggB0(16))+g
     #gAr(7)*dconjg(ggA0(7))+ggBr(7)*dconjg(ggB0(7))+ggAr(4)*dconjg(ggA0
     #(4))+ggBr(4)*dconjg(ggB0(4))+ggAr(3)*dconjg(ggA0(3))+ggBr(3)*dconj
     #g(ggB0(3))+ggAr(11)*dconjg(ggA0(11))+ggBr(11)*dconjg(ggB0(11))+ggA
     #r(12)*dconjg(ggA0(12))+ggBr(12)*dconjg(ggB0(12))+ggAr(15)*dconjg(g
     #gA0(15))+ggBr(15)*dconjg(ggB0(15))+ggBr(8)*dconjg(ggB0(8))+ggAr(6)
     #*dconjg(ggA0(6))+ggBr(6)*dconjg(ggB0(6))+ggAr(5)*dconjg(ggA0(5))+g
     #gAr(2)*dconjg(ggA0(2))+ggBr(2)*dconjg(ggB0(2))))

      sqAsl=4*((2D0)*(+ggAsl(10)*dconjg(ggA0(10))+ggAsl(10)*dconjg(ggB0(
     #10))+ggAsl(13)*dconjg(ggA0(13))+ggAsl(13)*dconjg(ggB0(16))+ggAsl(1
     #4)*dconjg(ggA0(14))+ggAsl(14)*dconjg(ggB0(7))+ggAsl(7)*dconjg(ggB0
     #(14))+ggAsl(7)*dconjg(ggA0(7))+ggAsl(4)*dconjg(ggA0(4))+ggAsl(4)*d
     #conjg(ggB0(4))+ggAsl(3)*dconjg(ggA0(3))+ggAsl(3)*dconjg(ggB0(3))+g
     #gAsl(11)*dconjg(ggA0(11))+ggAsl(11)*dconjg(ggB0(11))+ggAsl(12)*dco
     #njg(ggA0(12))+ggAsl(12)*dconjg(ggB0(12))+ggAsl(15)*dconjg(ggA0(15)
     #)+ggAsl(15)*dconjg(ggB0(6))+ggAsl(6)*dconjg(ggB0(15))+ggAsl(5)*dco
     #njg(ggB0(8))+ggAsl(6)*dconjg(ggA0(6))+ggAsl(5)*dconjg(ggA0(5))+ggA
     #sl(2)*dconjg(ggA0(2))+ggAsl(2)*dconjg(ggB0(2))))

      
