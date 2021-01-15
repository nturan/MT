/* $Modified: Mon May 14 21:00:31 2007 by puwer $ */
#ifdef __cplusplus
#include <complex>
#define XCOMPLEX std::complex<double>
extern "C" {
#else
#define XCOMPLEX _Complex double
#endif

// void boostx(p,q , pboost)
// void coup1x(sw2 , gw,gwwa,gwwz)
// void coup2x(sw2 , gal,gau,gad,gwf,gzn,gzl,gzu,gzd,g1)
// void coup3x(sw2,zmass,hmass , 
// void fsixxx(fi,sc,gc,fmass,fwidth , fsi)
// void fsoxxx(fo,sc,gc,fmass,fwidth , fso)
// void fvixxx(fi,vc,g,fmass,fwidth , fvi)
// void fvoxxx(fo,vc,g,fmass,fwidth , fvo)
// void ggggxx(wm,w31,wp,w32,g, vertex)
// void gggxxx(wm,wp,w3,g , vertex)
// void hioxxx(fi,fo,gc,smass,swidth , hio)                         
// void hvsxxx(vc,sc,g,smass,swidth , hvs)
// void hvvxxx(v1,v2,g,smass,swidth , hvv)
// void iovxxx(fi,fo,vc,g , vertex)
// void ixxxxx(
// void j3xxxx(fi,fo,gaf,gzf,zmass,zwidth , j3)
// void jgggxx(w1,w2,w3,g, jw3w)
// void jggxxx(v1,v2,g, jvv)
// void jioxxx(fi,fo,g,vmass,vwidth , jio)
// void jtioxx(fi,fo,g , jio)
// void jvsxxx(vc,sc,g,vmass,vwidth , jvs)
// void jvvxxx(v1,v2,g,vmass,vwidth , jvv)
// void jw3wxx(w1,w2,w3,g1,g2,wmass,wwidth,vmass,vwidth , jw3w)
// void jwwwxx(w1,w2,w3,gwwa,gwwz,zmass,zwidth,wmass,wwidth ,
// void oxxxxx(
// void rotxxx(p,q , prot)
// void vssxxx(vc,s1,s2,g , vertex)
// void vvsxxx(v1,v2,sc,g , vertex)
// void vvvxxx(wm,wp,w3,g , vertex)
  void vxxxxx_(const double p[4], const double & vmass, const int & nhel, 
	      const int & nsv, XCOMPLEX vc[6]);
// void w3w3xx(wm,w31,wp,w32,g31,g32,wmass,wwidth , vertex)
// void wwwwxx(wm1,wp1,wm2,wp2,gwwa,gwwz,zmass,zwidth , vertex)

#ifdef __cplusplus
}
#endif
