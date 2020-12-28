#include "parton_distribution_function.h"

using namespace std;

bool PartonDistFunc::hasInstance = false;
bool PartonDistFunc::isInitialized = false;
PartonDistFunc* PartonDistFunc::pdf = NULL;


PartonDistFunc::PartonDistFunc () {
  hasInstance = true;
}

PartonDistFunc::~PartonDistFunc () {
  hasInstance = false;
}

PartonDistFunc* PartonDistFunc::GetInstance(){
  if (!hasInstance) {
    pdf = new PartonDistFunc();
  }
  return pdf;
}

void PartonDistFunc::InitPdfSet(string pdfname){
  try {
    lhapdfset = LHAPDF::getPDFSet(pdfname);
  } catch (const LHAPDF::ReadError &) {
    cerr << "Cannot load PDF set from LHAPDF." << endl;
    cerr << "Please check the PDF name " << pdfname << " and"
         << endl;
    cerr << "ensure its presence in your LHAPDF setup." << endl;
    exit(1);
  }
  isInitialized = true;
  lhapdf = lhapdfset.mkPDF(0);
}

double PartonDistFunc::GetAlpha( double mur ){
  if (!isInitialized) { 
    cout << "PDF not initialized" << endl;
    exit(1);
  }
      return lhapdf->alphasQ( mur );
}	

double PartonDistFunc::gg( double x1, double x2, double muf2 ){
  if (!isInitialized) { 
    cout << "PDF not initialized" << endl;
    exit(1);
  }
  if ( x1>1.0 ){return 0.0; }
  if ( x2>1.0 ){return 0.0; }
  return lhapdf->xfxQ2(21, x1, muf2)*lhapdf->xfxQ2(21, x2, muf2)/x1/x2;
}

double PartonDistFunc::qqb( double x1, double x2, double muf2 ){
  if (!isInitialized) { 
    cout << "PDF not initialized" << endl;
    exit(1);
  }
  if ( x1>1.0 ){return 0.0; }
  if ( x2>1.0 ){return 0.0; }
  return ( lhapdf->xfxQ2( + 1, x1, muf2)*lhapdf->xfxQ2( - 1, x2, muf2) + 
           lhapdf->xfxQ2( + 2, x1, muf2)*lhapdf->xfxQ2( - 2, x2, muf2) + 
           lhapdf->xfxQ2( + 3, x1, muf2)*lhapdf->xfxQ2( - 3, x2, muf2) + 
           lhapdf->xfxQ2( + 4, x1, muf2)*lhapdf->xfxQ2( - 4, x2, muf2) + 
           lhapdf->xfxQ2( + 5, x1, muf2)*lhapdf->xfxQ2( - 5, x2, muf2) )/x1/x2;
}

double PartonDistFunc::qbq( double x1, double x2, double muf2 ){
  return qqb( x2, x1, muf2 );
}

double PartonDistFunc::gq( double x1, double x2, double muf2 ){
  if (!isInitialized) { 
    cout << "PDF not initialized" << endl;
    exit(1);
  }
  if ( x1>1.0 ){return 0.0; }
  if ( x2>1.0 ){return 0.0; }

  return ( lhapdf->xfxQ2(21, x1, muf2)*lhapdf->xfxQ2( + 1, x2, muf2) + 
           lhapdf->xfxQ2(21, x1, muf2)*lhapdf->xfxQ2( + 2, x2, muf2) + 
           lhapdf->xfxQ2(21, x1, muf2)*lhapdf->xfxQ2( + 3, x2, muf2) + 
           lhapdf->xfxQ2(21, x1, muf2)*lhapdf->xfxQ2( + 4, x2, muf2) + 
           lhapdf->xfxQ2(21, x1, muf2)*lhapdf->xfxQ2( + 5, x2, muf2) )/x1/x2;
}

double PartonDistFunc::qg( double x1, double x2, double muf2 ){
  return gq( x2, x1, muf2 );
}

double PartonDistFunc::gqb( double x1, double x2, double muf2 ){
  if (!isInitialized) { 
    cout << "PDF not initialized" << endl;
    exit(1);
  }
  if ( x1>1.0 ){return 0.0; }
  if ( x2>1.0 ){return 0.0; }
  return ( lhapdf->xfxQ2(21, x1, muf2)*lhapdf->xfxQ2( - 1, x2, muf2) + 
           lhapdf->xfxQ2(21, x1, muf2)*lhapdf->xfxQ2( - 2, x2, muf2) + 
           lhapdf->xfxQ2(21, x1, muf2)*lhapdf->xfxQ2( - 3, x2, muf2) + 
           lhapdf->xfxQ2(21, x1, muf2)*lhapdf->xfxQ2( - 4, x2, muf2) + 
           lhapdf->xfxQ2(21, x1, muf2)*lhapdf->xfxQ2( - 5, x2, muf2) )/x1/x2;
}

double PartonDistFunc::qbg( double x1, double x2, double muf2 ){
  return gqb( x2, x1, muf2 );
}
