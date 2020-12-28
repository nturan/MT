#include "parton_distribution_function.h"

using namespace std;

bool PartonDistributionFunction::hasInstance = false;
bool PartonDistributionFunction::isInitialized = false;
PartonDistributionFunction* PartonDistributionFunction::pdf = NULL;


PartonDistributionFunction::PartonDistributionFunction () {
  hasInstance = true;
}

PartonDistributionFunction::~PartonDistributionFunction () {
  hasInstance = false;
}

PartonDistributionFunction* PartonDistributionFunction::GetInstance(){
  if (!hasInstance) {
    pdf = new PartonDistributionFunction();
  }
  return pdf;
}

void PartonDistributionFunction::InitPdfSet(string pdfname){
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
  using namespace std::placeholders;
  Fs["gg"]  = std::bind(&PartonDistributionFunction::gg,  this,_1, _2, _3);
  Fs["qqb"] = std::bind(&PartonDistributionFunction::qqb, this, _1, _2, _3);
  Fs["qbq"] = std::bind(&PartonDistributionFunction::qbq, this, _1, _2, _3);
  Fs["gq"]  = std::bind(&PartonDistributionFunction::gq,  this, _1, _2, _3);
  Fs["qg"]  = std::bind(&PartonDistributionFunction::qg,  this, _1, _2, _3);
  Fs["gqb"] = std::bind(&PartonDistributionFunction::gqb, this, _1, _2, _3);
  Fs["qbg"] = std::bind(&PartonDistributionFunction::qbg, this, _1, _2, _3);
}


double PartonDistributionFunction::GetAlpha( double mur ){
  if (!isInitialized) { 
    cout << "PDF not initialized" << endl;
    exit(1);
  }
      return lhapdf->alphasQ( mur );
}	

double PartonDistributionFunction::gg( double x1, double x2, double muf2 ){
  if (!isInitialized) { 
    cout << "PDF not initialized" << endl;
    exit(1);
  }
  if ( x1>1.0 ){return 0.0; }
  if ( x2>1.0 ){return 0.0; }
  return lhapdf->xfxQ2(21, x1, muf2)*lhapdf->xfxQ2(21, x2, muf2)/x1/x2;
}

double PartonDistributionFunction::qqb( double x1, double x2, double muf2 ){
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

double PartonDistributionFunction::qbq( double x1, double x2, double muf2 ){
  return qqb( x2, x1, muf2 );
}

double PartonDistributionFunction::gq( double x1, double x2, double muf2 ){
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

double PartonDistributionFunction::qg( double x1, double x2, double muf2 ){
  return gq( x2, x1, muf2 );
}

double PartonDistributionFunction::gqb( double x1, double x2, double muf2 ){
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

double PartonDistributionFunction::qbg( double x1, double x2, double muf2 ){
  return gqb( x2, x1, muf2 );
}
