/***************************************************************************//**
 * @author: Turan Nuraliyev
 ******************************************************************************/

#ifndef PARTONDISTFUNC_H
#define PARTONDISTFUNC_H

#include "LHAPDF/LHAPDF.h"

/***************************************************************************//**
 * class for calling parton distribution sets
 * 
 * Designed as a singleton class.  Can only have one instance.  Maybe obselete
 * design.  It uses LHAPDF.
 ******************************************************************************/

class PartonDistFunc {
  private:
    static bool hasInstance;  
    static bool isInitialized; 
    static PartonDistFunc *pdf;
    PartonDistFunc();  /*!<  private class constructor                        */
    /***************************************************************************
     * see LHAPDF documentation page for following members. 
     **************************************************************************/
    LHAPDF::PDFSet lhapdfset;  
    std::vector<LHAPDF::PDF*> lhapdfs;
    LHAPDF::PDF* lhapdf;
    LHAPDF::AlphaS_Analytic alphas_analytic;
  public:
    /***********************************************************************//**
     * Getter for the singleton instance.
     *
     * returns new PartonDistFunc object, if it is being called for the first 
     * time.  Returns the existing instance, falls it has been called before. 
     **************************************************************************/
    static PartonDistFunc* GetInstance();
    
    void InitPdfSet( std::string pdfname );
    void GetPdfValue( double x, double h[13] );
    double GetPdfValueById( double x, int id );
    double GetAlpha( double mur );
    double GetAnalyticAlpha( double mT, double lambda, int order,
                             double mur );
    double gg ( double x1, double x2, double muf2 );
    double qqb( double x1, double x2, double muf2 );
    double qbq( double x1, double x2, double muf2 );
    double gq ( double x1, double x2, double muf2 );
    double gqb( double x1, double x2, double muf2 );
    double qg ( double x1, double x2, double muf2 );
    double qbg( double x1, double x2, double muf2 );
    ~PartonDistFunc();
};

#endif //PARTONDISTFUNC_H
