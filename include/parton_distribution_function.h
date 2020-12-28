/***************************************************************************//**
 * @author: Turan Nuraliyev
 ******************************************************************************/

#pragma once

#include "LHAPDF/LHAPDF.h"
#include <map>
#include <functional>
#include <string>

/***************************************************************************//**
 * class for calling parton distribution sets
 * 
 * Designed as a singleton class.  Can only have one instance.  Maybe obselete
 * design.  It uses LHAPDF.
 ******************************************************************************/

class PartonDistributionFunction {
  private:
    static bool hasInstance;  
    static bool isInitialized; 
    static PartonDistributionFunction *pdf;
    PartonDistributionFunction();  /*!<  private class constructor                        */
    /***************************************************************************
     * see LHAPDF documentation page for following members. 
     **************************************************************************/
    LHAPDF::PDFSet lhapdfset;
    LHAPDF::PDF* lhapdf;
  public:
    /***********************************************************************//**
     * Getter for the singleton instance.
     *
     * returns new PartonDistFunc object, if it is being called for the first 
     * time.  Returns the existing instance, falls it has been called before. 
     **************************************************************************/
    static PartonDistributionFunction* GetInstance();
    
    void InitPdfSet( std::string pdfname );
    std::map<std::string, std::function<double(double, double, double)>> Fs;
    double GetAlpha( double mur );
    double gg ( double x1, double x2, double muf2 );
    double qqb( double x1, double x2, double muf2 );
    double qbq( double x1, double x2, double muf2 );
    double gq ( double x1, double x2, double muf2 );
    double gqb( double x1, double x2, double muf2 );
    double qg ( double x1, double x2, double muf2 );
    double qbg( double x1, double x2, double muf2 );
    ~PartonDistributionFunction();
};

