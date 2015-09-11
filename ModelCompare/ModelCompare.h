//
//  ModelCompare.h
//  ModelCompare
//
//  Created by Christopher Jacobsen on 11/09/15.
//  Copyright (c) 2015 Christopher Jacobsen. All rights reserved.
//

#ifndef MODEL_COMPARE_H
#define MODEL_COMPARE_H

#include "common.h"
#include "RootUtil.h"

// Root includes
#include <Rtypes.h>

////////////////////////////////////////////////////////////////////////////////
// forward declarations

class TH1D;

namespace HepMC
{
class GenVertex;
}

////////////////////////////////////////////////////////////////////////////////

namespace ModelCompare
{

////////////////////////////////////////////////////////////////////////////////

typedef std::function<void(const HepMC::GenVertex & signal, TH1D & hist)> FillHistFunction;

struct Observable
{
    const char *        name;
    const char *        title;
    Int_t               nBins;
    Double_t            min;
    Double_t            max;
    const char *        xAxisTitle;
    const char *        yAxisTitle;
    FillHistFunction    fillFunction;
};

typedef std::vector<Observable> ObservableVector;

////////////////////////////////////////////////////////////////////////////////

struct ModelFile
{
    const char *    fileName;
    const char *    modelName;
    const char *    modelTitle;
    double          crossSection;  // in pb
};

typedef std::vector<ModelFile>  ModelFileVector;

////////////////////////////////////////////////////////////////////////////////

struct FigureSetup
{
  //const char *    name;   // TODO
  //const char *    title;  // TODO

    RootUtil::CStringVector     modelNames;
    double                      luminosity  = 0.0;     // in fb^-1
    RootUtil::ColorVector       colors      = DefaultColors;

    FigureSetup() = default;
    FigureSetup( const RootUtil::CStringVector & n )                                            : modelNames(n)                             {}
    FigureSetup( const RootUtil::CStringVector & n, double l )                                  : modelNames(n), luminosity(l)              {}
    FigureSetup( const RootUtil::CStringVector & n, double l, const RootUtil::ColorVector & c ) : modelNames(n), luminosity(l), colors(c)   {}

    static const RootUtil::ColorVector DefaultColors;
};

typedef std::vector<FigureSetup> FigureSetupVector;

////////////////////////////////////////////////////////////////////////////////

void WriteCompareFigure( const char * name, const char * title,
                         const RootUtil::ConstTH1DVector & data, RootUtil::ConstTH1DVector & compare,
                         const RootUtil::ColorVector & dataColors );

void LoadHistData( const ModelFileVector & models, const ObservableVector & observables, std::vector<RootUtil::TH1DVector> & hists );

void CalculateCompareHists( const Observable & obs, const RootUtil::ConstTH1DVector & data, RootUtil::TH1DVector & comp,
                            const ModelFileVector & models, const RootUtil::ColorVector & dataColors );

void ModelCompare( const char * outputFileName, const ModelFileVector & models, const ObservableVector & observables, const FigureSetupVector & figures );

////////////////////////////////////////////////////////////////////////////////

}  // namespace ModelCompare

////////////////////////////////////////////////////////////////////////////////

#endif // MODEL_COMPARE_H
