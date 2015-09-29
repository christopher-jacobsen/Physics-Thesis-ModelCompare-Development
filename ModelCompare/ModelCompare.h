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

struct Observable;

////////////////////////////////////////////////////////////////////////////////

typedef void FillHistFunctionType(TH1D & hist, double weight, const HepMC::GenVertex & signal);
typedef std::function< FillHistFunctionType  > FillHistFunction;

typedef TH1D * TH1DFactoryFunctionType( const Observable & obs, const char * name, const char * title );
typedef std::function< TH1DFactoryFunctionType > TH1DFactoryFunction;

TH1DFactoryFunctionType DefaultTH1DFactory;
TH1DFactoryFunctionType DefaultTProfileFactory;

struct Observable
{
    const char *            name;
    const char *            title;
    Int_t                   nBins;
    Double_t                xMin;
    Double_t                xMax;
    const char *            xAxisTitle;
    const char *            yAxisTitle;
    FillHistFunction        fillFunction;
    TH1DFactoryFunction     factoryFunction = DefaultTH1DFactory;

    // force required fields to be filled on construction
    Observable( const char * name, const char * title, Int_t nBins, Double_t xMin, Double_t xMax,
                const char * xAxisTitle, const char * yAxisTitle,
                const FillHistFunction & fillFunction )
      : name(name), title(title), nBins(nBins), xMin(xMin), xMax(xMax),
        xAxisTitle(xAxisTitle), yAxisTitle(yAxisTitle),
        fillFunction(fillFunction)
    {
    }

    Observable( const char * name, const char * title, Int_t nBins, Double_t xMin, Double_t xMax,
                const char * xAxisTitle, const char * yAxisTitle,
                const FillHistFunction & fillFunction,
                const TH1DFactoryFunction & factoryFunction )
      : name(name), title(title), nBins(nBins), xMin(xMin), xMax(xMax),
        xAxisTitle(xAxisTitle), yAxisTitle(yAxisTitle),
        fillFunction(fillFunction), factoryFunction(factoryFunction)
    {
    }

    TH1D * MakeHist( const char * namePrefix = nullptr, const char * titlePrefix = nullptr,
                     const char * nameSuffix = nullptr, const char * titleSuffix = nullptr ) const
    {
        std::string sName  = BuildHistName(  namePrefix,  nameSuffix );
        std::string sTitle = BuildHistTitle( titlePrefix, titleSuffix );

        if (factoryFunction == nullptr)
            return DefaultTH1DFactory( *this, sName.c_str(), sTitle.c_str() );

        return factoryFunction( *this, sName.c_str(), sTitle.c_str() );
    }

    std::string BuildHistName(  const char * namePrefix  = nullptr, const char * nameSuffix  = nullptr ) const;
    std::string BuildHistTitle( const char * titlePrefix = nullptr, const char * titleSuffix = nullptr ) const;
};

typedef std::vector<Observable> ObservableVector;

////////////////////////////////////////////////////////////////////////////////

struct ModelFile
{
    const char *    fileName;
    const char *    modelName;
    const char *    modelTitle;
    double          crossSection;       // in pb
    double          crossSectionError;  // in pb
    size_t          crossSectionEvents;

    // force all fields to be set on construction
    ModelFile( const char * fileName, const char * modelName, const char * modelTitle,
               double crossSection, double crossSectionError, size_t crossSectionEvents )
      : fileName(fileName), modelName(modelName), modelTitle(modelTitle),
        crossSection(crossSection), crossSectionError(crossSectionError), crossSectionEvents(crossSectionEvents)
    {
    }
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

struct GoodBadHists
{
    RootUtil::TH1DUniquePtr     good;
    RootUtil::TH1DUniquePtr     bad;
};

////////////////////////////////////////////////////////////////////////////////

void ScaleHistToLuminosity( double luminosity, const RootUtil::TH1DVector & hists, const ModelFile & eventFile, bool bApplyCrossSectionError = false );

GoodBadHists HistSplitGoodBadBins( const TH1D * pSource, const TH1D * pCompare = nullptr );
std::list<GoodBadHists> HistSplitGoodBadBins( const RootUtil::ConstTH1DVector & hists, const RootUtil::ConstTH1DVector & compare );

void WriteCompareFigure( const char * name, const char * title,
                         const RootUtil::ConstTH1DVector & data, const RootUtil::ConstTH1DVector & compare,
                         const RootUtil::ColorVector & dataColors,
                         const RootUtil::ConstTH1DVector & rawData );

void LoadHistData( const ModelFileVector & models, const ObservableVector & observables, std::vector<RootUtil::TH1DVector> & hists );

void CalculateCompareHists( const Observable & obs, const RootUtil::ConstTH1DVector & data, RootUtil::TH1DVector & comp,
                            const ModelFileVector & models, const RootUtil::ColorVector & dataColors );

void ModelCompare( const char * outputFileName, const ModelFileVector & models, const ObservableVector & observables, const FigureSetupVector & figures );

////////////////////////////////////////////////////////////////////////////////

}  // namespace ModelCompare

////////////////////////////////////////////////////////////////////////////////

#endif // MODEL_COMPARE_H
