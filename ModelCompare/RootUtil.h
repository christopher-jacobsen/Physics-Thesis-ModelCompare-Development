//
//  RootUtil.h
//
//  Created by Christopher Jacobsen on 06/09/15.
//  Copyright (c) 2015 Christopher Jacobsen. All rights reserved.
//

#ifndef ROOT_UTIL_H
#define ROOT_UTIL_H

#include "common.h"

// Root includes
#include <Rtypes.h>
#include <TLorentzVector.h>

// HepMC includes
#include <HepMC/SimpleVector.h>

////////////////////////////////////////////////////////////////////////////////
// forward declarations

class TH1D;
class TFile;

namespace HepMC
{
class GenVertex;
class GenParticle;
}

////////////////////////////////////////////////////////////////////////////////

namespace RootUtil
{

////////////////////////////////////////////////////////////////////////////////

typedef std::vector<const char *>   CStringVector;

typedef std::vector<const TH1D *>       ConstTH1DVector;
typedef std::vector<TH1D *>             TH1DVector;
typedef std::vector<Color_t>            ColorVector;

typedef std::unique_ptr<TH1D>           TH1DUniquePtr;
typedef std::unique_ptr<const TH1D>     ConstTH1DUniquePtr;

typedef std::vector< const HepMC::GenParticle * > ConstGenParticleVector;

////////////////////////////////////////////////////////////////////////////////

inline ConstTH1DVector ToConstTH1DVector( const TH1DVector & v )
{
    return ConstTH1DVector( v.cbegin(), v.cend() );
}

inline TLorentzVector ToLorentz( const HepMC::FourVector & v )
{
    return TLorentzVector( v.x(), v.y(), v.z(), v.t() );
}

////////////////////////////////////////////////////////////////////////////////

void LoadEvents( const char * eventFileName, std::function<void(const HepMC::GenVertex & signal)> EventFunc );

ConstGenParticleVector     FindOutgoingParticles(      const HepMC::GenVertex & signal, int pdg, bool bThrowNotFound = true );
const HepMC::GenParticle * FindSingleOutgoingParticle( const HepMC::GenVertex & signal, int pdg, bool bThrowNotFound = true );

double GetObsPT(   const HepMC::GenVertex & signal, int pdg );
double GetObsRap(  const HepMC::GenVertex & signal, int pdg );
double GetObsEta(  const HepMC::GenVertex & signal, int pdg );
double GetObsPhi(  const HepMC::GenVertex & signal, int pdg );
double GetObsMass( const HepMC::GenVertex & signal, int pdg1, int pdg2 );

void FillHistPT(   TH1D & hist, double weight, const HepMC::GenVertex & signal, int pdg );
void FillHistRap(  TH1D & hist, double weight, const HepMC::GenVertex & signal, int pdg );
void FillHistEta(  TH1D & hist, double weight, const HepMC::GenVertex & signal, int pdg );
void FillHistPhi(  TH1D & hist, double weight, const HepMC::GenVertex & signal, int pdg );
void FillHistMass( TH1D & hist, double weight, const HepMC::GenVertex & signal, int pdg1, int pdg2 );

////////////////////////////////////////////////////////////////////////////////

void LogMsgHistUnderOverflow( const TH1D & hist );
void LogMsgHistUnderOverflow( const ConstTH1DVector & hists );

void LogMsgHistStats( const TH1D & hist );

void LogMsgHistEffectiveEntries( const TH1D & hist );
void LogMsgHistEffectiveEntries( const ConstTH1DVector & hists );

void LogMsgHistDump( const TH1D & hist );

////////////////////////////////////////////////////////////////////////////////

TH1D * ConvertTProfileToTH1D( const TH1D * pProfile, bool bDeleteProfile );

bool IsHistSumw2Enabled( const TH1D & hist );

void SetupHist( TH1D & hist, const char * xAxisTitle = nullptr, const char * yAxisTitle = nullptr,
                Color_t lineColor = -1, Color_t markerColor = -1, Color_t fillColor = -1 );

void ScaleHistToLuminosity( double luminosity, const TH1D & hist, size_t nEvents, double crossSection,
                            double crossSectionError, bool bApplyCrossSectionError = false );

void ScaleHistToLuminosity( double luminosity, const TH1DVector & hist, size_t nEvents, double crossSection,
                            double crossSectionError, bool bApplyCrossSectionError = false );

void WriteHists( TFile * pFile, const TH1DVector & hists );

////////////////////////////////////////////////////////////////////////////////

void GetHistDrawMinMax( const TH1D & hist,             Double_t & ymin, Double_t & ymax );
void GetHistDrawMinMax( const ConstTH1DVector & hists, Double_t & ymin, Double_t & ymax );

TH1DVector DrawMultipleHist( const char * title, const ConstTH1DVector & hists, const ColorVector & colors = {}, const CStringVector drawOptions = {} );

////////////////////////////////////////////////////////////////////////////////

Double_t GetHistBinEffectiveEntries( const TH1D & hist, Int_t bin );

size_t HistNonEmptyBinCount( const TH1D & hist, bool bIncludeUnderOverflow = false );
size_t HistNonEmptyBinCount( const TH1D & h1, const TH1D & h2, bool bCountEitherNonEmpty = false, bool bIncludeUnderOverflow = false );

void ZeroHistBin( TH1D & hist, Int_t bin );
void ZeroHistEmptyBins( TH1D & h1, TH1D & h2 );

////////////////////////////////////////////////////////////////////////////////

Double_t KolmogorovTest_NonEmptyBins( const TH1D & h1, const TH1D & h2 );

struct Chi2Result
{
    Double_t chi2     = 0;
    Int_t    ndf      = 0;
    Int_t    igood    = 0;
    Double_t prob     = 0;
    Double_t chi2_ndf = 0;

    void Chi2Test( const TH1D & h1, const TH1D & h2 );

    std::string Label();
};

////////////////////////////////////////////////////////////////////////////////

}  // namespace RootUtil

#endif // ROOT_UTIL_H
