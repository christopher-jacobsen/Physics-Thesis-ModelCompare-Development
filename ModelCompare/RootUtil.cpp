//
//  RootUtil.cpp
//
//  Created by Christopher Jacobsen on 06/09/15.
//  Copyright (c) 2015 Christopher Jacobsen. All rights reserved.
//

#include "RootUtil.h"
#include "common.h"

// Root includes
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TProfile.h>
#include <TNtupleD.h>
#include <TObjArray.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TCanvas.h>
#include <THistPainter.h>

// HepMC includes
#include <HepMC/IO_GenEvent.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>

// Sherpa includes
#include "Gzip_Stream.H"

namespace RootUtil
{

////////////////////////////////////////////////////////////////////////////////

// gain access to TProfile fBinEntries
struct MyProfile : public TProfile
{
    using TProfile::fBinEntries;
};

////////////////////////////////////////////////////////////////////////////////
void LoadEvents( const char * eventFileName, std::function<void(const HepMC::GenVertex & signal)> EventFunc, size_t maxEvents /*= 0*/ )
{
    std::unique_ptr<ATOOLS::igzstream>  upStream;
    std::unique_ptr<HepMC::IO_GenEvent> upInput;

    try
    {
        LogMsgInfo( "Input file: %hs", FMT_HS(eventFileName) );

        upStream.reset( new ATOOLS::igzstream( eventFileName, std::ios::in ) );

        upInput.reset( new HepMC::IO_GenEvent( *upStream.get() ) );
    }
    catch (...)
    {
        LogMsgError( "Failed to construct HepMC IO object for file (%hs).", FMT_HS(eventFileName) );
        throw;
    }

    if (maxEvents == 0)
        maxEvents = std::numeric_limits<size_t>::max();

    HepMC::GenEvent genEvent;

    size_t nEvents = 0;
    while ((nEvents < maxEvents) && upInput->fill_next_event( &genEvent ))
    {
        const HepMC::GenVertex * pSignal = genEvent.signal_process_vertex();
        if (!pSignal)
            ThrowError( "Missing signal vertex for event." );

        EventFunc( *pSignal );

        ++nEvents;
    }
}

////////////////////////////////////////////////////////////////////////////////
ConstGenParticleVector FindOutgoingParticles( const HepMC::GenVertex & signal, int pdg, bool bThrowNotFound /*= true*/ )
{
    ConstGenParticleVector result;

    auto itrPart = signal.particles_out_const_begin();
    auto endPart = signal.particles_out_const_end();
    for ( ; itrPart != endPart; ++itrPart)
    {
        const HepMC::GenParticle * pPart = *itrPart;
        if (pPart && (pPart->pdg_id() == pdg))
            result.push_back( pPart );
    }

    if (bThrowNotFound && result.empty())
        ThrowError( "No outgoing particle with pdg code: " + std::to_string(pdg) );

    return result;
}

////////////////////////////////////////////////////////////////////////////////
const HepMC::GenParticle * FindSingleOutgoingParticle( const HepMC::GenVertex & signal, int pdg, bool bThrowNotFound /*= true*/ )
{
    const HepMC::GenParticle * result = nullptr;

    auto itrPart = signal.particles_out_const_begin();
    auto endPart = signal.particles_out_const_end();
    for ( ; itrPart != endPart; ++itrPart)
    {
        const HepMC::GenParticle * pPart = *itrPart;
        if (pPart && (pPart->pdg_id() == pdg))
        {
            if (result == nullptr)
                result = pPart;
            else
                ThrowError( "Multiple outgoing particles with pdg code: " + std::to_string(pdg) );
        }
    }

    if (bThrowNotFound && !result)
        ThrowError( "No outgoing particle with pdg code: " + std::to_string(pdg) );

    return result;
}

////////////////////////////////////////////////////////////////////////////////
double GetObsPT( const HepMC::GenVertex & signal, int pdg )
{
    const HepMC::GenParticle * pPart = FindSingleOutgoingParticle( signal, pdg );

    TLorentzVector vec = ToLorentz( pPart->momentum() );

    return vec.Pt();
}

////////////////////////////////////////////////////////////////////////////////
double GetObsRap( const HepMC::GenVertex & signal, int pdg )
{
    const HepMC::GenParticle * pPart = FindSingleOutgoingParticle( signal, pdg );

    TLorentzVector vec = ToLorentz( pPart->momentum() );

    return vec.Rapidity();
}

////////////////////////////////////////////////////////////////////////////////
double GetObsEta( const HepMC::GenVertex & signal, int pdg )
{
    const HepMC::GenParticle * pPart = FindSingleOutgoingParticle( signal, pdg );

    TLorentzVector vec = ToLorentz( pPart->momentum() );

    return vec.Eta();
}

////////////////////////////////////////////////////////////////////////////////
double GetObsPhi( const HepMC::GenVertex & signal, int pdg )
{
    const HepMC::GenParticle * pPart = FindSingleOutgoingParticle( signal, pdg );

    TLorentzVector vec = ToLorentz( pPart->momentum() );

    return vec.Phi();
}

////////////////////////////////////////////////////////////////////////////////
double GetObsMass( const HepMC::GenVertex & signal, int pdg1, int pdg2 )
{
    const HepMC::GenParticle * pPart1 = FindSingleOutgoingParticle( signal, pdg1 );
    const HepMC::GenParticle * pPart2 = FindSingleOutgoingParticle( signal, pdg2 );

    TLorentzVector vec1 = ToLorentz( pPart1->momentum() );
    TLorentzVector vec2 = ToLorentz( pPart2->momentum() );
    TLorentzVector vec  = vec1 + vec2;

    return vec.M();
}

////////////////////////////////////////////////////////////////////////////////
void FillHistPT( TH1D & hist, double weight, const HepMC::GenVertex & signal, int pdg )
{
    hist.Fill( GetObsPT( signal, pdg ), weight );
}

////////////////////////////////////////////////////////////////////////////////
void FillHistRap( TH1D & hist, double weight, const HepMC::GenVertex & signal, int pdg )
{
    hist.Fill( GetObsRap( signal, pdg ), weight );
}

////////////////////////////////////////////////////////////////////////////////
void FillHistEta( TH1D & hist, double weight, const HepMC::GenVertex & signal, int pdg )
{
    hist.Fill( GetObsEta( signal, pdg ), weight );
}

////////////////////////////////////////////////////////////////////////////////
void FillHistPhi( TH1D & hist, double weight, const HepMC::GenVertex & signal, int pdg )
{
    hist.Fill( GetObsPhi( signal, pdg ), weight );
}

////////////////////////////////////////////////////////////////////////////////
void FillHistMass( TH1D & hist, double weight, const HepMC::GenVertex & signal, int pdg1, int pdg2 )
{
    hist.Fill( GetObsMass( signal, pdg1, pdg2 ), weight );
}

////////////////////////////////////////////////////////////////////////////////
void LogMsgHistUnderOverflow( const TH1D & hist )
{
    Double_t underflow = hist.GetBinContent(0);
    Double_t overflow  = hist.GetBinContent( hist.GetNbinsX() + 1 );

    if ((underflow != 0) || (overflow != 0))
        LogMsgInfo( "%hs: under|overflow = %g | %g", FMT_HS(hist.GetName()), FMT_F(underflow), FMT_F(overflow) );
}

////////////////////////////////////////////////////////////////////////////////
void LogMsgHistUnderOverflow( const ConstTH1DVector & hists )
{
    for (const TH1D * pHist : hists)
    {
        if (pHist)
            LogMsgHistUnderOverflow( *pHist );
    }
}

////////////////////////////////////////////////////////////////////////////////
void LogMsgHistStats( const TH1D & hist )
{
   // size of statistics data (size of array used in GetStats()/ PutStats)
   // s[0]  = sumw       s[1]  = sumw2
   // s[2]  = sumwx      s[3]  = sumwx2                     
   // s[4]  = sumwy      s[5]  = sumwy2   s[6]  = sumwxy
   // s[7]  = sumwz      s[8]  = sumwz2   s[9]  = sumwxz   s[10]  = sumwyz  
   // s[11] = sumwt      s[12] = sumwt2  (11 and 12 used only by TProfile3D)

    Double_t stats[TH1::kNstat] = { };
    hist.GetStats( stats );

    LogMsgInfo( "%hs: sumw=%g sumw2=%g sumwx=%g sumwx2=%g sumwy=%g sumwy2=%g", FMT_HS(hist.GetName()),
                FMT_F(stats[0]), FMT_F(stats[1]), FMT_F(stats[2]), FMT_F(stats[3]),     // TH1D and TProfile
                FMT_F(stats[4]), FMT_F(stats[5]) );                                     // TProfile
}

////////////////////////////////////////////////////////////////////////////////
void LogMsgHistEffectiveEntries( const TH1D & hist )
{
    LogMsgInfo( "%hs: entries = %g, eff. entries = %g, sum bins = %g",
        FMT_HS(hist.GetName()),
        FMT_F(hist.GetEntries()), FMT_F(hist.GetEffectiveEntries()),
        FMT_F(hist.GetSumOfWeights()) );
}

////////////////////////////////////////////////////////////////////////////////
void LogMsgHistEffectiveEntries( const ConstTH1DVector & hists )
{
    for (const TH1D * pHist : hists)
    {
        if (pHist)
            LogMsgHistEffectiveEntries( *pHist );
    }
}

////////////////////////////////////////////////////////////////////////////////
void LogMsgHistBinCounts( const TH1D & hist )
{
    Int_t nBinsA     = hist.GetNbinsX();
    Int_t nBinsB     = hist.GetSize();
    Int_t nNonEmptyA = (Int_t)HistNonEmptyBinCount(hist, false);
    Int_t nNonEmptyB = (Int_t)HistNonEmptyBinCount(hist, true );
    Int_t nErrorA    = (Int_t)HistErrorBinCount(hist, false);
    Int_t nErrorB    = (Int_t)HistErrorBinCount(hist, true);

    LogMsgInfo( "%hs:\tbins=%i(%i)  non-empty=%i(%i)  errors=%i(%i)", FMT_HS(hist.GetName()),
                FMT_I(nBinsA),     FMT_I(nBinsB),
                FMT_I(nNonEmptyA), FMT_I(nNonEmptyB),
                FMT_I(nErrorA),    FMT_I(nErrorB) );
}

////////////////////////////////////////////////////////////////////////////////
void LogMsgHistBinCounts( const TH1D & hist1, const TH1D & hist2, bool bCountUnion /*= false*/ )
{
    Int_t nBinsA     = hist1.GetNbinsX();
    Int_t nBinsB     = hist1.GetSize();
    Int_t nNonEmptyA = (Int_t)HistNonEmptyBinCount( hist1, hist2, bCountUnion, false );
    Int_t nNonEmptyB = (Int_t)HistNonEmptyBinCount( hist1, hist2, bCountUnion, true  );
    Int_t nErrorA    = (Int_t)HistErrorBinCount(    hist1, hist2, bCountUnion, false );
    Int_t nErrorB    = (Int_t)HistErrorBinCount(    hist1, hist2, bCountUnion, true  );

    LogMsgInfo( "%hs %hs %hs:\tbins=%i(%i)  non-empty=%i(%i)  errors=%i(%i)",
                FMT_HS(hist1.GetName()), FMT_HS(bCountUnion ? "||" : "&&" ), FMT_HS(hist2.GetName()),
                FMT_I(nBinsA),     FMT_I(nBinsB),
                FMT_I(nNonEmptyA), FMT_I(nNonEmptyB),
                FMT_I(nErrorA),    FMT_I(nErrorB) );
}

////////////////////////////////////////////////////////////////////////////////
void LogMsgHistBinCounts( const ConstTH1DVector & hists )
{
    for (const TH1D * pHist : hists)
    {
        if (pHist)
            LogMsgHistBinCounts(*pHist);
    }
}

////////////////////////////////////////////////////////////////////////////////
void LogMsgHistBinCounts( const ConstTH1DVector & hists1, const ConstTH1DVector & hists2, bool bCountUnion /*= false*/ )
{
    auto       itr1 = hists1.cbegin();
    const auto end1 = hists1.cend();

    auto       itr2 = hists2.cbegin();
    const auto end2 = hists2.cend();

    for ( ; (itr1 != end1) && (itr2 != end2); ++itr1, ++itr2)
    {
        if (*itr1 && *itr2)
            LogMsgHistBinCounts( **itr1, **itr2, bCountUnion );
    }
}

////////////////////////////////////////////////////////////////////////////////
void LogMsgHistDump( const TH1D & hist )
{
    const MyProfile * pProf = hist.InheritsFrom(TProfile::Class()) ? static_cast<const MyProfile *>(&hist) : nullptr;

    const Double_t * pSumw        = hist.GetArray();
    const Double_t * pSumw2       = hist.GetSumw2()->GetArray(); // can be null
    const Double_t * pBinEntries  = pProf ? pProf->fBinEntries.GetArray()    : nullptr;
    const Double_t * pBinSumw2    = pProf ? pProf->GetBinSumw2()->GetArray() : nullptr;

    const Int_t nSize = hist.GetSize();
    for (Int_t bin = 0; bin < nSize; ++bin)
    {
        Double_t sumw  = pSumw[bin];
        Double_t sumw2 = pSumw2 ? pSumw2[bin] : sumw;

        if (!pProf)
        {
            Double_t error = hist.GetBinError(bin);
            Double_t nEff  = (sumw2 != 0 ? sumw * sumw / sumw2 : 0);

            LogMsgInfo( "%i: sumw=%.13E  sumw2=%.13E  error=%.13E  nEff=%.13E",
                        FMT_I(bin), FMT_F(sumw), FMT_F(sumw2), FMT_F(error), FMT_F(nEff) );
        }
        else
        {
            Double_t binEntries = pBinEntries[bin];
            Double_t binSumw2   = pBinSumw2 ? pBinSumw2[bin] : binEntries;

            Double_t content = pProf->GetBinContent(bin);
            Double_t error   = pProf->GetBinError(bin);
            Double_t nEff    = pProf->GetBinEffectiveEntries(bin);

            LogMsgInfo( "%i: sumw=%.13E  sumw2=%.13E  bEnt=%.13E  bSw2=%.13E  cnt=%.13E  err=%.13E  nEff=%.13E",
                        FMT_I(bin), FMT_F(sumw), FMT_F(sumw2), FMT_F(binEntries), FMT_F(binSumw2),
                        FMT_F(content), FMT_F(error), FMT_F(nEff) );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
TH1D * ConvertTProfileToTH1D( const TH1D * pProfile, bool bDeleteProfile )
{
    TH1D * pHist = nullptr;

    if (pProfile->InheritsFrom(TProfile::Class()))
    {
        // create a new TH1D from the TProfile
        pHist = ((TProfile *)pProfile)->ProjectionX( pProfile->GetName() );

        SetupHist( *pHist );

        if (bDeleteProfile)
            delete pProfile;
    }
    else
    {
        // clone the TH1D

        if (bDeleteProfile)
            pHist = const_cast<TH1D *>(pProfile);   // quick clone: just return the input object
        else
            pHist = (TH1D *)pProfile->Clone();
    }

    pHist->SetDirectory( nullptr );     // ensure not owned by any directory

    return pHist;
}

////////////////////////////////////////////////////////////////////////////////
bool IsHistSumw2Enabled( const TH1D & hist )
{
    if (hist.InheritsFrom(TProfile::Class()))
        return static_cast<const TProfile *>(&hist)->GetBinSumw2()->fN != 0;

    return hist.GetSumw2()->fN != 0;
}

////////////////////////////////////////////////////////////////////////////////
void SetupHist( TH1D & hist, const char * xAxisTitle, const char * yAxisTitle,
                Color_t lineColor /*= -1*/, Color_t markerColor /*= -1*/, Color_t fillColor /*= -1*/ )
{
    if (!IsHistSumw2Enabled(hist))
        hist.Sumw2();

    hist.SetStats( kFALSE );

    if (xAxisTitle)
        hist.GetXaxis()->SetTitle( xAxisTitle );
    if (yAxisTitle)
        hist.GetYaxis()->SetTitle( yAxisTitle );

    if (lineColor >= 0)
        hist.SetLineColor( lineColor );
    if (markerColor >= 0)
        hist.SetMarkerColor( markerColor );
    if (fillColor >= 0)
        hist.SetFillColor( fillColor );
}

////////////////////////////////////////////////////////////////////////////////
void ScaleHistToLuminosity( double luminosity, TH1D & hist, size_t nEvents, double crossSection,
                            double crossSectionError, bool bApplyCrossSectionError /*= false*/ )
{
    // Luminosity scaling is not the same as scaling by a constant.
    // Instead we want to change the number of effective entries:
    // All internal sums must be scaled in like fashion, as if the number
    // of entries in each sum was scaled.

    // --- TH1D ----
    // internal sums:       sumw, sumw2
    // effective entries:   sumw^2 / sumw2
    // scaling by s:        binContent *= s, binError *= sqrt(s)

    // --- TProfile ----
    // internal sums:       sumw, sumw2, binEntries, binSumw2
    // effective entries:   binEntries^2 / binSumw2
    // scaling by s:        binError *= 1/sqrt(s)

    double scale = luminosity * crossSection * 1000 / nEvents;

    LogMsgInfo( "Scaling %hs with %g", FMT_HS(hist.GetName()), FMT_F(scale) );

    //LogMsgInfo( "------ before luminosity scale ------" );
    //hist.Print("all");
    //LogMsgHistStats(hist);
    //LogMsgHistEffectiveEntries(hist);

    MyProfile * pProf = hist.InheritsFrom(TProfile::Class()) ? static_cast<MyProfile *>(&hist) : nullptr;

    Double_t * pSumw        = hist.GetArray();
    Double_t * pSumw2       = hist.GetSumw2()->GetArray(); // can be null
    Double_t * pBinEntries  = pProf ? pProf->fBinEntries.GetArray()    : nullptr;
    Double_t * pBinSumw2    = pProf ? pProf->GetBinSumw2()->GetArray() : nullptr;

    const Int_t nSize = hist.GetSize();
    for (Int_t bin = 0; bin < nSize; ++bin)
    {
        pSumw[bin] *= scale;

        if (pSumw2)
            pSumw2[bin] *= scale;

        if (pBinEntries)
            pBinEntries[bin] *= scale;

        if (pBinSumw2)
            pBinSumw2[bin] *= scale;
    }

    hist.ResetStats();

    //LogMsgInfo( "------ after luminosity scale ------" );
    //LogMsgHistStats(hist);
    //LogMsgHistEffectiveEntries(hist);
    //hist.Print("all");

    if (bApplyCrossSectionError)
    {
        double relError = crossSectionError / crossSection;

        for (Int_t bin = 0; bin <= hist.GetNbinsX() + 1; ++bin)
        {
            Double_t binContent = hist.GetBinContent(bin);
            Double_t addError   = binContent * relError;

            Double_t binError   = hist.GetBinError(bin);
            Double_t newError   = std::sqrt( binError * binError + addError * addError );

            hist.SetBinError( bin, newError );
        }

        hist.ResetStats();  // force recalculation of sumw2
    }
}

////////////////////////////////////////////////////////////////////////////////
void ScaleHistToLuminosity( double luminosity, const TH1DVector & hists, size_t nEvents, double crossSection,
                            double crossSectionError, bool bApplyCrossSectionError /*= false*/ )
{
    for (TH1D * pHist : hists)
    {
        if (pHist)
            ScaleHistToLuminosity( luminosity, *pHist, nEvents, crossSection, crossSectionError, bApplyCrossSectionError );
    }
}

////////////////////////////////////////////////////////////////////////////////
TH1D * LoadHist( const char * fileName, const char * histName )
{
    if (gSystem->AccessPathName( fileName ))
        return nullptr;

    struct Cleanup
    {
        TDirectory * oldDir = gDirectory;

        ~Cleanup()
        {
            if (oldDir)
                oldDir->cd();
        }

    } cleanup;

    TFile file( fileName, "READ" );
    if (file.IsZombie() || !file.IsOpen())    // IsZombie is true if constructor failed
        return nullptr;

    // try TProfile first
    {
        TProfile * pHist = nullptr;
        file.GetObject( histName, pHist );
        if (pHist)
        {
            pHist->SetDirectory( nullptr );
            return pHist;
        }
    }

    // try TH1D next
    {
        TH1D * pHist = nullptr;
        file.GetObject( histName, pHist );
        if (pHist)
        {
            pHist->SetDirectory( nullptr );
            return pHist;
        }
    }

    return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
void SaveHists( const char * fileName, const ConstTH1DVector & hists, const char * option )
{
    struct Cleanup
    {
        TDirectory * oldDir = gDirectory;

        ~Cleanup()
        {
            if (oldDir)
                oldDir->cd();
        }

    } cleanup;

    TFile file( fileName, option );
    if (file.IsZombie() || !file.IsOpen())    // IsZombie is true if constructor failed
    {
        LogMsgError( "Failed to create file (%hs).", FMT_HS(fileName) );
        ThrowError( std::invalid_argument( fileName ) );
    }

    for ( const TH1D * pHist : hists )
    {
        if (pHist)
        {
            TH1D * pClone = (TH1D *)pHist->Clone();
            pClone->SetDirectory( &file );  // owned by output file, which will call delete
            pClone->Write( 0, TObject::kOverwrite );
        }
    }

    file.Close();
}

////////////////////////////////////////////////////////////////////////////////
void WriteHists( TFile * pFile, const TH1DVector & hists )
{
    for ( TH1D * pHist : hists )
    {
        if (pHist)
        {
            pHist->SetDirectory( pFile );  // owned by output file, which will call delete
            pHist->Write();
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
TNtupleD * LoadTuple( const char * fileName, const char * tupleName )
{
    if (gSystem->AccessPathName( fileName ))
        return nullptr;

    struct Cleanup
    {
        TDirectory * oldDir = gDirectory;

        ~Cleanup()
        {
            if (oldDir)
                oldDir->cd();
        }

    } cleanup;

    std::unique_ptr<TFile> upFile( new TFile(fileName, "READ") );
    if (upFile->IsZombie() || !upFile->IsOpen())    // IsZombie is true if constructor failed
        return nullptr;

    TNtupleD * pTuple = nullptr;
    upFile->GetObject( tupleName, pTuple );
    if (!pTuple)
        return nullptr;

    TNtupleD * pClone = (TNtupleD *)pTuple->Clone();
    pClone->SetDirectory( nullptr );
    pClone->Reset();
    pClone->ResetBranchAddresses();

    Long64_t nEntries = pTuple->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry)
    {
        pTuple->GetEntry(entry);
        pClone->Fill( pTuple->GetArgs() );
    }

    return pClone;
}

////////////////////////////////////////////////////////////////////////////////
void SaveTuples( const char * fileName, const ConstTupleVector & tuples, const char * option /*= "UPDATE"*/ )
{
    struct Cleanup
    {
        TDirectory * oldDir = gDirectory;

        ~Cleanup()
        {
            if (oldDir)
                oldDir->cd();
        }

    } cleanup;

    TFile file( fileName, option );
    if (file.IsZombie() || !file.IsOpen())    // IsZombie is true if constructor failed
    {
        LogMsgError( "Failed to create file (%hs).", FMT_HS(fileName) );
        ThrowError( std::invalid_argument( fileName ) );
    }

    for ( const TNtupleD * pTuple : tuples )
    {
        if (pTuple)
        {
            TNtupleD * pClone = (TNtupleD *)pTuple->Clone();
            pClone->SetDirectory( &file );  // owned by output file, which will call delete
            pClone->Write( 0, TObject::kOverwrite );
        }
    }

    file.Close();
}

////////////////////////////////////////////////////////////////////////////////
bool LoadCacheTuple( const char * cacheFileName, TNtupleD * & pTuple )
{
    if (!cacheFileName || !cacheFileName[0] || !pTuple)
        return false;

    std::unique_ptr<TNtupleD> pTupleCache( LoadTuple( cacheFileName, pTuple->GetName() ) );
    if (!pTupleCache)
        return false;

    // compare cached object to current object
    {
        if (pTupleCache->Class() != pTuple->Class())
            return false;

        if (pTupleCache->GetNvar() != pTuple->GetNvar())
            return false;

        TObjArray * pBranches      = pTuple     ->GetListOfBranches();
        TObjArray * pBranchesCache = pTupleCache->GetListOfBranches();

        if (pBranchesCache->GetEntries() != pBranches->GetEntries())
            return false;

        Int_t nBranches = pBranches->GetEntries();
        for (Int_t branchIndex = 0; branchIndex < nBranches; ++branchIndex)
        {
            TBranch * pBranch      = (TBranch *)pBranches     ->At(branchIndex);
            TBranch * pBranchCache = (TBranch *)pBranchesCache->At(branchIndex);

            if (strcmp( pBranchCache->GetName(),  pBranch->GetName()  ) != 0)
                return false;
            if (strcmp( pBranchCache->GetTitle(), pBranch->GetTitle() ) != 0)
                return false;

            TObjArray * pLeaves      = pBranch     ->GetListOfLeaves();
            TObjArray * pLeavesCache = pBranchCache->GetListOfLeaves();

            if (pLeavesCache->GetEntries() != pLeaves->GetEntries())
                return false;

            Int_t nLeaves = pLeaves->GetEntries();
            for (Int_t leafIndex = 0; leafIndex < nLeaves; ++leafIndex)
            {
                TLeaf * pLeaf      = (TLeaf *)pLeaves     ->At(leafIndex);
                TLeaf * pLeafCache = (TLeaf *)pLeavesCache->At(leafIndex);

                if (strcmp( pLeafCache->GetName(),  pLeaf->GetName()  ) != 0)
                    return false;
                if (strcmp( pLeafCache->GetTitle(), pLeaf->GetTitle() ) != 0)
                    return false;
            }
        }
    }

    // replace current tuple object with cached one
    delete pTuple;
    pTuple = pTupleCache.release();
    return true;
}

////////////////////////////////////////////////////////////////////////////////
void GetHistDrawMinMax( const TH1D & hist, Double_t & ymin, Double_t & ymax )
{
    ymin = std::numeric_limits<Double_t>::max();
    ymax = -ymin;

    TVirtualPad * oldpad = gPad;

    {
        TCanvas canvas;

        hist.DrawCopy();  // histogram copy is owned by the current pad (i.e. canvas)

        canvas.Update();  // calculate new ranges

        Double_t xmin, xmax;
        canvas.GetRangeAxis( xmin, ymin, xmax, ymax );
    }

    if (oldpad) oldpad->cd();
}

////////////////////////////////////////////////////////////////////////////////
void GetHistDrawMinMax( const ConstTH1DVector & hists, Double_t & ymin, Double_t & ymax )
{
    ymin = std::numeric_limits<Double_t>::max();
    ymax = -ymin;

    for (const TH1D * pHist : hists)
    {
        Double_t hist_ymin, hist_ymax;
        GetHistDrawMinMax( *pHist, hist_ymin, hist_ymax );

        ymin = std::min( ymin, hist_ymin );
        ymax = std::max( ymax, hist_ymax );
    }
}

////////////////////////////////////////////////////////////////////////////////
TH1DVector DrawMultipleHist( const char * title, const ConstTH1DVector & hists, const ColorVector & colors /*= {}*/, const CStringVector drawOptions /*= {}*/ )
{
    TH1DVector drawHists;

    Double_t yAxisMin, yAxisMax;
    GetHistDrawMinMax( hists, yAxisMin, yAxisMax );

    for ( size_t i = 0; i < hists.size(); ++i )
    {
        std::string drawOption = (i < drawOptions.size()) ? drawOptions[i] : "";

        if (i != 0)
            drawOption += " SAME";

        TH1D * pHist = reinterpret_cast<TH1D *>( hists[i]->DrawCopy(drawOption.c_str()) );     // histogram copy is owned by the current pad
        if (!pHist)
            ThrowError( "DrawCopy failed" );

        drawHists.push_back(pHist);

        if (i < colors.size())
        {
            Color_t color = colors[i];
            pHist->SetLineColor(   color );
            pHist->SetMarkerColor( color );
        }

        pHist->SetBit( TH1::kNoTitle );  // disable title from histogram

        // set the y-axis min/max (Note: do not use TCanvas::RangeAxis as this only works if TCanvas::Range is also set appropriately).
        pHist->SetMinimum( yAxisMin );
        pHist->SetMaximum( yAxisMax );
    }

    // add the title, if defined
    if (title && title[0])
    {
        TH1D            dummyHist;
        THistPainter    painter;

        dummyHist.SetDirectory( nullptr );  // ensure not owned by any directory
        dummyHist.SetTitle(title);
        
        painter.SetHistogram( &dummyHist );
        
        painter.PaintTitle();  // creates a TPaveLabel with the name "title" which is owned by the current pad
    }
    
    return drawHists;
}

////////////////////////////////////////////////////////////////////////////////
Double_t GetHistBinEffectiveEntries( const TH1D & hist, Int_t bin )
{
    if (hist.InheritsFrom(TProfile::Class()))
        return static_cast<const TProfile &>(hist).GetBinEffectiveEntries(bin);

    if ((bin < 0) || (bin >= hist.GetSize()))
        return 0;

    Double_t sumW = hist.GetArray()[bin];

    if (bin >= hist.GetSumw2()->GetSize())
        return sumW;

    Double_t sumW2 = hist.GetSumw2()->GetArray()[bin];

    Double_t nEff = (sumW2 > 0 ? sumW * sumW / sumW2 : sumW);

    return nEff;
}

////////////////////////////////////////////////////////////////////////////////
size_t HistNonEmptyBinCount( const TH1D & hist, bool bIncludeUnderOverflow /*= false*/ )
{
    size_t nNonEmpty(0);

    const Int_t first = bIncludeUnderOverflow ? 0 : 1;
    const Int_t last  = hist.GetSize() - 1 - first;

    for (Int_t bin = first; bin <= last; ++bin)
    {
        if (GetHistBinEffectiveEntries( hist, bin ) != 0)
            ++nNonEmpty;
    }

    return nNonEmpty;
}

////////////////////////////////////////////////////////////////////////////////
size_t HistErrorBinCount( const TH1D & hist, bool bIncludeUnderOverflow /*= false*/ )
{
    size_t nError(0);

    const Int_t first = bIncludeUnderOverflow ? 0 : 1;
    const Int_t last  = hist.GetSize() - 1 - first;

    for (Int_t bin = first; bin <= last; ++bin)
    {
        if (hist.GetBinError(bin) != 0)
            ++nError;
    }

    return nError;
}

////////////////////////////////////////////////////////////////////////////////
size_t HistNonEmptyBinCount( const TH1D & h1, const TH1D & h2, bool bCountUnion /*= false*/, bool bIncludeUnderOverflow /*= false*/ )
{
    if (h1.GetSize() != h2.GetSize())
        ThrowError( "HistNonEmptyBinCount: histogram size mismatch." );

    size_t nNonEmpty(0);

    const Int_t first = bIncludeUnderOverflow ? 0 : 1;
    const Int_t last  = h1.GetSize() - 1 - first;

    for (Int_t bin = first; bin <= last; ++bin)
    {
        bool ne1 = (GetHistBinEffectiveEntries( h1, bin ) != 0);
        bool ne2 = (GetHistBinEffectiveEntries( h2, bin ) != 0);

        if ((ne1 && ne2) || (bCountUnion && (ne1 || ne2)))
            ++nNonEmpty;
    }

    return nNonEmpty;
}

////////////////////////////////////////////////////////////////////////////////
size_t HistErrorBinCount( const TH1D & h1, const TH1D & h2, bool bCountUnion /*= false*/, bool bIncludeUnderOverflow /*= false*/ )
{
    if (h1.GetSize() != h2.GetSize())
        ThrowError( "HistErrorBinCount: histogram size mismatch." );

    size_t nError(0);

    const Int_t first = bIncludeUnderOverflow ? 0 : 1;
    const Int_t last  = h1.GetSize() - 1 - first;

    for (Int_t bin = first; bin <= last; ++bin)
    {
        bool err1 = (h1.GetBinError(bin) != 0);
        bool err2 = (h2.GetBinError(bin) != 0);

        if ((err1 && err2) || (bCountUnion && (err1 || err2)))
            ++nError;
    }

    return nError;
}

////////////////////////////////////////////////////////////////////////////////
void ZeroHistBin( TH1D & hist, Int_t bin )
{
    hist.SetBinContent( bin, 0 );
    hist.SetBinError(   bin, 0 );

    if (hist.InheritsFrom(TProfile::Class()))
        static_cast<TProfile &>(hist).SetBinEntries( bin, 0 );
}

////////////////////////////////////////////////////////////////////////////////
void ZeroHistEmptyBins( TH1D & h1, TH1D & h2 )
{
    // zero bins if either are zero

    if (h1.GetSize() != h2.GetSize())
        ThrowError( "ZeroHistEmptyBins: histogram size mismatch." );

    const Int_t nSize = h1.GetSize();
    for (Int_t bin = 0; bin < nSize; ++bin)  // include under/overflow bins
    {
        bool empty1 = (GetHistBinEffectiveEntries( h1, bin ) == 0);
        bool empty2 = (GetHistBinEffectiveEntries( h2, bin ) == 0);

        if (empty1 == empty2)  // both true or both false
            continue;

        ZeroHistBin( empty1 ? h2 : h1, bin );
    }
}

////////////////////////////////////////////////////////////////////////////////
Double_t KolmogorovTest_NonEmptyBins( const TH1D & h1, const TH1D & h2 )
{
    TH1DUniquePtr p1( (TH1D *)h1.Clone() );
    TH1DUniquePtr p2( (TH1D *)h2.Clone() );

    ZeroHistEmptyBins( *p1, *p2 );  // zero bins if either are zero

    // KolmogorovTest ignores bins where BOTH are zero
    return p1->KolmogorovTest( p2.get() );
}

////////////////////////////////////////////////////////////////////////////////
Double_t HistPointChi2Test( const TH1D & p1, const TH1D & p2, Double_t & chi2, Int_t & ndf )
{
    // Use this with TProfile, or as an alternative chi2test to those implemented in TH1D
    // Calculate chi2 as if the histograms are just measurement plots.

    if (p1.GetSize() != p2.GetSize())
        ThrowError( "HistValueChi2Test: profile size mismatch." );

    const Int_t nBins = p1.GetSize() - 2;

    chi2 = 0;
    ndf  = nBins - 1;

    for (Int_t bin = 1; bin <= nBins; ++bin)  // do not include under/overflow
    {
        Double_t n1 = GetHistBinEffectiveEntries(p1, bin);
        Double_t n2 = GetHistBinEffectiveEntries(p2, bin);

        if ((n1 == 0) || (n2 == 0))
        {
            // skip this bin, and reduce the ndf
            --ndf;
            continue;
        }

        Double_t v1 = p1.GetBinContent(bin);
        Double_t v2 = p2.GetBinContent(bin);

        Double_t e1 = p1.GetBinError(bin);
        Double_t e2 = p2.GetBinError(bin);

        Double_t delta  = v1 - v2;
        Double_t errsqr = e1*e1 + e2*e2;

        if (errsqr == 0)
        {
            // skip this bin, and reduce the ndf
            --ndf;
            continue;
        }

        Double_t chisqr = delta * delta / errsqr;

        chi2 += chisqr;
    }

    if (ndf < 0) ndf = 0; 

    Double_t prob = TMath::Prob( chi2, ndf );   // can handle ndf <= 0 (see TMath.cxx)

    return prob;
}

////////////////////////////////////////////////////////////////////////////////
void Chi2Result::Chi2Test( const TH1D & h1, const TH1D & h2 )
{
    if (!h1.InheritsFrom(TProfile::Class()) && !h2.InheritsFrom(TProfile::Class()))
    {
        // both are TH1D

        // We want the Chi2Test to skip bins if either are empty,
        // however, Chi2TestX only skips bins if both are empty.
        // To accomplish the desired behavior we make sure that if a bin is empty
        // in one, it is also empty in the other.

        TH1DUniquePtr p1( (TH1D *)h1.Clone() );
        TH1DUniquePtr p2( (TH1D *)h2.Clone() );

        ZeroHistEmptyBins( *p1, *p2 );  // zero bins if either are zero

        LogMsgInfo( "Chi2Test(%hs, %hs): %u -> %u non-empty bins", FMT_HS(h1.GetName()), FMT_HS(h2.GetName()),
                    FMT_U(HistNonEmptyBinCount(h1,h2,true)), FMT_U(HistNonEmptyBinCount(*p1,*p2,true)) );

        // perform chi2test

        prob     = p1->Chi2TestX( p2.get(), chi2, ndf, igood, "WW" );
        chi2_ndf = (ndf > 0 ? chi2 / ndf : 0.0);
    }
    else
    {
        if (!h1.InheritsFrom(TProfile::Class()) || !h2.InheritsFrom(TProfile::Class()))
            ThrowError( "Chi2Test: both histograms must inherit from TProfile" );

        // both are TProfile

        const TProfile & p1 = static_cast<const TProfile &>(h1);
        const TProfile & p2 = static_cast<const TProfile &>(h2);

        prob     = HistPointChi2Test( p1, p2, chi2, ndf );
        chi2_ndf = (ndf > 0 ? chi2 / ndf : 0.0);
        igood    = 0;
    }
}

////////////////////////////////////////////////////////////////////////////////
std::string Chi2Result::Label()
{
    char label[200];
    sprintf( label, "#chi^{2}/ndf = %.4g/%i = %.4g  p-value = %.4g", FMT_F(chi2), FMT_I(ndf), FMT_F(chi2_ndf), FMT_F(prob) );
    return label;
}

////////////////////////////////////////////////////////////////////////////////

}  // namespace RootUtil
