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
#include <TFile.h>
#include <TH1.h>
#include <TProfile.h>
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
void LoadEvents( const char * eventFileName, std::function<void(const HepMC::GenVertex & signal)> EventFunc )
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

    HepMC::GenEvent genEvent;

    while (upInput->fill_next_event( &genEvent ))
    {
        const HepMC::GenVertex * pSignal = genEvent.signal_process_vertex();
        if (!pSignal)
            ThrowError( "Missing signal vertex for event." );

        EventFunc( *pSignal );
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

    Double_t stats[TH1::kNstat];
    hist.GetStats( stats );

    LogMsgInfo( "%hs: sumw=%g sumw2=%g sumwx=%g sumwx2=%g sumwy=%g sumwy2=%g", FMT_HS(hist.GetName()),
                FMT_F(stats[0]), FMT_F(stats[1]), FMT_F(stats[2]), FMT_F(stats[3]),     // TH1D and TProfile
                FMT_F(stats[4]), FMT_F(stats[5]) );                                     // TProfile
}

////////////////////////////////////////////////////////////////////////////////
TH1D * ConvertTProfileToTH1D( const TH1D * pProfile, bool bDeleteProfile )
{
    TH1D * pHist = nullptr;

    if (pProfile->InheritsFrom(TProfile::Class()))
    {
        // create a new TH1D from the TProfile
        pHist = ((TProfile *)pProfile)->ProjectionX("");  // do not add a suffix

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
void SetupHist( TH1D & hist, const char * xAxisTitle, const char * yAxisTitle,
                Color_t lineColor /*= -1*/, Color_t markerColor /*= -1*/, Color_t fillColor /*= -1*/ )
{
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
void WriteHists( TFile * pFile, const TH1DVector & hists )
{
    for ( TH1D * pHist : hists )
    {
        pHist->SetDirectory( pFile );  // owned by output file, which will call delete
        pHist->Write();
    }
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

}