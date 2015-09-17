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
void FillHistPT( const HepMC::GenVertex & signal, TH1D & hist, double weight, int pdg )
{
    auto itrPart = signal.particles_out_const_begin();
    auto endPart = signal.particles_out_const_end();
    for ( ; itrPart != endPart; ++itrPart)
    {
        const HepMC::GenParticle & part = **itrPart;

        if (part.pdg_id() == pdg)
        {
            TLorentzVector vec = ToLorentz( part.momentum() );

            double pT = vec.Pt();
            hist.Fill( pT, weight );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void FillHistRap( const HepMC::GenVertex & signal, TH1D & hist, double weight, int pdg )
{
    auto itrPart = signal.particles_out_const_begin();
    auto endPart = signal.particles_out_const_end();
    for ( ; itrPart != endPart; ++itrPart)
    {
        const HepMC::GenParticle & part = **itrPart;

        if (part.pdg_id() == pdg)
        {
            TLorentzVector vec = ToLorentz( part.momentum() );

            double eta = vec.Rapidity();
            hist.Fill( eta, weight );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void FillHistEta( const HepMC::GenVertex & signal, TH1D & hist, double weight, int pdg )
{
    auto itrPart = signal.particles_out_const_begin();
    auto endPart = signal.particles_out_const_end();
    for ( ; itrPart != endPart; ++itrPart)
    {
        const HepMC::GenParticle & part = **itrPart;

        if (part.pdg_id() == pdg)
        {
            TLorentzVector vec = ToLorentz( part.momentum() );

            double eta = vec.Eta();
            hist.Fill( eta, weight );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void FillHistPhi( const HepMC::GenVertex & signal, TH1D & hist, double weight, int pdg )
{
    auto itrPart = signal.particles_out_const_begin();
    auto endPart = signal.particles_out_const_end();
    for ( ; itrPart != endPart; ++itrPart)
    {
        const HepMC::GenParticle & part = **itrPart;

        if (part.pdg_id() == pdg)
        {
            TLorentzVector vec = ToLorentz( part.momentum() );

            double phi = vec.Phi();
            hist.Fill( phi, weight );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void FillHistMass( const HepMC::GenVertex & signal, TH1D & hist, double weight, int pdg1, int pdg2 )
{
    const HepMC::GenParticle * pPart1 = nullptr;
    const HepMC::GenParticle * pPart2 = nullptr;

    auto itrPart = signal.particles_out_const_begin();
    auto endPart = signal.particles_out_const_end();
    for ( ; itrPart != endPart; ++itrPart)
    {
        const HepMC::GenParticle * pPart = *itrPart;

        int part_pdg = pPart->pdg_id();

        if (part_pdg == pdg1)
        {
            pPart1 = pPart;
            continue;
        }
        if (part_pdg == pdg2)
        {
            pPart2 = pPart;
            continue;
        }
    }

    if (pPart1 && pPart2)
    {
        TLorentzVector vec1 = ToLorentz( pPart1->momentum() );
        TLorentzVector vec2 = ToLorentz( pPart2->momentum() );
        TLorentzVector vec  = vec1 + vec2;

        double mass = vec.M();
        hist.Fill( mass, weight );
    }
}

////////////////////////////////////////////////////////////////////////////////
void LogMsgUnderOverflow( const TH1D & hist )
{
    Double_t underflow = hist.GetBinContent(0);
    Double_t overflow  = hist.GetBinContent( hist.GetNbinsX() + 1 );

    if (underflow != 0)
        LogMsgInfo( "Underflow in %hs of %g", FMT_HS(hist.GetName()), FMT_F(underflow) );

    if (overflow != 0)
        LogMsgInfo( "Overflow in %hs of %g", FMT_HS(hist.GetName()), FMT_F(overflow) );
}

////////////////////////////////////////////////////////////////////////////////
void LogMsgUnderOverflow( const ConstTH1DVector & hists )
{
    for (const TH1D * pHist : hists)
    {
        if (pHist)
            LogMsgUnderOverflow( *pHist );
    }
}

////////////////////////////////////////////////////////////////////////////////
void SetupHist( TH1D & hist, const char * xAxisTitle, const char * yAxisTitle,
                Color_t lineColor /*= -1*/, Color_t markerColor /*= -1*/, Color_t fillColor /*= -1*/ )
{
    hist.Sumw2();
    hist.SetStats( kFALSE );

    hist.GetXaxis()->SetTitle( xAxisTitle );
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