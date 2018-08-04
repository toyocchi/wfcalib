//--/----|----/----|----/----|----/----|----/----|----/----|----/----|----/----|
#include "Waveform_WFCalib.h"
#include <random>

void Init(Double_t tauTiming, Double_t tauDecay, Double_t trueGain);

void WFCalib_sample() {
	Double_t analysisStart = -230., analysisEnd = -80.;
	Double_t samplingRate = 1.6;
	Int_t defPnt = samplingRate * (analysisEnd - analysisStart);

	Double_t tauTiming = 45.;
	Double_t tauDecay = 30.;
	Double_t trueGain=0.04534;

	Int_t nEvent = 40000;
	Double_t noise = 0.;

	Int_t nSample = 5, nMu = 5;

	const Int_t N = 256;
	Double_t canvX = 1150., canvY = 650;
	Int_t factor = 1;

	// prepare for random number
	std::random_device seed;
	std::mt19937 rnd(seed());
	
	Init(tauTiming, tauDecay, trueGain);
	
	TClonesArray* cagrGainNoiseCharge =
		new TClonesArray("TGraphErrors", nSample);
	TClonesArray* cagrGainDiffVar =
		new TClonesArray("TGraphErrors", nSample);
	
	TMultiGraph* mgGainNoiseCharge = new TMultiGraph;
	TMultiGraph* mgGainDiffVar = new TMultiGraph;
	TLegend* legGainNoiseCharge = new TLegend(0.1, 0.7, 0.48, 0.9);
	TLegend* legGainDiffVar = new TLegend(0.1, 0.7, 0.48, 0.9);
	
	TCanvas* canvWF = new TCanvas("canvWF", "Waveform", canvX, canvY);
	canvWF->Divide(2,1);
	TCanvas* canvNoiseCharge = 
		new TCanvas("canvNoiseCharge", "noise charge", canvX, canvY);
	canvNoiseCharge->Divide(3,2);
	TCanvas* canvDiffVar = 
		new TCanvas("canvDiffVar", "diffvar", canvX, canvY);
	canvDiffVar->Divide(3,2);
		
	for (Int_t iSample = 0; iSample < nSample; iSample++) {
		Int_t nPnt = defPnt*(iSample+1)*factor;
		cout << "sampling rate: " << samplingRate*(iSample+1)*factor 
			  << " GHz" << endl;
		new ((TGraphErrors*)(*cagrGainNoiseCharge)[iSample]) TGraphErrors();
		new ((TGraphErrors*)(*cagrGainDiffVar)[iSample]) TGraphErrors();
		
		TClonesArray* cagrNoiseCharge =
			new TClonesArray("TGraphErrors", nSample);
		TClonesArray* cagrDiffVar =
			new TClonesArray("TGraphErrors", nSample);
		
		TMultiGraph* mgNoiseCharge = new TMultiGraph;
		TMultiGraph* mgDiffVar = new TMultiGraph;
		Double_t gainNoiseChargeBasis, gainDiffVarBasis;

		for(Int_t iMu=0; iMu<nMu; iMu++) {
			Int_t mu = 10 * (iMu+1);
			cout << "mu: " << mu << endl;
			std::poisson_distribution<> poisson(mu);
			
			new ((TGraphErrors*)(*cagrNoiseCharge)[iMu]) TGraphErrors();
			new ((TGraphErrors*)(*cagrDiffVar)[iMu]) TGraphErrors();
				
			new ((TGraphErrors*)(*cagrNoiseCharge)[1]) TGraphErrors();
			new ((TGraphErrors*)(*cagrDiffVar)[1]) TGraphErrors();
			
			for(int iEvent=0; iEvent<nEvent; iEvent++){
				Waveform* wf = 
					new Waveform(nPnt, analysisStart, analysisEnd);
				Int_t nPhoton = poisson(rnd);
				wf->MakeEvent(nPhoton);
				wf->SetNoiseAmp(noise);
					
				if((iSample==nSample-1 && iMu==nMu-1) && iEvent==nEvent-1) {
					canvWF->cd(1);
					TGraph* grNoDiff = new TGraph;
					for(Int_t iPnt=0; iPnt<nPnt; iPnt++) {
						Double_t time = wf->GetTimeAt(iPnt);
						Double_t ampSignal = wf->GetSignalAmpAt(iPnt);
						Double_t ampNoise = wf->GetNoiseAmpAt(iPnt);
						Double_t ampTotal = ampSignal + ampNoise;
						grNoDiff->SetPoint(iPnt, time, ampTotal);
					}
					grNoDiff
						->SetTitle("Waveform before Diffentiation; time [ns]; Amplitude [mV]");
					
					grNoDiff->Draw("al");
				}
				
				Double_t charge = wf->GetCharge();
				if(charge<.5*trueGain) {
					continue;
				}
				
				wf->Differentiate(2);
				
				if((iSample==nSample-1 && iMu==nMu-1) && iEvent==nEvent-1) {
					canvWF->cd(2);
					TGraph* grTwoDiff = new TGraph;
					for(Int_t iPnt=0; iPnt<nPnt; iPnt++) {
						Double_t time = wf->GetTimeAt(iPnt);
						Double_t ampSignal = wf->GetSignalAmpAt(iPnt);
						Double_t ampNoise = wf->GetNoiseAmpAt(iPnt);
						Double_t ampTotal = ampSignal + ampNoise;
						grTwoDiff->SetPoint(iPnt, time, ampTotal);
					}
					grTwoDiff->SetTitle("Second Derivative of Waveform; time [ns]; Amplitude [mV]");
					grTwoDiff->Draw("al");
				}

				Double_t diffVar = wf->GetTotalVariance();
				Double_t noiseCharge = charge*charge/diffVar;
				
				((TGraphErrors*)(*cagrNoiseCharge)[iMu])
					->SetPoint(iEvent, noiseCharge, charge);
				((TGraphErrors*)(*cagrNoiseCharge)[iMu])->
					SetPointError(iEvent, 0., 0.);
				((TGraphErrors*)(*cagrDiffVar)[iMu])->
					SetPoint(iEvent, charge, diffVar);
				((TGraphErrors*)(*cagrDiffVar)[iMu])
					->SetPointError(iEvent, 0., 0.);
				delete wf;
			}
				
			TF1* funcLinear = new TF1("funcLinear", "pol1");
			funcLinear->SetLineColor(iMu+2);
			TF1* funcQuad = new TF1("funcQuad", "pol2");
			funcQuad->SetLineColor(iMu+2);
			
			Double_t parNoiseCharge[2];
			Double_t errNoiseCharge[2];
			Double_t parDiffVar[3];
			Double_t errDiffVar[3];
			((TGraphErrors*)(*cagrNoiseCharge)[iMu])->Fit(funcLinear,"Q");
			((TGraphErrors*)(*cagrDiffVar)[iMu])->Fit(funcQuad,"Q");
			for (int iPar = 0; iPar < 3; iPar++) {
				if(iPar!=2) {
					parNoiseCharge[iPar]=funcLinear->GetParameter(iPar);
					errNoiseCharge[iPar]=funcLinear->GetParError(iPar);
					}
				parDiffVar[iPar]=funcQuad->GetParameter(iPar);
				errDiffVar[iPar]=funcQuad->GetParError(iPar);
			}
			Double_t gainNoiseCharge = parNoiseCharge[1];
			Double_t gainNoiseChargeErr = errNoiseCharge[1];
			Double_t gainDiffVar = parDiffVar[1]/(1-parDiffVar[2]);
			Double_t gainDiffVarErr = 
				TMath::Sqrt(TMath::Power(errDiffVar[1]/(1-parDiffVar[2]), 2)
								+ TMath::Power(parDiffVar[1]*errDiffVar[2]
													/(1-parDiffVar[2])/(1-parDiffVar[2])
													, 2));

			if(iMu==0) {
				gainNoiseChargeBasis = gainNoiseCharge;
				gainDiffVarBasis = gainDiffVar;
			}

			((TGraphErrors*)(*cagrGainNoiseCharge)[iSample])->
				SetPoint(iMu, mu, gainNoiseCharge/gainNoiseChargeBasis);
			((TGraphErrors*)(*cagrGainNoiseCharge)[iSample])->
				SetPointError(iMu, 0., gainNoiseChargeErr/gainNoiseChargeBasis);
			((TGraphErrors*)(*cagrGainDiffVar)[iSample])->
				SetPoint(iMu, mu, gainDiffVar/gainDiffVarBasis);
			((TGraphErrors*)(*cagrGainDiffVar)[iSample])->
				SetPointError(iMu, 0., gainDiffVarErr/gainDiffVarBasis);
			
			((TGraphErrors*)(*cagrDiffVar)[iMu])->SetMinimum(0);
			((TGraphErrors*)(*cagrNoiseCharge)[iMu])->SetMinimum(0);
			((TGraphErrors*)(*cagrDiffVar)[iMu])->SetMarkerStyle(8);
			((TGraphErrors*)(*cagrNoiseCharge)[iMu])->SetMarkerStyle(8);
			((TGraphErrors*)(*cagrDiffVar)[iMu])
				->SetMarkerColor(2+iMu);
			((TGraphErrors*)(*cagrNoiseCharge)[iMu])
				->SetMarkerColor(2+iMu);
			mgNoiseCharge->Add(((TGraphErrors*)(*cagrNoiseCharge)[iMu]));
			mgDiffVar->Add(((TGraphErrors*)(*cagrDiffVar)[iMu]));
		}
		
		Char_t nameNoiseCharge[N];
		sprintf(nameNoiseCharge, 
				  "Charge vs Noise Charge (%.1f GHz);Noise Charge [s^4];Charge [10^10e]", 
				  samplingRate*(iSample+1)*factor);
		Char_t nameDiffVar[N];
		sprintf(nameDiffVar, 
				  "DiffVar vs Charge (%.1f GHz);Charge [10^10e];DiffVar[(10^10e/s/s)^2]",
					  samplingRate*(iSample+1)*factor);
		
		mgNoiseCharge->SetTitle(nameNoiseCharge);
		mgDiffVar->SetTitle(nameDiffVar);
		canvNoiseCharge->cd(iSample+1);
		mgNoiseCharge->Draw("ap");
		Double_t xMaxNoiseCharge = mgNoiseCharge->GetXaxis()->GetXmax();
		mgNoiseCharge->GetXaxis()->SetLimits(0., xMaxNoiseCharge);
		canvDiffVar->cd(iSample+1);
		mgDiffVar->Draw("ap");
		Double_t xMaxDiffVar = mgDiffVar->GetXaxis()->GetXmax();
		mgNoiseCharge->GetXaxis()->SetLimits(0., xMaxDiffVar);
		
		((TGraphErrors*)(*cagrGainNoiseCharge)[iSample])->SetMarkerStyle(8);
		((TGraphErrors*)(*cagrGainNoiseCharge)[iSample])->
			SetMarkerColor(iSample+2);
		((TGraphErrors*)(*cagrGainDiffVar)[iSample])->SetMarkerStyle(8);
		((TGraphErrors*)(*cagrGainDiffVar)[iSample])->
			SetMarkerColor(iSample+2);
		
		mgGainNoiseCharge->Add(((TGraphErrors*)(*cagrGainNoiseCharge)[iSample]));
		mgGainDiffVar->Add(((TGraphErrors*)(*cagrGainDiffVar)[iSample]));
		Char_t name[N];
		sprintf(name, "%.1f Hz", samplingRate*(iSample+1)*factor);
		legGainNoiseCharge->
			AddEntry(((TGraphErrors*)(*cagrGainNoiseCharge)[iSample])
						, name, "pe");
		legGainDiffVar->
			AddEntry(((TGraphErrors*)(*cagrGainDiffVar)[iSample])
						, name, "pe");
	}
	
	mgGainNoiseCharge->SetTitle("gain from noise charge vs mu;mu [photon];gain relative to smallest mu");
	mgGainDiffVar->SetTitle("gain from diffvar vs mu;mu [photon];gain relative to smallest mu");
	TCanvas* canvGainNoiseCharge =
		new TCanvas("canvGainNoiseCharge", "gain from noise charge", 
						canvX, canvY);
	mgGainNoiseCharge->Draw("ap");
	legGainNoiseCharge->Draw();
	TCanvas* canvGainDiffVar =
		new TCanvas("canvGainDiffVar", "gain from diffvar", canvX, canvY);
	mgGainDiffVar->Draw("ap");
	legGainDiffVar->Draw("ap");
}
	

	
void Init(Double_t tauTiming, Double_t tauDecay, Double_t trueGain) {
	// sctime: PDF of pulse timing (exponential decay)
   funcTimingPDF->FixParameter(0,tauTiming);
	// single: 0:decay const, 1:start time, 2:single charge
   funcSingle->FixParameter(0,tauDecay);
   funcSingle->FixParameter(2,trueGain);
}
