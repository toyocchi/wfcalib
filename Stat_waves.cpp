// Simulation of Statistical Method by Waveforms
//--/----|----/----|----/----|----/----|----/----|----/----|----/----|----/----|
#include "Waveform_stat.h"
#include <random>

void Init(Double_t tauTiming, Double_t tauDecay, Double_t trueGain);

void Stat_waves() {
	const Int_t N = 256;
	Double_t canvX = 1150., canvY = 650;
	Int_t nBin = 1000;
	Double_t histMin = -.08;
	Double_t histMax = .6;
//	Double_t histMax = .1;
	Double_t samplingRate = 1.6;

	std::random_device seed;
	std::mt19937 rnd(seed());

	Double_t analysisStart = -230., analysisEnd = -80.;

	Int_t nPnt = samplingRate * (analysisEnd - analysisStart);

	Double_t tauTiming = .001;
	Double_t tauDecay = 30;
	Double_t trueGain = 0.04534;

/*	Double_t lambda = 0.01;
	Double_t alpha = 0.134;
	Double_t dcr = 500; // kcps
	Double_t noise = .0005;
	Double_t noise1 = .0003;
*/
	Double_t lambda = .03;
	Double_t alpha = .2;
	Double_t alphaCT = .2;
	Double_t dcr = 500;
	Double_t noise = .00085;
	Double_t noise1 = .0004;

	Int_t nEventWeak = 40;
	Int_t nEventStrong = 40000;

	Int_t nMu = 10;
	
   Init(tauTiming, tauDecay, trueGain);

	// weak
	Double_t muWeak = 2.74;
	std::poisson_distribution<> poisson(muWeak);

	TH1F* histWeak = new TH1F("histWeak", "weak", nBin, histMin, histMax);
	TH1F* histPed = new TH1F("histPed", "Pedestal", 300, histMin, trueGain);

	cout << "Processing weak run..." << endl;
	for(int iEvent=0; iEvent<nEventWeak; iEvent++){
		Waveform* wf = 
			new Waveform(nPnt, analysisStart, analysisEnd, lambda, alpha, alphaCT,  dcr);
		Int_t nPhoton = poisson(rnd);
		Int_t nHit = wf->MakeEvent(nPhoton);
		wf->SetNoiseAmp(noise+noise1*nHit);
		Double_t charge = wf->GetCharge();
/*		if(charge>.045 && charge<.055) {
			cout << funcSingle->GetParameter(1) << endl;
		}
*/
		histWeak->Fill(charge);
		histPed->Fill(charge);
		delete wf;
	}

	TCanvas* canvHistWeak = 
		new TCanvas("canvHistWeak", "weak hist", canvX, canvY);
//	canvHistWeak->SetLogy();
	histWeak->Draw();

	TF1* gauss = new TF1("gauss", "[0]*exp(-(x-[1])^2/2/[2]^2)");
	gauss->SetParameters(40, 0, .001);
	gauss->SetParLimits(1, -.01, .01);
	TCanvas* canvPed = new TCanvas("canvPed", "Pedestal", canvX, canvY);
	histPed->Fit(gauss, "q", "", -.1, trueGain/2);

	Double_t mean = histWeak->GetMean();
	Double_t meanErr = histWeak->GetMean(11);
	Double_t stdDev = histWeak->GetStdDev();
	Double_t stdDevErr = histWeak->GetStdDev(11);
	Double_t var = stdDev*stdDev;
	Double_t varErr = 2*stdDev*stdDevErr;
	
	Double_t par0 = gauss->GetParameter(0);
	Double_t par0Err = gauss->GetParError(0);
	Double_t par2 = gauss->GetParameter(2);
	Double_t par2Err = gauss->GetParError(2);

	Double_t area = par0*sqrt(2*TMath::Pi())*abs(par2);
	Double_t areaErr = area*sqrt(pow(par0Err/par0,2)
										  + pow(par2Err/par2,2));
	Double_t f0 = area/nEventWeak*nBin/(trueGain-histMin);
	Double_t f0Err = sqrt(f0*(1-f0)/nEventWeak + pow(areaErr/nEventWeak,2));
	Double_t mu = -log(f0);
	Double_t muErr = f0Err/f0;

	Double_t eqf = mean/trueGain/mu;
	Double_t eqfErr = sqrt(pow(meanErr/mu,2)+pow(mean/mu/mu*muErr,2));
	Double_t enf = mu*var/mean/mean;
	Double_t enfErr = sqrt(pow(enf*muErr/mu,2) + pow(enf*varErr/var,2)
								  + pow(2*mu*var/mean/mean/mean*meanErr,2));

	cout << "par0: " << par0 << ", par2: " << par2 << ", area: " << area << endl;
	cout << "mean: " << mean << ", mu: " << mu << endl;
	cout << "f0: " << f0 << ", eqf: " << eqf << ", enf: " << enf << endl;

	// strong
	TGraphErrors* gr = new TGraphErrors;
	TClonesArray* cahStrong = new TClonesArray("TH1F", nMu);
	TCanvas* canvHistStrong = new TCanvas("canvHistStrong", "strong hist",
													  canvX, canvY);
	canvHistStrong->Divide(4,3);

	for(Int_t iMu=0; iMu<nMu; iMu++) {
		Double_t mu = 15 + 3 * (iMu+1);
		std::poisson_distribution<> poisson(mu);
		Char_t histName[N];
		sprintf(histName, "hist (mu=%2.f)", mu);
		new ((TH1F*)(*cahStrong)[iMu]) TH1F(histName, histName, nBin, -.1, 15);

		cout << "Processing run mu = " << mu << "..." << endl;
		for(int iEvent=
				 0; iEvent<nEventStrong; iEvent++){
			Waveform* wf = 
				new Waveform(nPnt, analysisStart, analysisEnd, lambda, alpha, alphaCT, dcr);
			Int_t nPhoton = poisson(rnd);
			wf->MakeEvent(nPhoton);
			wf->SetNoiseAmp(noise);
			Double_t charge = wf->GetCharge();
			((TH1F*)(*cahStrong)[iMu])->Fill(charge);
			delete wf;
		}
		
		canvHistStrong->cd(iMu+1);
		((TH1F*)(*cahStrong)[iMu])->Draw();
		
		Double_t mean = ((TH1F*)(*cahStrong)[iMu])->GetMean();
		Double_t meanErr = ((TH1F*)(*cahStrong)[iMu])->GetMean(11);
		Double_t stdDev = ((TH1F*)(*cahStrong)[iMu])->GetStdDev();
		Double_t stdDevErr = ((TH1F*)(*cahStrong)[iMu])->GetStdDev(11);
		Double_t var = stdDev * stdDev;
		Double_t varErr = 2*stdDev*stdDevErr;
		
		gr->SetPoint(iMu, mean, var);
		gr->SetPointError(iMu, meanErr, varErr);
	}
	
	TF1* linear = new TF1("linear", "pol1");
	linear->SetLineColor(kRed);
	TCanvas* canvGr = new TCanvas("canvGr", "var vs mean", canvX, canvY);
	gr->SetTitle("var vs mean;mean [10^10e];var [(10^10e)^2]");
	gr->SetMinimum(0.);
	gr->GetXaxis()->SetLimits(0., 2.5);
	gr->SetMarkerStyle(8);
	gr->SetMarkerColor(1);
	gr->Draw("ap");
	gr->Fit(linear, "Q", "same");

	Double_t gain = linear->GetParameter(1);
	Double_t gainErr = linear->GetParError(1);
	Double_t gainCorr = gain/eqf/enf;
	Double_t gainCorrErr = sqrt(pow(gainErr/eqf/enf,2) 
										 + pow(gain/eqf/eqf/enf*eqfErr,2)
										 + pow(gain/eqf/enf/enf*enfErr,2));

	cout << "true gain: " << trueGain <<
		", gain: " << gain << " +- " << gainErr <<
		", corrected gain: " << gainCorr << " +- " << gainCorrErr << endl;

}

void Init(Double_t tauTiming, Double_t tauDecay, Double_t trueGain) {
	// sctime: PDF of pulse timing (exponential decay)
   funcExpPDF->SetParameter(0,tauTiming);
	// single: 0:decay const, 1:start time, 2:single charge
   funcSingle->FixParameter(0,tauDecay);
   funcSingle->FixParameter(2,trueGain);

}
