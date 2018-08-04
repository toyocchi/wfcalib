// Simulation of Statistical Method by charge spectrum
//--/----|----/----|----/----|----/----|----/----|----/----|----/----|----/----|
#include "myphFunc.h"

Double_t ExpectedSpectrum(Double_t *x,Double_t *par);

void Stat_func() {
   const Int_t N = 256;
	Double_t canvX = 1150, canvY = 650;
	Double_t min = -.5, max = .6;

	const Int_t nPar = 10;
//	Double_t parIni[nPar] = {.137, 2.74, .01, .00507, .00271, .0252,
	Double_t parIni[nPar] = {.137, 2.74, .01, .00507, .05271, .0252,
									 .04534, -.00065, 18.88};
	
	Int_t nMu = 10;

	TF1 *gpap 
		= new TF1("gpap", ExpectedSpectrum, min, max, nPar);
	gpap->SetTitle("gpap;charge;# of Events");
	gpap->SetLineColor(kRed);

	gpap->SetParNames("alpha", "mu", "lambda", "sigma0", "sigma1", "beta", 
                        "gain", "ped", "norm");

	for(Int_t iPar=0; iPar<nPar; iPar++) {
		gpap->SetParameter(iPar, parIni[iPar]);
	}

//	gpap->SetParameter(5, 10*parIni[6]);

	gpap->FixParameter(9, 15); 

	const Int_t nRandom = 40000;
	Double_t rand[nRandom];
	
	Double_t mean = 0;
	for(Int_t i=0; i<nRandom; i++) {
		rand[i] = gpap->GetRandom();
		mean += rand[i];
	}
	mean /= nRandom;

	Double_t var = 0;
	for(Int_t i=0; i<nRandom; i++) {
		var += pow(rand[i]-mean,2);
	}
	var /= nRandom;

	Double_t meanErr = sqrt(var/nRandom);
	Double_t varErr = sqrt(var/2/nRandom);

	mean -= parIni[7];
	mean /= parIni[6];
	var -= pow(parIni[3],2);
	var /= pow(parIni[6],2);

	Double_t mu = gpap->GetParameter(1);
	
	Double_t eqf = mean/mu;
	Double_t eqfErr = meanErr/mu;
	Double_t enf = mu*var/mean/mean;
	Double_t enfErr = sqrt(pow(mu*varErr/mean/mean,2) +
								  pow(2*mu*var/pow(mean,3)*meanErr,2));
	cout << "eqf: " << eqf << ", enf: " << enf << endl;

	TCanvas* canvWeak = new TCanvas("canvHist", "Weak", canvX, canvY);
	gpap->SetNpx(1000);
	gpap->Draw();
	

// strong
	TF1 *gpapStrong 
		= new TF1("gpapStrong", ExpectedSpectrum, min, 3, nPar);
	gpapStrong->SetNpx(1000);
	for(Int_t iPar=0; iPar<nPar; iPar++) {
		gpapStrong->SetParameter(iPar, parIni[iPar]);
	}
	TGraphErrors* gr = new TGraphErrors;
	TCanvas* canvStrong = new TCanvas("canvStrong", "Strong", canvX, canvY);
	canvStrong->Divide(4,3);
		
	for(Int_t iMu=0; iMu<nMu; iMu++) {
		Double_t mu = 3*(iMu+1);
		cout << "mu: " << mu << endl;
		gpapStrong->SetParameter(1, mu);
		gpapStrong->FixParameter(9, 15+5*iMu);
		canvStrong->cd(iMu+1);
		gStyle->SetOptFit();
		gpapStrong->Draw();

		Double_t mean = 0;
		for(Int_t i=0; i<nRandom; i++) {
			rand[i] = gpapStrong->GetRandom();
			mean += rand[i];
		}
		mean /= nRandom;
		
		Double_t var = 0;
		for(Int_t i=0; i<nRandom; i++) {
			var += pow(rand[i]-mean,2);
		}
		var /= nRandom;
		
		Double_t meanErr = sqrt(var/nRandom);
		Double_t varErr = sqrt(var/2/nRandom);

		mean -= parIni[7];
		var -= pow(parIni[3],2);
		gr->SetPoint(iMu, mean, var);
		gr->SetPointError(iMu, meanErr, varErr);
	}

   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kRed);
	TF1* linear = new TF1("linear", "[0]+[1]*x", 0, 5.);   
	TCanvas* canvGR = new TCanvas("canvGR", "var vs mean", 1150, 650);
	gr->Draw("ap");
	gr->Fit(linear);
	linear->Draw("l same");
	
	Double_t gain = linear->GetParameter(1);
	Double_t gainErr = linear->GetParError(1);
	Double_t gainCorr = gain/eqf/enf;
	Double_t gainCorrErr = sqrt(pow(gainErr/eqf/enf,2) +
										 pow(gain/eqf/eqf/enf*eqfErr,2) +
										 pow(gain/eqf/enf/enf*enfErr,2));

	cout << "trueGain: " << parIni[6] << ", gain: " << gain << " +- " << gainErr
		  << ", calculatedGain: " << gainCorr << " +- " << gainCorrErr << endl;

}

//////////////////////////////////////////////////////////////////////////////
Double_t ExpectedSpectrum(Double_t *x, Double_t *par) {
   Int_t kmax = par[9];
   Double_t xx = x[0];
   Double_t alpha = par[0];
   Double_t mu = par[1];
   Double_t lambda = par[2];
   Double_t sigma0 = par[3];
   Double_t sigma1 = par[4];
   Double_t beta = par[5];
   Double_t gain = par[6];
   Double_t ped = par[7];
   Double_t norm = par[8];

	Double_t sumK = 0;
   for(Int_t k=1; k<kmax+1; k++) {
		Double_t sumI = 0;
		Double_t phK = PHk(k,ped,gain);
		Double_t sigmaK = SigmaK(k, sigma0, sigma1);
		for(Int_t i=2; i<k+1; i++) {
			sumI += Binominal(k,i,alpha)*AP(xx,k,i,phK,beta);
		}
		sumK += GP(k,mu,lambda) *
			(Binominal(k,0,alpha)*Gauss(xx,phK,sigmaK) +
			 Binominal(k,1,alpha)*APSingle(xx,k,phK,sigmaK,beta) +
			 sumI);
	}
	return GP(0,mu,lambda)*Gauss(xx,PHk(0,ped,gain),lambda) + sumK;
}
