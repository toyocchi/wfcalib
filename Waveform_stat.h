//--/----|----/----|----/----|----/----|----/----|----/----|----/----|----/----|
// make tf1 object for template waveform
Int_t funcStart = -230;
Int_t funcEnd = -80;
Int_t funcRange = funcEnd - funcStart;

Double_t cfuncSingle(Double_t *x, Double_t *par);
Double_t cfuncExpPDF(Double_t *x, Double_t *par);

TF1 *funcSingle =
	new TF1("funcSingle", cfuncSingle, funcStart, funcEnd, 3);
// TF1(Char_t name, Double_t func, Double_t min, Double_t max, Int_t nPar);
// single pulse

TF1 *funcExpPDF = 
	new TF1("funcExpPDF", cfuncExpPDF, 0., 500, 1);
// exponential decay (x>=0)

// define cfunctions used in the template
Double_t cfuncSingle(Double_t *x, Double_t *par) {
   //par 0: decay time const
   //par 1: start time of pulse
   //par 2: scale of the single pulse
   //Area normalized to par[2]
   if (x[0] < par[1]) {
		return 0;
	}
   else {
		return par[2]/par[0]*exp(-(x[0]-par[1])/par[0]);
	}
			  
}

Double_t cfuncExpPDF(Double_t *x, Double_t *par) {
   //par 0: decay time const
   if (x[0] < 0) return 0;
   else return TMath::Exp(-x[0]/par[0]);
}

// define class Waveform
class Waveform{
protected:
   Int_t    fNPnt;
	Double_t fPntSize;
   Double_t *fNoiseAmp;
   Double_t *fSignalAmp;
   Double_t *fTime;
   Double_t fAnalysisStart;
   Double_t fAnalysisEnd;
	Double_t fLambda;
	Double_t fAlpha;
	Double_t fAlphaCT;
	Double_t fDCR;
	Int_t fHit;

public:
   Waveform(Int_t nPnt, Double_t analysisStart, Double_t analysisEnd,
				Double_t lambda, Double_t alpha, Double_t alphaCT, Double_t dcr);
   ~Waveform();
   Int_t MakeEvent(Int_t nPhoton);
	void Hit(Double_t pulseTime);
	void Differentiate(Int_t nDiff);
   void SetSignalAmp(Double_t* signalAmp);
	void SetSignalAmpAt(Int_t iPnt,Double_t signalAmp);
	void SetNoiseAmp(Double_t noiseLevel); 
	Double_t GetSignalAmpAt(Int_t iPnt);
	Double_t GetNoiseAmpAt(Int_t iPnt);
   Double_t GetTimeAt(Int_t iPnt);
   Double_t GetCharge();
   Double_t GetTotalVariance();
   Double_t GetSignalVariance();
   Double_t GetNoiseVariance();
   Double_t GetInterference();
};


// implement of class Waveform
Waveform::Waveform(Int_t nPnt, Double_t analysisStart, Double_t analysisEnd,
						 Double_t lambda, Double_t alpha, Double_t alphaCT,
						 Double_t dcr) {
   fNPnt   = nPnt;
   fTime      = new Double_t[nPnt];
   fSignalAmp = new Double_t[nPnt];
   fNoiseAmp  = new Double_t[nPnt];
   fAnalysisStart = analysisStart;
   fAnalysisEnd   = analysisEnd;
	fLambda = lambda;
	fAlpha = alpha;
	fAlphaCT = alphaCT;
	fDCR = dcr;
   memset(fSignalAmp, 0, sizeof(Double_t)*nPnt);
   memset(fNoiseAmp,  0, sizeof(Double_t)*nPnt);
	// initialize fSignalAmp and fNoiseAmp with 0
	// needed memory size is (size of Double_t)*nPnt

	fPntSize = (fAnalysisEnd-fAnalysisStart)/fNPnt;
   for (int iPnt = 0; iPnt < nPnt; iPnt++) {
      fTime[iPnt] = fAnalysisStart + fPntSize*iPnt;
   }
}

Waveform::~Waveform(){
   delete[] fTime;
   delete[] fSignalAmp;
   delete[] fNoiseAmp;
}

void Waveform::Hit(Double_t pulseTime) {
	fHit++;
	funcSingle->SetParameter(1, pulseTime);
	for(Int_t iPnt=0; iPnt<fNPnt; iPnt++) {
		fSignalAmp[iPnt] += funcSingle->Eval(fTime[iPnt]);
	}

	// prompt cross talk
	std::random_device seed;
	std::mt19937 rnd(seed());
	std::poisson_distribution<> poisson(fLambda);
	Int_t numCT = poisson(rnd);
	for(Int_t iCT=0; iCT<numCT; iCT++) {
		Hit(pulseTime);
	}

	// delayed cross talk
	TRandom *gen = new TRandom();
	gen->SetSeed();
	if(gen->Rndm() < fAlphaCT) {
		Double_t tauTiming = funcExpPDF->GetParameter(0);

		Double_t tauAP = 20;
		funcExpPDF->SetParameter(0, tauAP);
		Double_t apTime= funcExpPDF->GetRandom();
		pulseTime += apTime;

		funcExpPDF->SetParameter(0, tauTiming);

		if(pulseTime < fAnalysisEnd) {
			funcSingle->SetParameter(1, pulseTime);
			for(Int_t iPnt=0; iPnt<fNPnt; iPnt++) {
				fSignalAmp[iPnt] += funcSingle->Eval(fTime[iPnt]);
			}
			// Crosstalk of delayed crosstalk
			Int_t numCT = poisson(rnd);
			for(Int_t iCT=0; iCT<numCT; iCT++) {
				Hit(pulseTime);
			}
		}
	}

   // Afterpulsing
	if(gen->Rndm() < fAlpha) {
		Double_t tauTiming = funcExpPDF->GetParameter(0);

		Double_t tauAP = 20;
		funcExpPDF->SetParameter(0, tauAP);
		Double_t apTime= funcExpPDF->GetRandom();
		pulseTime += apTime;

		funcExpPDF->SetParameter(0, tauTiming);

		Double_t tauQuench = funcSingle->GetParameter(0);
		if(pulseTime < fAnalysisEnd) {
			funcSingle->SetParameter(1, pulseTime);
			for(Int_t iPnt=0; iPnt<fNPnt; iPnt++) {
				fSignalAmp[iPnt] 
					+= (1-exp(-apTime/tauQuench)) * funcSingle->Eval(fTime[iPnt]);
			}
			// Cross Talk of After Pulse
			Int_t numCT = poisson(rnd);
			for(Int_t iCT=0; iCT<numCT; iCT++) {
				Hit(pulseTime);
			}
		}
	}
	delete gen;
}

Int_t Waveform::MakeEvent(Int_t nPhoton) {
	fHit = 0;
	Double_t triggerTime = fAnalysisStart + 10.;
   for (size_t iPhoton = 0; iPhoton < nPhoton; iPhoton++) {
      Double_t pulseTime = funcExpPDF->GetRandom();
		pulseTime += triggerTime;
		if(pulseTime < fAnalysisEnd) {
			Hit(pulseTime);
		}
   }

	// dark
	if(fDCR>0) {
		std::random_device seed;
		std::mt19937 rnd(seed());
		Double_t nExp = fDCR*1000/pow(10,9) * (fAnalysisEnd - fAnalysisStart);
		Int_t nDC = gRandom->Poisson(nExp);
		for(Int_t iDC=0; iDC<nDC; iDC++) {
			Double_t pulseTime = gRandom->Uniform(fAnalysisStart, fAnalysisEnd);
			Hit(pulseTime);
		}
	}
	return fHit;
}

void Waveform::Differentiate(Int_t nDiff) {
   if (nDiff>0) {
      Double_t* diffSignalAmp= new Double_t[fNPnt];
      Double_t* diffNoiseAmp= new Double_t[fNPnt];

		for(int iDiff=0; iDiff<nDiff; iDiff++){
			for (int iPnt = 0; iPnt < fNPnt-1; iPnt++) {
				diffSignalAmp[iPnt]=(fSignalAmp[iPnt+1]-fSignalAmp[iPnt])/fPntSize;
				diffNoiseAmp[iPnt]=(fNoiseAmp[iPnt+1]-fNoiseAmp[iPnt])/fPntSize;
			}
			//the last point of diff amp is 0
			diffSignalAmp[fNPnt-1-iDiff] = 0;
			diffNoiseAmp[fNPnt-1-iDiff]  = 0;

			memcpy(fSignalAmp,diffSignalAmp,fNPnt*sizeof(Double_t) );
			memcpy(fNoiseAmp,diffNoiseAmp,fNPnt*sizeof(Double_t) );
		}
   }
}

void Waveform::SetSignalAmp(Double_t* signalAmp) {
	memcpy(fSignalAmp, signalAmp, sizeof(Double_t)*fNPnt);
}

void Waveform::SetSignalAmpAt(Int_t iPnt,Double_t amp) {
	fSignalAmp[iPnt] = amp;
}

void Waveform::SetNoiseAmp(Double_t noiseLevel) {
   for (int iPnt = 0; iPnt < fNPnt; iPnt++) {
      fNoiseAmp[iPnt] = gRandom->Gaus(0, noiseLevel);
   }
}

Double_t Waveform::GetSignalAmpAt(Int_t iPnt) {
   return fSignalAmp[iPnt];
}

Double_t Waveform::GetNoiseAmpAt(Int_t iPnt) {
	return fNoiseAmp[iPnt];
}

Double_t Waveform::GetTimeAt(Int_t iPnt) {
   return fTime[iPnt];
}

Double_t Waveform::GetCharge() {
   Double_t charge = 0;
   for (int iPnt = 0; iPnt < fNPnt; iPnt++) {
      charge += (fSignalAmp[iPnt]+fNoiseAmp[iPnt])*fPntSize;
      // std::cout<<"pnt: "<<fPntSize<<" amp: "<<fSignalAmp[iPnt]<<std::endl;
   }
   return charge;
}

Double_t Waveform::GetTotalVariance() {
   Double_t variance=0;
   for (int iPnt = 0; iPnt < fNPnt; iPnt++) {
      variance += TMath::Power((fSignalAmp[iPnt]+fNoiseAmp[iPnt]),2) * fPntSize;
   }
   return variance;
}

Double_t Waveform::GetSignalVariance(){
   Double_t variance=0;
   for (int iPnt = 0; iPnt < fNPnt; iPnt++) {
      variance += TMath::Power(fSignalAmp[iPnt],2) * fPntSize;
   }
   return variance;
}

Double_t Waveform::GetNoiseVariance(){
   Double_t variance=0;
   for (int iPnt = 0; iPnt < fNPnt; iPnt++) {
      variance += TMath::Power(fNoiseAmp[iPnt],2) * fPntSize;
	}
   return variance;
}

Double_t Waveform::GetInterference(){
   Double_t interference=0;
   for (int iPnt = 0; iPnt < fNPnt; iPnt++) {
      interference += 2 * fSignalAmp[iPnt] * fNoiseAmp[iPnt] * fPntSize;
   }
   return interference;
}

