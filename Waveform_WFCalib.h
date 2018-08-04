//--/----|----/----|----/----|----/----|----/----|----/----|----/----|----/----|
// make tf1 object for template waveform
Int_t funcStart = -230;
Int_t funcEnd = -80;
Int_t funcRange = funcEnd - funcStart;

Double_t cfuncSingle(Double_t *x, Double_t *par);
Double_t cfuncTimingPDF(Double_t *x, Double_t *par);

TF1 *funcSingle =
	new TF1("funcSingle", cfuncSingle, funcStart, funcEnd, 3);
// TF1(Char_t name, Double_t func, Double_t min, Double_t max, Int_t nPar);
// single pulse

TF1 *funcTimingPDF = 
	new TF1("funcTimingPDF", cfuncTimingPDF, 0., 500, 1);
// exponential decay (x>=0)

// define cfunctions used in the template
Double_t cfuncSingle(Double_t *x, Double_t *par) {
   //par 0: decay time const
   //par 1: start time of pulse
   //par 2: scale of the single pulse
   //Area normalized to par[2]
   if(x[0] < par[1]) {
		return 0;
	} else {
		return par[2]*TMath::Exp(-(x[0]-par[1])/par[0]) 
			/ (par[0]*(1-TMath::Exp(-4)));
	}
}

Double_t cfuncTimingPDF(Double_t *x, Double_t *par) {
   //par 0: decay time const
   if (x[0] < 0) return 0;
   else return TMath::Exp(-x[0]/par[0]);
}

// define class Waveform
class Waveform{
protected:
   Int_t           fNPnt;
   Double_t        fPntSize;
   Double_t        *fNoiseAmp;
   Double_t        *fSignalAmp;
   Double_t        *fTime;
   Double_t        fAnalysisStart;
   Double_t        fAnalysisEnd;

public:
   Waveform(Int_t nPnt,Double_t analysisStart,Double_t analysisEnd);
   ~Waveform();
   void MakeEvent(Int_t nPhoton);
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
Waveform::Waveform(Int_t nPnt,Double_t analysisStart,Double_t analysisEnd){
   fNPnt   = nPnt;
   fTime      = new Double_t[nPnt];
   fSignalAmp = new Double_t[nPnt];
   fNoiseAmp  = new Double_t[nPnt];
   fAnalysisStart = analysisStart;
   fAnalysisEnd   = analysisEnd;
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

void Waveform::MakeEvent(Int_t nPhoton) {
   for (size_t iPhoton = 0; iPhoton < nPhoton; iPhoton++) {
      Double_t pulseTime = funcTimingPDF->GetRandom();
		Double_t triggerTime = fAnalysisStart + 10.;
		pulseTime += triggerTime + 10.;
      funcSingle->SetParameter(1, pulseTime);

      for (int iPnt = 0; iPnt < fNPnt; iPnt++) {
         fSignalAmp[iPnt] += funcSingle->Eval(fTime[iPnt]);
			// eval is just the calculated value
      }
   }
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

