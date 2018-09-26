//--/----|----/----|----/----|----/----|----/----|----/----|----/----|----/----|
#include "Waveform_LED.h"

Waveform::Waveform() {
}

Waveform::~Waveform(){
   delete[] fTime;
   delete[] fAmp;
   delete fFuncSingle;
}

void Waveform::SetNoise(Double_t noise) {
   TRandom3 rndm(0);
   for (int iPnt = 0; iPnt < fNPnt; iPnt++) {
      fAmp[iPnt] += rndm.Gaus(0, noise);
   }
}

void Waveform::MakeTemplate() {
   fNPnt = fSamplingRate * fIntegRange;
   fAnalysisEnd = fAnalysisStart + fIntegRange;

   fTime = new Double_t[fNPnt];
   fAmp = new Double_t[fNPnt];

   fFuncSingle = new TF1("funcSingle", "[0]/[1]*exp(-x/[1])",
                         fAnalysisStart, fAnalysisEnd);
   fFuncSingle->SetParameters(fTrueGain, fTauDecay);
   fFuncTiming = new TF1("funcTiming", "exp(-x/[0])", 0, 99999);
   fFuncTiming->SetParameter(0, fTauTiming);
   fFuncAPTiming = new TF1("funcAPTiming", "exp(-x/[0])", 0, 99999);
   fFuncAPTiming->SetParameter(0, fTauAP);

   memset(fAmp, 0, sizeof(Double_t)*fNPnt);

   fPntSize = (fAnalysisEnd-fAnalysisStart)/fNPnt;
   for (int iPnt = 0; iPnt < fNPnt; iPnt++) {
      fTime[iPnt] = fAnalysisStart + fPntSize*iPnt;
   }
}


Int_t Waveform::MakeEvent(Int_t nPhoton) {
   TRandom3 rndm(0);
	
   for (size_t iPhoton = 0; iPhoton < nPhoton; iPhoton++) {
      Double_t pulseTime;
      if(fLED) {
         pulseTime = rndm.Uniform(fTauTiming);
      } else {
         pulseTime = fFuncTiming->GetRandom();
      }

      pulseTime += fTriggerTime;

      if(pulseTime < fAnalysisEnd) {
         Hit(pulseTime);
      }
   }

   // dark
   if(fDCR>0) {
      Double_t nExp = fDCR*pow(10,-6) * fIntegRange;
      Int_t nDC = rndm.Poisson(nExp);
      for(Int_t iDC=0; iDC<nDC; iDC++) {
         Double_t pulseTime = rndm.Uniform(fAnalysisStart, fAnalysisEnd);
         Hit(pulseTime);
      }
   }
	
   SetNoise(fNoise);

   return fHit;
}


void Waveform::Hit(Double_t pulseTime) {
   TRandom3 rndm(0);

   fHit++;
   for(Int_t iPnt=0; iPnt<fNPnt; iPnt++) {
		if(fTime[iPnt]<pulseTime) continue;
      fAmp[iPnt] += fFuncSingle->Eval(fTime[iPnt]-pulseTime);
   }

   // prompt crosstalk
   if(fLambda>0) {
      Int_t numCT = rndm.Poisson(fLambda);
      for(Int_t iCT=0; iCT<numCT; iCT++) {
         Hit(pulseTime);
         fNCT++;
      }
   }

   // delayed cross talk                                                      
   if(fAlphaCT>0) {
      if(rndm.Uniform() < fAlphaCT) {
         Double_t dCTTime = fFuncAPTiming->GetRandom();
         Double_t dCTPulseTime = pulseTime + dCTTime;

         if(pulseTime < fAnalysisEnd) {
            for(Int_t iPnt=0; iPnt<fNPnt; iPnt++) {
               if(fTime[iPnt]<dCTPulseTime) continue;
               fAmp[iPnt] += fFuncSingle->Eval(fTime[iPnt]-dCTPulseTime);
            }

            // Crosstalk of delayed crosstalk      
            if(fLambda>0) {
               Int_t numCT = rndm.Poisson(fLambda);
               for(Int_t iCT=0; iCT<numCT; iCT++) {
                  Hit(dCTPulseTime);
               }
            }
         }
      }
   }

   // Afterpulsing
   if(fAlpha>0) {
      if(rndm.Uniform() < fAlpha) {
         Double_t apTime= fFuncAPTiming->GetRandom();
         Double_t apPulseTime = pulseTime + apTime;

         if(pulseTime < fAnalysisEnd) {
            for(Int_t iPnt=0; iPnt<fNPnt; iPnt++) {
					if(fTime[iPnt]<apPulseTime) continue;
               fAmp[iPnt] +=
                  (1-exp(-apTime/fTauDecay)) *
                  fFuncSingle->Eval(fTime[iPnt]-apPulseTime);
            }

            // Cross Talk of After Pulse
            Int_t numCT = rndm.Poisson(fLambda);
            for(Int_t iCT=0; iCT<numCT; iCT++) {
               Hit(pulseTime);
            }
         }
      }
   }
}


void Waveform::Differentiate(Int_t nDiff) {
   if (nDiff>0) {
      Double_t* diffAmp= new Double_t[fNPnt];

      for(int iDiff=0; iDiff<nDiff; iDiff++){
         for (int iPnt = 1; iPnt < fNPnt-1; iPnt++) {
            diffAmp[iPnt] =
                  (fAmp[iPnt+1]-fAmp[iPnt-1]) / fPntSize;
         }
         diffAmp[0+iDiff] = 0;
         diffAmp[fNPnt-1-iDiff] = 0;

         memcpy(fAmp,diffAmp,fNPnt*sizeof(Double_t) );
      }
   }
}


void Waveform::MakeGraph() {
	fGraph = new TGraph;
	for(Int_t iPnt=0; iPnt<fNPnt; iPnt++) {
		fGraph->SetPoint(iPnt, fTime[iPnt], fAmp[iPnt]);
	}
	fGraph->SetMarkerStyle(8);
	fGraph->SetMarkerSize(.5);
	fGraph->Draw("ap");
}

void Waveform::DeleteGraph() {
	delete fGraph;
}

void Waveform::SetSamplingRate(Double_t samplingRate) {
	fSamplingRate = samplingRate;
}

void Waveform::SetIntegRange(Double_t integRange) {
	fIntegRange = integRange;
}

void Waveform::SetTrueGain(Double_t trueGain) {
	fTrueGain = trueGain;
}

void Waveform::SetTauTiming(Double_t tauTiming) {
   fTauTiming = tauTiming;
}

void Waveform::SetNoises(Double_t lambda, Double_t alpha,
                                Double_t alphaCT, Double_t dcr,
                                Double_t noise) {
   fLambda = lambda;
   fAlpha = alpha;
   fAlphaCT = alphaCT;
   fDCR = dcr;
   fNoise = noise;
}

Double_t Waveform::GetCharge() {
   Double_t charge = 0;
   for (int iPnt = 0; iPnt < fNPnt; iPnt++) {
      charge += (fAmp[iPnt])*fPntSize;
   }
   return charge;
}

Double_t Waveform::GetVariance() {
   Double_t variance=0;
   for (int iPnt = 0; iPnt < fNPnt; iPnt++) {
      variance += TMath::Power(fAmp[iPnt],2) * fPntSize;
   }
   return variance;
}

Int_t Waveform::GetNCT() {
   return fNCT;
}
