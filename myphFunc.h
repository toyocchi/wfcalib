Double_t GP(Int_t k, Double_t mu, Double_t lambda) {
	return mu * TMath::Power(mu+k*lambda, k-1) *
		TMath::Exp(-(mu+k*lambda)) / TMath::Factorial(k);
}

Int_t Combination(Int_t n, Int_t r) {
	Int_t num = 1;
	for(Int_t i = 1; i <= r; i++) {
		num = num * (n - i + 1) / i;
		//if *= is used, (Double_t) is necessary
	}
	return num;
}

Double_t Binominal(Int_t k, Int_t i, Double_t alpha) {
	return Combination(k, i)*TMath::Power(alpha,i)*TMath::Power(1-alpha,k-i);
}

Double_t Gauss(Double_t x, Double_t mu, Double_t sigma) {
	return 1./TMath::Sqrt(2*TMath::Pi())/sigma *
			TMath::Exp(-TMath::Power(x-mu,2)/2/sigma/sigma);
}

Double_t PHk(Int_t k, Double_t ped, Double_t gain) {
	return ped + k*gain;
}

// sigma^2+k*sigma_1^2
Double_t SigmaK(Int_t k, Double_t sigma0, Double_t sigma1) {
	return TMath::Sqrt(sigma0*sigma0+k*sigma1*sigma1);
}

Double_t AP(Double_t x, Int_t k, Int_t i, Double_t phK, Double_t beta) {
	if(x>=phK) {
		return TMath::Power(x-phK,i-1)/TMath::Factorial(i-1)/TMath::Power(beta,i) *
			TMath::Exp(-(x-phK)/beta);
	} else {
		return 0;
	}
}

Double_t APSingle(Double_t x, Int_t k,
						Double_t phK, Double_t sigmaK, Double_t beta) {
	Double_t integ = TMath::Sqrt(TMath::Pi()/2.)*sigmaK *
		TMath::Erfc(-(x-phK)/TMath::Sqrt(2)/sigmaK);
	return TMath::Exp(-(x-phK)/beta) /
		(TMath::Sqrt(2*TMath::Pi())*sigmaK*beta) * integ;
}
