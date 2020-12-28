extern "C" {
	void bsyppttinit_(double* mT, const int* Nf, const int* SCHEME);

	void bsyggttsq_(double p1[4], double p2[4], double p3[4], double p4[4],
		double* mur2, std::complex<double> amp[3]);
	void bsyqqttsq_(double p1[4], double p2[4], double p3[4], double p4[4],
		double* mur2, std::complex<double> amp[3]);
}