#pragma once
struct coupqcd {
	double gg[2];
	double g;
};
struct fermions {
	double fmass[12];
	double fwidth[12];
};
extern struct coupqcd coupqcd_;
extern struct fermions fermions_;
extern "C" {
	double sgg_ttb_(double p1[4], double p2[4], double p3[4], double p4[4]);
}