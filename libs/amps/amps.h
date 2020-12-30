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
	double suub_ttb_(double p1[4], double p2[4], double p3[4], double p4[4]);

	double sgg_ttbg_(double p1[4], double p2[4], double p3[4], double p4[4], double p5[4]);
	double suub_ttbg_(double p1[4], double p2[4], double p3[4], double p4[4], double p5[4]);
	double subu_ttbg_(double p1[4], double p2[4], double p3[4], double p4[4], double p5[4]);
	double sgu_ttbu_(double p1[4], double p2[4], double p3[4], double p4[4], double p5[4]);
	double sgub_ttbub_(double p1[4], double p2[4], double p3[4], double p4[4], double p5[4]);
	double sug_ttbu_(double p1[4], double p2[4], double p3[4], double p4[4], double p5[4]);
	double subg_ttbub_(double p1[4], double p2[4], double p3[4], double p4[4], double p5[4]);
}