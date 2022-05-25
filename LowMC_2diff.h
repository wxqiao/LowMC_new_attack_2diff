#ifndef _LOWMC_2DIFF_H
#define _LOWMC_2DIFF_H

#include<vector>
using namespace std;
typedef unsigned long long UINT64;

//s-inverse ddt(2-diff)
const int ddt[64][8] = {
	{0,0,0,0,0,0,0,0},
	{1,5,3,7,0,0,0,0},
	{7,3,6,2,0,0,0,0},
	{2,6,1,5,0,0,0,0},
	{5,7,4,6,0,0,0,0},
	{6,4,3,1,0,0,0,0},
	{3,5,2,4,0,0,0,0},
	{4,2,1,7,0,0,0,0},
	{8,40,24,56,0,0,0,0},
	{9,45,27,63,0,0,0,0},
	{15,11,47,43,30,26,62,58},
	{10,14,46,42,25,29,61,57},
	{13,15,44,46,29,31,60,62},
	{14,12,43,41,28,30,57,59},
	{11,13,42,44,26,28,59,61},
	{12,10,41,47,31,25,58,60},
	{56,24,48,16,0,0,0,0},
	{57,25,61,29,51,19,55,23},
	{63,27,54,18,0,0,0,0},
	{58,30,62,26,49,21,53,17},
	{61,31,60,30,53,23,52,22},
	{62,28,59,25,52,22,49,19},
	{59,29,58,28,50,20,51,21},
	{60,26,57,31,55,17,50,20},
	{16,48,8,40,0,0,0,0},
	{17,49,53,21,11,43,47,15},
	{23,51,55,19,14,42,46,10},
	{18,54,9,45,0,0,0,0},
	{21,55,52,22,13,47,44,14},
	{22,52,51,17,12,46,41,11},
	{19,53,50,20,10,44,43,13},
	{20,50,49,23,15,41,42,12},
	{40,56,32,48,0,0,0,0},
	{41,57,37,53,43,59,39,55},
	{47,59,39,51,46,58,38,50},
	{42,62,38,50,41,61,37,49},
	{45,63,36,54,0,0,0,0},
	{46,60,35,49,44,62,33,51},
	{43,61,34,52,42,60,35,53},
	{44,58,33,55,47,57,34,52},
	{48,32,24,8,0,0,0,0},
	{49,33,29,13,35,51,15,31},
	{55,35,31,11,38,50,14,26},
	{50,38,30,10,33,53,13,25},
	{53,39,28,14,37,55,12,30},
	{54,36,27,9,0,0,0,0},
	{51,37,26,12,34,52,11,29},
	{52,34,25,15,39,49,10,28},
	{24,40,16,32,0,0,0,0},
	{25,41,21,37,19,35,31,47},
	{31,43,23,35,22,34,30,42},
	{26,46,22,34,17,37,29,41},
	{29,47,20,38,21,39,28,46},
	{30,44,19,33,20,38,25,43},
	{27,45,18,36,0,0,0,0},
	{28,42,17,39,23,33,26,44},
	{32,16,8,56,0,0,0,0},
	{33,17,13,61,59,11,23,39},
	{39,19,15,59,62,10,22,34},
	{34,22,14,58,57,13,21,33},
	{37,23,12,62,61,15,20,38},
	{38,20,11,57,60,14,17,35},
	{35,21,10,60,58,12,19,37},
	{36,18,9,63,0,0,0,0}
};

struct matrix {
	int r;
	int c;
	bool ma[256][256];
};

class LowMC {
private:
	int bs;//blocksize
	int ks;//keysize
	int m;//#sbox
	int r;//rounds
	matrix *linear;
	matrix *invLinear;
	matrix constant;
	matrix *keyMa;



public:
	LowMC(int bs, int ks, int m, int r);
	void loadFileFull();
	void encryptFull(bool p[], bool k[], int psize, int ksize, bool c[], int rounds, bool eachRoundOut[][21], bool eachRoundSBoxOut[][21]);
	void matrixMul(matrix& m, bool x[], bool y[]);
	void matrixMul(matrix& m1, matrix& m2);
	int getInactiveNum(bool state1[], bool state2[], int sboxCnt);
	int getActive1Num(bool state[], bool state2[], int sboxCnt);
	int startTestingFullSBoxLayer(bool g_diffr0Sout_1[], bool g_diffr0Sout_2[], bool cdiff1[], bool cdiff2[], bool eachRoundDiff1_0[], bool eachRoundDiff2_0[],
		bool eachRoundDiff1_1[], bool eachRoundDiff2_1[], bool eachRoundDiff1_2[], bool eachRoundDiff2_2[], int t, int j);
	void constructExpressions(bool diff_1[], bool diff_2[], matrix& eq);
	void clearMatrix(matrix &ma);
	void findNext(int boxIndex, bool csdiff1[], bool csdiff2[], bool result1[], bool result2[], matrix& expression, int& total, int& compact, int& totaltime, int z, bool iscorrect_0,
		bool eachRoundDiff1_1[], bool eachRoundDiff2_1[],
		bool eachRoundDiff1_2[], bool eachRoundDiff2_2[]);
	int constructDiffEquations(bool test1[], bool test2[], matrix& expression, bool iscorrect_0,
		bool eachRoundDiff1_1[], bool eachRoundDiff2_1[], bool eachRoundSBoxDiff1_2[], bool eachRoundSBoxDiff2_2[]);
	void gauss(matrix& eqSys);
	void storeSolutions(vector<vector<bool> >& sol, matrix& eqSys, int& solNum);
};

#endif