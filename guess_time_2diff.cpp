#include"LowMC_2diff.h"
#include <iostream>
#include <ctime>

using namespace std;

void testLayerWithFullSBox() {
	LowMC lowmc(21, 21, 7, 4);
	bool pdiff1[21] = { 1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };//choose it
	bool pdiff2[21] = { 1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };//choose it

	bool g1_diffr0Sout_1[21] = { 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	bool g1_diffr0Sout_2[21] = { 0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	bool g2_diffr0Sout_1[21] = { 1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	bool g2_diffr0Sout_2[21] = { 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	bool g3_diffr0Sout_1[21] = { 0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	bool g3_diffr0Sout_2[21] = { 1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	bool g4_diffr0Sout_1[21] = { 0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	bool g4_diffr0Sout_2[21] = { 0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	bool g5_diffr0Sout_1[21] = { 1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	bool g5_diffr0Sout_2[21] = { 0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

	bool g6_diffr0Sout_1[21] = { 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	bool g6_diffr0Sout_2[21] = { 1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

	bool g7_diffr0Sout_1[21] = { 0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	bool g7_diffr0Sout_2[21] = { 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

	bool g8_diffr0Sout_1[21] = { 0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	bool g8_diffr0Sout_2[21] = { 0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	bool p0[21], c0[21], p1[21], c1[21], p2[21], c2[21], k[21];
	bool cdiff1[21], cdiff2[21];
	bool eachRoundOut0[4][21], eachRoundOut1[4][21], eachRoundOut2[4][21], eachRoundDiff1[4][21], eachRoundDiff2[4][21];
	bool eachRoundSBoxOut0[4][21], eachRoundSBoxOut1[4][21], eachRoundSBoxOut2[4][21], eachRoundSBoxDiff1[4][21], eachRoundSBoxDiff2[4][21];
	int psize = 21, ksize = 21;
	int smallerTimes = 0, largerTimes = 0, testTimes = 100;
	for (int test = 0; test < testTimes; test++) {
		cout << endl << "--------Experiments:Times" << dec << test + 1 << "--------"<<endl;

		for (int i = 0; i < psize; i++) {
				p0[i] = rand() % 2;
				k[i] = rand() % 2;
				p1[i] = p0[i] ^ pdiff1[i];
				p2[i] = p0[i] ^ pdiff2[i];
		}

		lowmc.encryptFull(p0, k, psize, ksize, c0, 4, eachRoundOut0, eachRoundSBoxOut0);
		lowmc.encryptFull(p1, k, psize, ksize, c1, 4, eachRoundOut1, eachRoundSBoxOut1);
		lowmc.encryptFull(p2, k, psize, ksize, c2, 4, eachRoundOut2, eachRoundSBoxOut2);
			

		//the internal state diff
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 21; j++) {
				eachRoundDiff1[i][j] = eachRoundOut0[i][j] ^ eachRoundOut1[i][j];
				eachRoundDiff2[i][j] = eachRoundOut0[i][j] ^ eachRoundOut2[i][j];
				eachRoundSBoxDiff1[i][j] = eachRoundSBoxOut0[i][j] ^ eachRoundSBoxOut1[i][j];
				eachRoundSBoxDiff2[i][j] = eachRoundSBoxOut0[i][j] ^ eachRoundSBoxOut2[i][j];
			}
		}


		for (int i = 0; i < psize; i++) {
			cdiff1[i] = c0[i] ^ c1[i];
			cdiff2[i] = c0[i] ^ c2[i];
		}
		
		int t = lowmc.getInactiveNum(eachRoundSBoxDiff1[3],eachRoundSBoxDiff2[3], 7);
		int j = lowmc.getActive1Num(eachRoundSBoxDiff1[3],eachRoundSBoxDiff2[3], 7);



		cout<<endl<<"Guess 1th:"<<endl;
		int f1 = lowmc.startTestingFullSBoxLayer(g1_diffr0Sout_1, g1_diffr0Sout_2, cdiff1, cdiff2,eachRoundDiff1[0], eachRoundDiff2[0],
		 eachRoundDiff1[1], eachRoundDiff2[1],eachRoundSBoxDiff1[2], eachRoundSBoxDiff2[2], t, j);
		if (f1 == 1) {
			smallerTimes++;
		}
		else {
			largerTimes++;
		}


		cout<<endl<<"Guess 2th:"<<endl;
		int f2 = lowmc.startTestingFullSBoxLayer(g2_diffr0Sout_1, g2_diffr0Sout_2, cdiff1, cdiff2,eachRoundDiff1[0], eachRoundDiff2[0],
		 eachRoundDiff1[1], eachRoundDiff2[1],eachRoundSBoxDiff1[2], eachRoundSBoxDiff2[2], t, j);
		if (f2 == 1) {
			smallerTimes++;
		}
		else {
			largerTimes++;
		}


		cout<<endl<<"Guess 3th:"<<endl;
		int f3 = lowmc.startTestingFullSBoxLayer(g3_diffr0Sout_1, g3_diffr0Sout_2, cdiff1, cdiff2,eachRoundDiff1[0], eachRoundDiff2[0],
		 eachRoundDiff1[1], eachRoundDiff2[1],eachRoundSBoxDiff1[2], eachRoundSBoxDiff2[2], t, j);
		if (f3 == 1) {
			smallerTimes++;
		}
		else {
			largerTimes++;
		}


		cout<<endl<<"Guess 4th:"<<endl;
		int f4 = lowmc.startTestingFullSBoxLayer(g4_diffr0Sout_1, g4_diffr0Sout_2, cdiff1, cdiff2,eachRoundDiff1[0], eachRoundDiff2[0],
		 eachRoundDiff1[1], eachRoundDiff2[1],eachRoundSBoxDiff1[2], eachRoundSBoxDiff2[2], t, j);
		if (f4 == 1) {
			smallerTimes++;
		}
		else {
			largerTimes++;
		}


		cout<<endl<<"Guess 5th:"<<endl;
		int f5 = lowmc.startTestingFullSBoxLayer(g5_diffr0Sout_1, g5_diffr0Sout_2, cdiff1, cdiff2,eachRoundDiff1[0], eachRoundDiff2[0],
		 eachRoundDiff1[1], eachRoundDiff2[1],eachRoundSBoxDiff1[2], eachRoundSBoxDiff2[2], t, j);
		if (f5 == 1) {
			smallerTimes++;
		}
		else {
			largerTimes++;
		}


		cout<<endl<<"Guess 6th:"<<endl;
		int f6 = lowmc.startTestingFullSBoxLayer(g6_diffr0Sout_1, g6_diffr0Sout_2, cdiff1, cdiff2,eachRoundDiff1[0], eachRoundDiff2[0],
		 eachRoundDiff1[1], eachRoundDiff2[1],eachRoundSBoxDiff1[2], eachRoundSBoxDiff2[2], t, j);
		if (f6 == 1) {
			smallerTimes++;
		}
		else {
			largerTimes++;
		}


		cout<<endl<<"Guess 7th:"<<endl;
		int f7 = lowmc.startTestingFullSBoxLayer(g7_diffr0Sout_1, g7_diffr0Sout_2, cdiff1, cdiff2,eachRoundDiff1[0], eachRoundDiff2[0],
		 eachRoundDiff1[1], eachRoundDiff2[1],eachRoundSBoxDiff1[2], eachRoundSBoxDiff2[2], t, j);
		if (f7 == 1) {
			smallerTimes++;
		}
		else {
			largerTimes++;
		}


		cout<<endl<<"Guess 8th:"<<endl;
		int f8 = lowmc.startTestingFullSBoxLayer(g8_diffr0Sout_1, g8_diffr0Sout_2, cdiff1, cdiff2,eachRoundDiff1[0], eachRoundDiff2[0],
		 eachRoundDiff1[1], eachRoundDiff2[1],eachRoundSBoxDiff1[2], eachRoundSBoxDiff2[2], t, j);
		if (f8 == 1) {
			smallerTimes++;
		}
		else {
			largerTimes++;
		}
	}
	cout << "testTimes:" << dec << testTimes << endl;
	cout << "Time to #compact diff is smaller than expected:" << dec << smallerTimes << endl;
	cout << "Time to #compact diff is larger than expected:" << dec << largerTimes << endl;
}

int main() {
	testLayerWithFullSBox();
}