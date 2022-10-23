#include"LowMC_2diff.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

using namespace std;

LowMC::LowMC(int BS, int KS, int M, int R) {
	bs = BS;
	ks = KS;
	m = M;
	r = R;

	linear = new matrix[r];
	invLinear = new matrix[r];
	keyMa = new matrix[r + 1];

	loadFileFull();

	srand(time(NULL));
}

void LowMC::loadFileFull() {
	ifstream linearFile(".\\linearFull.txt");
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < bs; j++) {
			for (int k = 0; k < bs; k++) {
					linearFile >> linear[i].ma[j][k];
			}
		}
		linear[i].r = bs;
		linear[i].c = bs;
	}
	linearFile.close();

	//inverse
	ifstream invLinearFile(".\\inverseFull.txt");
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < bs; j++) {
			for (int k = 0; k < bs; k++) {
				invLinearFile >> invLinear[i].ma[j][k];
			}
		}
		invLinear[i].r = bs;
		invLinear[i].c = bs;
	}
	invLinearFile.close();
	//constant
	ifstream constantFile(".\\constantFull.txt");
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < bs; j++) {
			constantFile >> constant.ma[i][j];
		}
	}
	constant.r = r;
	constant.c = bs;
	constantFile.close();
	//key
	ifstream keyFile(".\\keyFull.txt");
	for (int i = 0; i < r + 1; i++) {
		for (int j = 0; j < bs; j++) {
			for (int k = 0; k < bs; k++) {
				keyFile >> keyMa[i].ma[j][k];
			}
		}
		keyMa[i].r = bs;
		keyMa[i].c = ks;
	}
	keyFile.close();
}

void LowMC::encryptFull(bool p0[], bool k[], int psize, int ksize, bool c[], int rounds, bool eachRoundOut[][21], bool eachRoundSBoxOut[][21]) {
	matrix rk;
	rk.r = r + 1;
	rk.c = bs;

	bool* p;
	p = new bool[psize];
	for (int i = 0; i < psize; i++) {
		p[i] = p0[i];
	}

	for (int i = 0; i < r + 1; i++) {
		matrixMul(keyMa[i], k, rk.ma[i]);
	}
	//whitening key
	for (int i = 0; i < psize; i++) {
		p[i] = p[i] ^ rk.ma[0][i];
	}

	int t[3];
	bool* pt = new bool[psize];
	for (int i = 0; i < rounds; i++) {
		//s-box (only works for the first 3m bits)
		for (int j = 0; j < m; j++) {
			t[0] = p[3 * j] ^ (p[1 + 3 * j] & p[2 + 3 * j]);
			t[1] = p[3 * j] ^ p[1 + 3 * j] ^ (p[3 * j] & p[2 + 3 * j]);
			t[2] = p[3 * j] ^ p[1 + 3 * j] ^ p[2 + 3 * j] ^ (p[3 * j] & p[1 + 3 * j]);
			p[3 * j] = t[0];
			p[3 * j + 1] = t[1];
			p[3 * j + 2] = t[2];
		}

		//the value after SBox
		for (int j = 0; j < psize; j++) {
			eachRoundSBoxOut[i][j] = p[j];
		}
		//linear matrix
		matrixMul(linear[i], p, pt);
		
		//constant addition
		for (int j = 0; j < psize; j++) {
			p[j] = pt[j] ^ constant.ma[i][j];
		}

		//key addition
		for (int j = 0; j < psize; j++) {
			p[j] = p[j] ^ rk.ma[i + 1][j];
		}

		for (int j = 0; j < psize; j++) {
			eachRoundOut[i][j] = p[j];
		}
	}
	for (int i = 0; i < psize; i++) {
		c[i] = p[i];
	}
	delete[]pt;
	delete[]p;
}

void LowMC::matrixMul(matrix& m, bool x[], bool y[]) {
	for (int i = 0; i < m.r; i++) {
		y[i] = 0;
		for (int j = 0; j < m.c; j++) {
			y[i] = y[i] ^ (m.ma[i][j] & x[j]);
		}
	}
}

void LowMC::matrixMul(matrix& m1, matrix& m2) {
	matrix m3;
	for (int i = 0; i < m1.r; i++) {
		for (int j = 0; j < m2.c; j++) {
			m3.ma[i][j] = 0;
			for (int k = 0; k < m1.c; k++) {
				m3.ma[i][j] = m3.ma[i][j] ^ (m1.ma[i][k] & m2.ma[k][j]);
			}
		}
	}
	for (int i = 0; i < m1.r; i++) {
		for (int j = 0; j < m2.c; j++) {
			m2.ma[i][j] = m3.ma[i][j];
		}
	}
}

int LowMC::getInactiveNum(bool state1[], bool state2[], int sboxCnt) {
	int cnt = 0;
	for (int i = 0; i < sboxCnt; i++) {
		if (state1[i * 3] == 0
			&& state2[i * 3] == 0
			&& state1[i * 3 + 1] == 0
			&& state2[i * 3 + 1] == 0
			&& state1[i * 3 + 2] == 0
			&& state2[i * 3 + 2] == 0) {
			cnt+=1;
		}

	}
	return cnt;
}

int LowMC::getActive1Num(bool state1[], bool state2[], int sboxCnt) {
	int cnt = 0;
	for (int i = 0; i < sboxCnt; i++) {
		int sum = 32 * state1[3 * i] + 16 * state1[3 * i + 1] + 8 * state1[3 * i + 2] + 4 * state2[3 * i] + 2 * state2[3 * i + 1] + state2[3 * i + 2];
		if (sum == 1 || sum == 2 || sum == 3 || sum == 4 || sum == 5 || sum == 6 || sum == 7 ||
			sum == 8 || sum == 16 || sum == 24 || sum == 32 || sum == 40 || sum == 48 || sum == 56 ||
			sum == 9 || sum == 18 || sum == 36 || sum == 27 || sum == 45 || sum == 54 || sum == 63){
			cnt+=1;
		}
}
	return cnt;
}

int LowMC::startTestingFullSBoxLayer(bool g_diffr0Sout_1[], bool g_diffr0Sout_2[], bool cdiff1[], bool cdiff2[], bool eachRoundDiff1_0[], bool eachRoundDiff2_0[],
	bool eachRoundDiff1_1[], bool eachRoundDiff2_1[],bool eachRoundSBoxDiff1_2[], bool eachRoundSBoxDiff2_2[], int t, int j) {
	matrix expression;
	bool diff_1[21], diff_2[21];
	matrixMul(linear[0], g_diffr0Sout_1, diff_1);
	matrixMul(linear[0], g_diffr0Sout_2, diff_2);

	bool iscorrect_0=true;
	for(int t2=0;t2<bs;t2++){
		if((diff_1[t2]!=eachRoundDiff1_0[t2])||(diff_2[t2]!=eachRoundDiff2_0[t2])){
			iscorrect_0=false;
			break;
		}
	}

	if (iscorrect_0) {
			cout << "The find0:";
			for (int i = 0; i < bs; i++) {
				cout << diff_1[i];
			}
			cout << " ";
			for (int i = 0; i < bs; i++) {
				cout << diff_2[i];
			}
			cout << endl;

			cout << "The corr0:";
		    for (int i = 0; i < bs; i++) {
				cout << eachRoundDiff1_0[i];
			}
			cout << " ";
			for (int i = 0; i < bs; i++) {
				cout << eachRoundDiff2_0[i];
			}
			cout << endl;
	}



	int c = getInactiveNum(diff_1, diff_2, 7);
	int d = getActive1Num(diff_1, diff_2, 7);
	constructExpressions(diff_1, diff_2, expression);

	bool result1[21], result2[21];
	bool csdiff1[21], csdiff2[21];
	matrixMul(invLinear[3], cdiff1, csdiff1);
	matrixMul(invLinear[3], cdiff2, csdiff2);

	int total = 0;//Time to  enumerate diff
	int compact = 0;//the number of iterations to enumerate all  compact diff trails
	int totaltime = 0;//Time to recover key
	int boxIndex = 0;

	int z = (6 * t + j) - (5 * 7 + 3 - 3 * c - d);


	findNext(boxIndex, csdiff1, csdiff2, result1, result2, expression, total, compact, totaltime, z, iscorrect_0,
		eachRoundDiff1_1, eachRoundDiff2_1,
		eachRoundSBoxDiff1_2, eachRoundSBoxDiff2_2);


	int exp0 = pow(2, 3 * m - 3 * t - j);//diff
	int exp1 = 0;
	if (z <= 0) {
		exp1 = pow(2, 2.73 * m - 3 * t - j - 3 * c - d);
	}
	else {
		exp1 = pow(2, 3*t-2.27*m-3);
	}

	
	cout << "Time to enumerate diff backwards:0x" << hex << total << ", expected: 0x" << exp0<<endl;

	cout << "the number of iterations to enumerate all compact diff trails:0x" << hex << compact;
	
	int isSmaller = 1;
	if (compact <= exp0) {
		cout << ", smaller than expected" << endl;
	}
	else {
		isSmaller = 0;
		cout << ", larger than expected" << endl;
	}

	cout << "Time to recover key : 0x" << hex << totaltime << ", expected:0x" << exp1<<endl;
	return isSmaller;
}

void LowMC::constructExpressions(bool diff_1[], bool diff_2[], matrix& eq) {
	int varNum = 0;
	eq.r = 2*bs;
	eq.c = 3 * m;
	bool* con = new bool[2 * bs];
	for (int i = 0; i < 2 * bs; i++) {
		con[i] = 0;
	}
	clearMatrix(eq);

	for (int i = 0; i < m; i++) {
		int a = diff_1[3 * i];
		int b = diff_1[3 * i + 1];
		int c = diff_1[3 * i + 2];
		int d = diff_2[3 * i];
		int e = diff_2[3 * i + 1];
		int f = diff_2[3 * i + 2];
		int s = 32 * a + 16 * b + 8 * c + 4 * d + 2 * e + f;
		if (s == 1) {//(0,0,0,0,0,1)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			con[bs + 3 * i + 2] = 1;
			varNum += 2;
		}
		if (s == 2) {//(0,0,0,0,1,0)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			con[bs + 3 * i + 1] = 1;
			varNum += 2;
		}
		if (s == 3) {//(0,0,0,0,1,1)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum+1] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			con[bs + 3 * i + 2] = 1;
			varNum += 2;
		}
		if (s == 4) {//(0,0,0,1,0,0)
			eq.ma[bs + 3 * i + 1][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			con[bs + 3 * i] = 1;
			varNum += 2;
		}
		if (s == 5) {//(0,0,0,1,0,1)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum] = 1;
			con[bs + 3 * i + 2] = 1;
			varNum += 2;
		}
		if (s == 6) {//(0,0,0,1,1,0)
			eq.ma[bs + 3 * i + 1][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum+1] = 1;
			eq.ma[bs + 3 * i][varNum] = 1;
			con[bs + 3 * i] = 1;
			varNum += 2;
		}
		if (s == 7) {//(0,0,0,1,1,1)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			con[bs + 3 * i + 2] = 1;
			varNum += 2;
		}
		if (s == 8) {//(0,0,1,0,0,0)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			varNum += 2;
		}
		if (s == 9) {//(0,0,1,0,0,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			con[bs + 3 * i + 2] = 1;
			varNum += 2;
		}
		if (s == 10) {//(0,0,1,0,1,0)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i][varNum + 2] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			con[bs + 3 * i + 1] = 1;
			con[bs + 3 * i + 2] = 1;
			varNum += 3;
		}
		if (s == 11) {//(0,0,1,0,1,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i][varNum + 2] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			con[bs + 3 * i + 1] = 1;
			varNum += 3;
		}
		if (s == 12) {//(0,0,1,1,0,0)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 2] = 1;
			eq.ma[bs + 3 * i + 2][varNum] = 1;
			con[3 * i + 2] = 1;
			con[bs + 3 * i] = 1;
			con[bs + 3 * i + 2] = 1;
			varNum += 3;
		}
		if (s == 13) {//(0,0,1,1,0,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 2] = 1;
			con[3 * i + 2] = 1;
			eq.ma[bs + 3 * i][varNum] = 1;
			con[bs + 3 * i] = 1;
			eq.ma[bs + 3 * i + 2][varNum] = 1;
			varNum += 3;
		}
		if (s == 14) {//(0,0,1,1,1,0)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i][varNum + 2] = 1;
			con[3 * i + 2] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 2] = 1;
			con[bs + 3 * i + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			con[bs + 3 * i + 2] = 1;
			varNum += 3;
		}
		if (s == 15) {//(0,0,1,1,1,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i][varNum + 2] = 1;
			con[3 * i + 2] = 1;
			eq.ma[bs + 3 * i + 2][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 1][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 2] = 1;
			con[bs + 3 * i + 1] = 1;
			varNum += 3;
		}
		if (s == 16) {//(0,1,0,0,0,0)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[3 * i + 1] = 1;
			varNum += 2;
		}
		if (s == 17) {//(0,1,0,0,0,1)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i][varNum + 2] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[bs + 3 * i + 2] = 1;
			con[3 * i + 1] = 1;
			con[3 * i + 2] = 1;
			varNum += 3;
		}
		if (s == 18) {//(0,1,0,0,1,0)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[3 * i + 1] = 1;
			con[bs + 3 * i + 1] = 1;
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			varNum += 2;
		}
		if (s == 19) {//(0,1,0,0,1,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			eq.ma[bs + 3 * i][varNum + 2] = 1;
			con[3 * i + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			con[bs + 3 * i + 2] = 1;
			eq.ma[bs + 3 * i + 1][varNum+1] = 1;
			varNum += 3;
		}
		if (s == 20) {//(0,1,0,1,0,0)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 2] = 1;
			con[3 * i + 1] = 1;
			con[bs + 3 * i] = 1;
			eq.ma[bs + 3 * i + 1][varNum] = 1;
			con[bs + 3 * i + 1] = 1;
			varNum += 3;
		}
		if (s == 21) {//(0,1,0,1,0,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			eq.ma[bs + 3 * i][varNum + 2] = 1;
			con[3 * i + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 2] = 1;
			con[bs + 3 * i + 2] = 1;
			eq.ma[bs + 3 * i + 1][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			con[bs + 3 * i + 1] = 1;
			varNum += 3;
		}
		if (s == 22) {//(0,1,0,1,1,0)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 2] = 1;
			con[3 * i + 1] = 1;
			eq.ma[bs + 3 * i][varNum] = 1;
			con[bs + 3 * i] = 1;
			eq.ma[bs + 3 * i + 1][varNum] = 1;
			varNum += 3;
		}
		if (s == 23) {//(0,1,0,1,1,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			eq.ma[bs + 3 * i][varNum + 2] = 1;
			con[3 * i + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 2] = 1;
			con[bs + 3 * i + 2] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 1][varNum] = 1;
			varNum += 3;
		}
		if (s == 24) {//(0,1,1,0,0,0)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			varNum += 2;
		}
		if (s == 25) {//(0,1,1,0,0,1)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i][varNum + 2] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[bs + 3 * i + 2] = 1;
			con[3 * i + 1] = 1;
			varNum += 3;
		}
		if (s == 26) {//(0,1,1,0,1,0)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			eq.ma[3 * i][varNum + 2] = 1;
			con[bs + 3 * i + 1] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			varNum += 3;
		}
		if (s == 27) {//(0,1,1,0,1,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			con[bs + 3 * i + 2] = 1;
			varNum += 2;
		}
		if (s == 28) {//(0,1,1,1,0,0)
			eq.ma[bs + 3 * i + 1][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			eq.ma[3 * i + 1][varNum + 2] = 1;
			con[bs + 3 * i] = 1;
			eq.ma[3 * i + 2][varNum + 2] = 1;
			con[3 * i + 2] = 1;
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i][varNum + 1] = 1;
			con[3 * i] = 1;
			varNum += 3;
		}
		if (s == 29) {//(0,1,1,1,0,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i][varNum + 2] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 2] = 1;
			con[bs + 3 * i + 2] = 1;
			eq.ma[bs + 3 * i + 1][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 2] = 1;
			con[bs + 3 * i + 1] = 1;
			varNum += 3;
		}
		if (s == 30) {//(0,1,1,1,1,0)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i][varNum + 2] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 2] = 1;
			con[bs + 3 * i + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 2] = 1;
			varNum += 3;
		}
		if (s == 31) {//(0,1,1,1,1,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 2] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 2] = 1;
			eq.ma[bs + 3 * i + 2][varNum] = 1;
			eq.ma[bs + 3 * i][varNum] = 1;
			con[bs + 3 * i] = 1;
			varNum += 3;
		}
		if (s == 32) {//(1,0,0,0,0,0)
			eq.ma[3 * i + 1][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[3 * i] = 1;
			varNum += 2;
		}
		if (s == 33) {//(1,0,0,0,0,1)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 1][varNum + 2] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			con[bs + 3 * i + 2] = 1;
			con[3 * i] = 1;
			con[3 * i + 2] = 1;
			varNum += 3;
		}
		if (s == 34) {//(1,0,0,0,1,0)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum + 2] = 1;
			con[bs + 3 * i + 1] = 1;
			con[3 * i] = 1;
			eq.ma[3 * i + 1][varNum] = 1;
			con[3 * i + 1] = 1;
			varNum += 3;
		}
		if (s == 35) {//(1,0,0,0,1,1)
			eq.ma[3 * i + 1][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 2] = 1;
			con[3 * i] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 2] = 1;
			con[bs + 3 * i + 2] = 1;
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i][varNum + 1] = 1;
			con[bs + 3 * i] = 1;
			varNum += 3;
		}
		if (s == 36) {//(1,0,0,1,0,0)
			eq.ma[3 * i + 1][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[3 * i] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 1][varNum] = 1;
			con[bs + 3 * i] = 1;
			varNum += 2;
		}
		if (s == 37) {//(1,0,0,1,0,1)
			eq.ma[3 * i + 1][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 2] = 1;
			con[3 * i] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			con[bs + 3 * i + 2] = 1;
			eq.ma[bs + 3 * i][varNum + 1] = 1;
			varNum += 3;
		}
		if (s == 38) {//(1,0,0,1,1,0)
			eq.ma[3 * i + 1][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 2] = 1;
			con[3 * i] = 1;
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum] = 1;
			con[bs + 3 * i + 1] = 1;
			varNum += 3;
		}
		if (s == 39) {//(1,0,0,1,1,1)
			eq.ma[3 * i + 1][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 2] = 1;
			con[3 * i] = 1;
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 2] = 1;
			con[bs + 3 * i + 2] = 1;
			varNum += 3;
		}
		if (s == 40) {//(1,0,1,0,0,0)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			con[3 * i + 2] = 1;
			varNum += 2;
		}
		if (s == 41) {//(1,0,1,0,0,1)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 1][varNum + 2] = 1;
			con[bs + 3 * i + 2] = 1;
			eq.ma[3 * i][varNum] = 1;
			con[3 * i] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			varNum += 3;
		}
		if (s == 42) {//(1,0,1,0,1,0)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			eq.ma[3 * i][varNum + 2] = 1;
			con[bs + 3 * i + 1] = 1;
			eq.ma[3 * i + 2][varNum + 2] = 1;
			con[3 * i + 2] = 1;
			eq.ma[3 * i + 1][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			con[3 * i + 1] = 1;
			varNum += 3;
		}
		if (s == 43) {//(1,0,1,0,1,1)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i][varNum + 2] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			con[bs + 3 * i + 2] = 1;
			eq.ma[3 * i + 2][varNum + 2] = 1;
			con[3 * i + 2] = 1;
			eq.ma[3 * i + 1][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 1][varNum + 2] = 1;
			con[3 * i + 1] = 1;
			varNum += 3;
		}
		if (s == 44) {//(1,0,1,1,0,0)
			eq.ma[bs + 3 * i + 1][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			eq.ma[3 * i + 1][varNum + 2] = 1;
			con[bs + 3 * i] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			eq.ma[3 * i][varNum + 1] = 1;
			varNum += 3;
		}
		if (s == 45) {//(1,0,1,1,0,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			con[3 * i + 2] = 1;
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum] = 1;
			con[bs + 3 * i + 2] = 1;
			varNum += 2;
		}
		if (s == 46) {//(1,0,1,1,1,0)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i][varNum + 2] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			con[3 * i + 2] = 1;
			eq.ma[bs + 3 * i + 2][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 2] = 1;
			con[bs + 3 * i + 2] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 2] = 1;
			con[bs + 3 * i + 1] = 1;
			varNum += 3;
		}
		if (s == 47) {//(1,0,1,1,1,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i][varNum + 2] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			con[3 * i + 2] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			con[bs + 3 * i + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 2] = 1;
			eq.ma[bs + 3 * i + 2][varNum+1] = 1;
			varNum += 3;
		}
		if (s == 48) {//(1,1,0,0,0,0)
			eq.ma[3 * i + 1][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			eq.ma[3 * i][varNum] = 1;
			con[3 * i] = 1;
			varNum += 2;
		}
		if (s == 49) {//(1,1,0,0,0,1)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i][varNum + 2] = 1;
			con[bs + 3 * i + 2] = 1;
			eq.ma[3 * i + 1][varNum + 2] = 1;
			con[3 * i + 1] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			varNum += 3;
		}
		if (s == 50) {//(1,1,0,0,1,0)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum + 2] = 1;
			con[bs + 3 * i + 1] = 1;
			eq.ma[3 * i][varNum] = 1;
			con[3 * i] = 1;
			eq.ma[3 * i + 1][varNum] = 1;
			varNum += 3;
		}
		if (s == 51) {//(1,1,0,0,1,1)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i][varNum + 2] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			con[bs + 3 * i + 2] = 1;
			eq.ma[3 * i + 1][varNum + 2] = 1;
			con[3 * i + 1] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum + 2] = 1;
			varNum += 3;
		}
		if (s == 52) {//(1,1,0,1,0,0)
			eq.ma[bs + 3 * i + 1][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum + 2] = 1;
			con[bs + 3 * i] = 1;
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum] = 1;
			con[3 * i + 1] = 1;
			varNum += 3;
		}
		if (s == 53) {//(1,1,0,1,0,1)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i][varNum + 2] = 1;
			eq.ma[bs + 3 * i + 2][varNum] = 1;
			con[bs + 3 * i + 2] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum + 2] = 1;
			con[3 * i + 2] = 1;
			eq.ma[3 * i + 1][varNum + 2] = 1;
			con[3 * i + 1] = 1;
			varNum += 3;
		}
		if (s == 54) {//(1,1,0,1,1,0)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			eq.ma[3 * i + 1][varNum] = 1;
			con[3 * i + 1] = 1;
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum] = 1;
			con[bs + 3 * i + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			varNum += 2;
		}
		if (s == 55) {//(1,1,0,1,1,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			eq.ma[bs + 3 * i][varNum + 2] = 1;
			eq.ma[3 * i + 1][varNum] = 1;
			con[3 * i + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			con[bs + 3 * i + 2] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 2] = 1;
			varNum += 3;
		}
		if (s == 56) {//(1,1,1,0,0,0)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			varNum += 2;
		}
		if (s == 57) {//(1,1,1,0,0,1)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i][varNum + 2] = 1;
			con[bs + 3 * i + 2] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			eq.ma[3 * i + 1][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 1][varNum + 2] = 1;
			con[3 * i + 1] = 1;
			varNum += 3;
		}
		if (s == 58) {//(1,1,1,0,1,0)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			eq.ma[3 * i][varNum + 2] = 1;
			con[bs + 3 * i + 1] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum + 2] = 1;
			con[3 * i + 2] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 1][varNum] = 1;
			varNum += 3;
		}
		if (s == 59) {//(1,1,1,0,1,1)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 1][varNum + 2] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			con[bs + 3 * i + 2] = 1;
			eq.ma[3 * i + 2][varNum + 2] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			eq.ma[3 * i][varNum] = 1;
			con[3 * i] = 1;
			varNum += 3;
		}
		if (s == 60) {//(1,1,1,1,0,0)
			eq.ma[bs + 3 * i + 1][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			eq.ma[3 * i + 1][varNum + 2] = 1;
			con[bs + 3 * i] = 1;
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum + 2] = 1;
			con[3 * i + 2] = 1;
			varNum += 3;
		}
		if (s == 61) {//(1,1,1,1,0,1)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i][varNum + 2] = 1;
			eq.ma[bs + 3 * i + 2][varNum] = 1;
			con[bs + 3 * i + 2] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			con[3 * i + 1] = 1;
			eq.ma[3 * i + 2][varNum + 2] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			varNum += 3;
		}
		if (s == 62) {//(1,1,1,1,1,0)
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			eq.ma[3 * i][varNum + 2] = 1;
			eq.ma[bs + 3 * i + 1][varNum] = 1;
			con[bs + 3 * i + 1] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 1][varNum + 2] = 1;
			varNum += 3;
		}
		if (s == 63) {//(1,1,1,1,1,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i][varNum] = 1;
			eq.ma[bs + 3 * i + 1][varNum + 1] = 1;
			eq.ma[bs + 3 * i + 2][varNum] = 1;
			eq.ma[bs + 3 * i + 2][varNum + 1] = 1;
			con[bs + 3 * i + 2] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			varNum += 2;
		}
	}
	
	eq.c = varNum + 1;
	for (int i = 0; i < 2 * bs; i++) {
		eq.ma[i][eq.c - 1] = con[i];
	}
	matrix eq1, eq2;
	eq1.r = bs;
	eq1.c = eq.c;
	eq2.r = bs;
	eq2.c = eq.c;
	for (int i = 0; i < bs; i++) {
		for (int j = 0; j < eq1.c; j++) {
			eq1.ma[i][j] = eq.ma[i][j];
		}
	}
	for (int i = 0; i < bs; i++) {
		for (int j = 0; j < eq2.c; j++) {
			eq2.ma[i][j] = eq.ma[bs + i][j];
		}
	}

	matrixMul(linear[1], eq1);

	matrixMul(linear[1], eq2);

	for (int i = 0; i < bs; i++) {
		for (int j = 0; j < eq.c; j++) {
			eq.ma[i][j] = eq1.ma[i][j];
		}
	}
	for (int i = 0; i < bs; i++) {
		for (int j = 0; j < eq.c; j++) {
			eq.ma[bs+i][j] = eq2.ma[i][j];
		}
	}

	delete[]con;
}

void LowMC::clearMatrix(matrix& ma) {
	for (int i = 0; i < ma.r; i++) {
		for (int j = 0; j < ma.c; j++) {
			ma.ma[i][j] = 0;
		}
	}
}

void LowMC::findNext(int boxIndex, bool csdiff1[], bool csdiff2[],bool result1[],bool result2[], matrix& expression, int &total, int &compact, int &totaltime, int z,bool iscorrect_0,
	bool eachRoundDiff1_1[], bool eachRoundDiff2_1[],
	bool eachRoundSBoxDiff1_2[], bool eachRoundSBoxDiff2_2[]) {
	if (boxIndex == m) {
		bool test1[21], test2[21];
		matrixMul(invLinear[2], result1, test1);
		matrixMul(invLinear[2], result2, test2);
		int a = getInactiveNum(test1, test2, 7);
		int b = getActive1Num(test1, test2, 7);
		int cnt = constructDiffEquations(test1, test2, expression,iscorrect_0,
			eachRoundDiff1_1, eachRoundDiff2_1,eachRoundSBoxDiff1_2, eachRoundSBoxDiff2_2);
		int time = 0;
		time = pow(2, 2 * a);
		if (z > 0) {
			int v = pow(2,z);
			time = time * v;
		}
		time = time * cnt;

		total += 1;

		compact += cnt;
		
		totaltime += time;

		return;
	}
	else {
		int t = boxIndex * 3;
		int x0 = csdiff1[t];
		int x1 = csdiff1[t + 1];
		int x2 = csdiff1[t + 2];
		int x3 = csdiff2[t];
		int x4 = csdiff2[t + 1];
		int x5 = csdiff2[t + 2];
		int sum = 32 * x0 + 16 * x1 + 8 * x2 + 4 * x3 + 2 * x4 + x5;
		if (sum == 0) {
			result1[t] = 0;
			result1[t + 1] = 0;
			result1[t + 2] = 0;
			result2[t] = 0;
			result2[t + 1] = 0;
			result2[t + 2] = 0;
			findNext(boxIndex + 1, csdiff1, csdiff2, result1, result2, expression, total, compact, totaltime, z,iscorrect_0,
				eachRoundDiff1_1,eachRoundDiff2_1,
				eachRoundSBoxDiff1_2,eachRoundSBoxDiff2_2);
		}
		else {
			if (sum == 1 || sum == 2 || sum == 3 || sum == 4 || sum == 5 || sum == 6 || sum == 7 ||
				sum == 8 || sum == 16 || sum == 24 || sum == 32 || sum == 40 || sum == 48 || sum == 56 ||
				sum == 9 || sum == 18 || sum == 36 || sum == 27 || sum == 45 || sum == 54 || sum == 63) {
				for (int i = 0; i < 4; i++) {
					result1[t] = (ddt[sum][i] >> 5) & 0x1;
					result1[t + 1] = (ddt[sum][i] >> 4) & 0x1;
					result1[t + 2] = (ddt[sum][i] >> 3) & 0x1;
					result2[t] = (ddt[sum][i] >> 2) & 0x1;
					result2[t + 1] = (ddt[sum][i] >> 1) & 0x1;
					result2[t + 2] = (ddt[sum][i]) & 0x1;
					findNext(boxIndex + 1, csdiff1, csdiff2, result1, result2, expression, total, compact, totaltime, z,iscorrect_0,
						eachRoundDiff1_1, eachRoundDiff2_1,
						eachRoundSBoxDiff1_2, eachRoundSBoxDiff2_2);
				}
			}
			else {
				for (int i = 0; i < 8; i++) {
					result1[t] = (ddt[sum][i] >> 5) & 0x1;
					result1[t + 1] = (ddt[sum][i] >> 4) & 0x1;
					result1[t + 2] = (ddt[sum][i] >> 3) & 0x1;
					result2[t] = (ddt[sum][i] >> 2) & 0x1;
					result2[t + 1] = (ddt[sum][i] >> 1) & 0x1;
					result2[t + 2] = (ddt[sum][i]) & 0x1;
					findNext(boxIndex + 1, csdiff1, csdiff2, result1, result2, expression, total, compact, totaltime, z,iscorrect_0,
						eachRoundDiff1_1, eachRoundDiff2_1,
						eachRoundSBoxDiff1_2, eachRoundSBoxDiff2_2);
				}

			}
		}
	}
}

int LowMC::constructDiffEquations(bool test1[], bool test2[], matrix& expression,bool iscorrect_0,
	bool eachRoundDiff1_1[], bool eachRoundDiff2_1[], bool eachRoundSBoxDiff1_2[], bool eachRoundSBoxDiff2_2[]) {
	matrix eq;
	int currentRow = 0;
	eq.c = expression.c;
	eq.r = 6 * bs;
	clearMatrix(eq);
	for (int i = 0; i < m; i++) {
		int z0 = test1[3 * i];
		int z1 = test1[3 * i + 1];
		int z2 = test1[3 * i + 2];
		int z3 = test2[3 * i];
		int z4 = test2[3 * i + 1];
		int z5 = test2[3 * i + 2];
		int w = 32 * z0 + 16 * z1 + 8 * z2 + 4 * z3 + 2 * z4 + z5;
		if (w == 0) {//(0,0,0,0,0,0)
			for (int j = 3 * i; j < 3 * i + 3; j++) {
				for (int t = 0; t < expression.c; t++) {
					eq.ma[currentRow][t] ^= expression.ma[j][t];
					eq.ma[currentRow + 1][t] ^= expression.ma[bs + j][t];
				}
				currentRow += 2;
			}
		}
		if (w == 1) {//(0,0,0,0,0,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] ^= expression.ma[3 * i][t];
				eq.ma[currentRow + 1][t] ^= expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 2][t] ^= expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[bs + 3 * i + 2][t];	
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 2) {//(0,0,0,0,1,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] ^= expression.ma[3 * i][t];
				eq.ma[currentRow + 1][t] ^= expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 2][t] ^= expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[bs + 3 * i + 1][t];		
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 3) {//(0,0,0,0,1,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] ^= expression.ma[3 * i][t];
				eq.ma[currentRow + 1][t] ^= expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 2][t] ^= expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[bs + 3 * i + 2][t];	
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 4) {//(0,0,0,1,0,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] ^= expression.ma[3 * i][t];
				eq.ma[currentRow + 1][t] ^= expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 2][t] ^= expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[bs + 3 * i][t];		
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 5) {//(0,0,0,1,0,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] ^= expression.ma[3 * i][t];
				eq.ma[currentRow + 1][t] ^= expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 2][t] ^= expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[bs + 3 * i + 2][t];		
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 6) {//(0,0,0,1,1,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] ^= expression.ma[3 * i][t];
				eq.ma[currentRow + 1][t] ^= expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 2][t] ^= expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[bs + 3 * i + 1][t];	
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 7) {//(0,0,0,1,1,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] ^= expression.ma[3 * i][t];
				eq.ma[currentRow + 1][t] ^= expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 2][t] ^= expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[bs + 3 * i + 2][t];
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 8) {//(0,0,1,0,0,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] ^= expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] ^= expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 2][t] ^= expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[3 * i + 2][t];	
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 9) {//(0,0,1,0,0,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i][t]
					^ expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 2][t]
					^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 3][t] = eq.ma[currentRow + 3][t] ^ expression.ma[3 * i + 2][t];
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 10) {//(0,0,1,0,1,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i + 2][t]
					^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 2][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 11) {//(0,0,1,0,1,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 2][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 12) {//(0,0,1,1,0,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i + 2][t]
					^ expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i][t]
					^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 2][t];	
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 13) {//(0,0,1,1,0,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i][t]
					^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i][t]
					^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 2][t];	
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 14) {//(0,0,1,1,1,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i + 2][t]
					^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i][t]
					^ expression.ma[3 * i + 1][t] ^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 2][t];	
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 15) {//(0,0,1,1,1,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 2][t];	
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 16) {//(0,1,0,0,0,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] ^= expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] ^= expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 2][t] ^= expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[3 * i + 1][t];	
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 17) {//(0,1,0,0,0,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i + 2][t]
					^ expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i + 1][t]
					^ expression.ma[bs + 3 * i + 2][t] ^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i + 2][t];	
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 18) {//(0,1,0,0,1,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i][t]
					^ expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 2][t]
					^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 3][t] = eq.ma[currentRow + 3][t] ^ expression.ma[3 * i + 1][t];
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 19) {//(0,1,0,0,1,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i + 2][t]
					^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 1][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 20) {//(0,1,0,1,0,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i][t]
					^ expression.ma[3 * i + 1][t] ^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 1][t];	
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 21) {//(0,1,0,1,0,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i][t]
					^ expression.ma[3 * i + 1][t] ^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 1][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 22) {//(0,1,0,1,1,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i][t]
					^ expression.ma[3 * i + 1][t] ^ expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i][t]
					^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 1][t];	
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 23) {//(0,1,0,1,1,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i][t]
					^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i][t]
					^ expression.ma[3 * i + 1][t] ^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 1][t];	
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 24) {//(0,1,1,0,0,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] ^= expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] ^= expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 2][t] ^= expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[3 * i + 2][t];	
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 25) {//(0,1,1,0,0,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i + 1][t]
					^ expression.ma[bs + 3 * i + 2][t] ^ expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i + 1][t]
					^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i + 2][t];	
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 26) {//(0,1,1,0,1,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i + 2][t]
					^ expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i + 1][t]
					^ expression.ma[bs + 3 * i + 2][t] ^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i + 1][t];	
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 27) {//(0,1,1,0,1,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i][t]
					^ expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 2][t]
					^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 3][t] = eq.ma[currentRow + 3][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[3 * i + 2][t];	
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 28) {//(0,1,1,1,0,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 1][t]
					^ expression.ma[bs + 3 * i + 2][t] ^ expression.ma[3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[3 * i + 1][t] ^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i][t];	
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 29) {//(0,1,1,1,0,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 2][t]
					^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i + 1][t] ^ expression.ma[3 * i + 2][t]
					^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 1][t] ^ expression.ma[3 * i + 2][t];	
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 30) {//(0,1,1,1,1,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i + 1][t] ^ expression.ma[3 * i + 2][t]
					^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 1][t] ^ expression.ma[3 * i + 2][t];	
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 31) {//(0,1,1,1,1,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i][t]
					^ expression.ma[bs + 3 * i + 1][t] ^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[3 * i + 2][t];	
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 32) {//(1,0,0,0,0,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] ^= expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] ^= expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 2][t] ^= expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[3 * i][t];
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 33) {//(1,0,0,0,0,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i + 2][t]
					^ expression.ma[3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[bs + 3 * i + 2][t] ^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i + 2][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 34) {//(1,0,0,0,1,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i + 1][t]
					^ expression.ma[3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[bs + 3 * i + 1][t] ^ expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i + 1][t];	
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;

		}
		if (w == 35) {//(1,0,0,0,1,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i][t]
					^ expression.ma[bs + 3 * i + 1][t] ^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i][t];	
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 36) {//(1,0,0,1,0,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i][t] ^ expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i + 1][t] ^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 3][t] = eq.ma[currentRow + 3][t] ^ expression.ma[3 * i][t];
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 37) {//(1,0,0,1,0,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 2][t]
					^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 38) {//(1,0,0,1,1,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i][t]
					^ expression.ma[3 * i + 1][t] ^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 39) {//(1,0,0,1,1,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i][t]
					^ expression.ma[3 * i + 1][t] ^ expression.ma[3 * i + 2][t]
					^ expression.ma[bs + 3 * i + 1][t] ^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 40) {//(1,0,1,0,0,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] ^= expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] ^= expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 2][t] ^= expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[3 * i][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[3 * i + 2][t];
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 41) {//(1,0,1,0,0,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[bs + 3 * i + 2][t] ^ expression.ma[3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i + 2][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 42) {//(1,0,1,0,1,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[bs + 3 * i + 1][t] ^ expression.ma[bs + 3 * i + 2][t] ^ expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i + 1][t]
					^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i + 1][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 43) {//(1,0,1,0,1,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 2][t]
					^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i + 1][t] ^ expression.ma[bs + 3 * i + 2][t]
					^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i + 1][t] ^ expression.ma[bs + 3 * i + 2][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 44) {//(1,0,1,1,0,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i + 2][t] ^ expression.ma[3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 2][t]
					^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 45) {//(1,0,1,1,0,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i][t]
					^ expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 2][t]
					^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 3][t] = eq.ma[currentRow + 3][t] ^ expression.ma[3 * i][t]
					^ expression.ma[3 * i + 2][t];
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 46) {//(1,0,1,1,1,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i][t]
					^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i + 1][t] ^ expression.ma[3 * i + 2][t]
					^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 2][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 47) {//(1,0,1,1,1,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i][t]
					^ expression.ma[3 * i + 1][t] ^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i + 1][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 2][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 48) {//(1,1,0,0,0,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] ^= expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] ^= expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 2][t] ^= expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[3 * i][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[3 * i + 1][t];
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 49) {//(1,1,0,0,0,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i + 2][t]
					^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[bs + 3 * i + 1][t] ^ expression.ma[bs + 3 * i + 2][t] ^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i + 2][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 50) {//(1,1,0,0,1,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[bs + 3 * i + 1][t] ^ expression.ma[3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i + 1][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 51) {//(1,1,0,0,1,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i + 1][t] ^ expression.ma[bs + 3 * i + 2][t]
					^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 1][t]
					^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i + 1][t] ^ expression.ma[bs + 3 * i + 2][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 52) {//(1,1,0,1,0,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i + 1][t]
					^ expression.ma[3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[bs + 3 * i + 1][t] ^ expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 53) {//(1,1,0,1,0,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[bs + 3 * i + 2][t] ^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i + 1][t] ^ expression.ma[bs + 3 * i + 2][t]
					^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 2][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 54) {//(1,1,0,1,1,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i][t] ^ expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 2][t]
					^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 3][t] = eq.ma[currentRow + 3][t] ^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 1][t];
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 55) {//(1,1,0,1,1,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i + 2][t]
					^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 1][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 56) {//(1,1,1,0,0,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] ^= expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] ^= expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 2][t] ^= expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[3 * i][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 3][t] ^= expression.ma[3 * i + 2][t];
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
		if (w == 57) {//(1,1,1,0,0,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 1][t]
					^ expression.ma[bs + 3 * i + 2][t] ^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 1][t]
					^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i + 2][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 58) {//(1,1,1,0,1,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[bs + 3 * i + 2][t] ^ expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[bs + 3 * i + 1][t] ^ expression.ma[bs + 3 * i + 2][t] ^ expression.ma[3 * i][t]
					^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i + 1][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 59) {//(1,1,1,0,1,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 1][t]
					^ expression.ma[bs + 3 * i + 2][t] ^ expression.ma[3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[3 * i + 1][t] ^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i + 1][t]
					^ expression.ma[bs + 3 * i + 2][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 60) {//(1,1,1,1,0,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i + 1][t]
					^ expression.ma[bs + 3 * i + 2][t] ^ expression.ma[3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[bs + 3 * i + 1][t] ^ expression.ma[bs + 3 * i + 2][t]
					^ expression.ma[3 * i + 1][t] ^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 61) {//(1,1,1,1,0,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i][t]
					^ expression.ma[bs + 3 * i + 1][t] ^ expression.ma[bs + 3 * i + 2][t] ^ expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i + 1][t] ^ expression.ma[3 * i][t]
					^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 2][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 62) {//(1,1,1,1,1,0)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[bs + 3 * i + 2][t]
					^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 1][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 1][t]
					^ expression.ma[bs + 3 * i + 2][t] ^ expression.ma[3 * i + 2][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[bs + 3 * i][t] ^ expression.ma[bs + 3 * i + 1][t];
			}
			eq.ma[currentRow + 2][eq.c - 1] ^= 1;
			currentRow += 3;
		}
		if (w == 63) {//(1,1,1,1,1,1)
			for (int t = 0; t < expression.c; t++) {
				eq.ma[currentRow][t] = eq.ma[currentRow][t] ^ expression.ma[3 * i][t] ^ expression.ma[bs + 3 * i][t];
				eq.ma[currentRow + 1][t] = eq.ma[currentRow + 1][t] ^ expression.ma[3 * i + 1][t] ^ expression.ma[bs + 3 * i + 1][t];
				eq.ma[currentRow + 2][t] = eq.ma[currentRow + 2][t] ^ expression.ma[3 * i + 2][t] ^ expression.ma[bs + 3 * i + 2][t];
				eq.ma[currentRow + 3][t] = eq.ma[currentRow + 3][t] ^ expression.ma[3 * i][t] ^ expression.ma[3 * i + 1][t]
					^ expression.ma[3 * i + 2][t];
			}
			eq.ma[currentRow + 3][eq.c - 1] ^= 1;
			currentRow += 4;
		}
	}
	eq.r = currentRow;

	gauss(eq);

	vector<vector<bool> >sol;
	sol.clear();
	int solNum = 0;
	storeSolutions(sol, eq, solNum);

	if (iscorrect_0) {
		bool iscorrect_2 = true;
		for (int t2 = 0; t2 < bs; t2++) {
			if (test1[t2] != eachRoundSBoxDiff1_2[t2] || test2[t2] != eachRoundSBoxDiff2_2[t2]) {
				iscorrect_2 = false;
				break;
			}
		}

		if (iscorrect_2) {
			bool diff1_1[21], diff1_2[21];
			bool iscorrect_1 = true;
			for (int i = 0; i < solNum; i++) {
				for (int k = 0; k < bs; k++) {
					diff1_1[k] = expression.ma[k][expression.c - 1];
					diff1_2[k] = expression.ma[bs + k][expression.c - 1];
					for (int ti = 0; ti < expression.c - 1; ti++) {
						if (expression.ma[k][ti]) {
							diff1_1[k] = diff1_1[k] ^ sol[i][ti];
						}
						if (expression.ma[bs + k][ti]) {
							diff1_2[k] = diff1_2[k] ^ sol[i][ti];
						}
					}
				}
				iscorrect_1 = true;
				for (int k = 0; k < bs; k++) {
					if (diff1_1[k] != eachRoundDiff1_1[k] || diff1_2[k] != eachRoundDiff2_1[k]) {
						iscorrect_1 = false;
						break;
					}
				}
				if (iscorrect_1) {
					
					cout << "The find1:";
					for (int i = 0; i < bs; i++) {
						cout << diff1_1[i];
					}
					cout << " ";
					for (int i = 0; i < bs; i++) {
						cout << diff1_2[i];
					}
					cout << endl;
					cout << "The corr1:";
					for (int i = 0; i < bs; i++) {
						cout << eachRoundDiff1_1[i];
					}
					cout << " ";
					for (int i = 0; i < bs; i++) {
						cout << eachRoundDiff2_1[i];
					}
					cout << endl;

					cout << "The findSbox2:";
					for (int i = 0; i < bs; i++) {
						cout << test1[i];
					}
					cout << " ";
					for (int i = 0; i < bs; i++) {
						cout << test2[i];
					}
					cout << endl;
					cout << "The corrSbox1:";
					for (int i = 0; i < bs; i++) {
						cout << eachRoundSBoxDiff1_2[i];
					}
					cout << " ";
					for (int i = 0; i < bs; i++) {
						cout << eachRoundSBoxDiff2_2[i];
					}
					cout << endl;

					cout << "The correct one is being enumerated!" << endl;
				}

			}

		}	
	}

	for (int i = 0; i < solNum; i++) {
		sol[i].clear();
	}
	sol.clear();

	return solNum;
}

void LowMC::gauss(matrix& eqSys) {
	int variableNum = eqSys.c - 1;
	bool isFirst = false;
	int targetRow = 0;

	for (int i = 0; i < variableNum; i++) {
		isFirst = true;
		for (int j = targetRow; j < eqSys.r; j++) {
			if (isFirst && eqSys.ma[j][i]) {
				isFirst = false;
				swap(eqSys.ma[j], eqSys.ma[targetRow]);
				targetRow++;
			}
			else {
				if (eqSys.ma[j][i]) {//apply Gauss
					for (int k = i; k < eqSys.c; k++) {
						eqSys.ma[j][k] ^= eqSys.ma[targetRow - 1][k];
					}
				}
			}
		}
	}
}

void LowMC::storeSolutions(vector<vector<bool> >& sol, matrix& eqSys, int& solNum) {
	vector<int> lead;
	vector<int> freebits;
	freebits.clear();
	lead.clear();
	bool* isFree;
	isFree = new bool[eqSys.c - 1];
	memset(isFree, 1, eqSys.c - 1);

	int start = 0;
	for (int r = 0; r < eqSys.r; r++) {
		while (start < eqSys.c - 1 && eqSys.ma[r][start] == 0) {
			start++;
		}
		if (start == eqSys.c - 1) {
			break;
		}
		lead.push_back(start);
		isFree[start] = false;
		start++;
	}

	if (lead.size() < eqSys.r) {
		for (int j = lead.size(); j < eqSys.r; j++) {
			if (eqSys.ma[j][eqSys.c - 1] != 0) {
				solNum = 0;
				return;
			}
		}
	}

	for (int i = 0; i < eqSys.c - 1; i++) {
		if (isFree[i]) {
			freebits.push_back(i);
		}
	}
	//cout << "free size:" << freebits.size() << endl;
	//cout << "lead size:" << lead.size() << endl;*/

	vector<bool> eachsol;
	eachsol.clear();
	eachsol.resize(eqSys.c - 1);
	int solSize = 1 << freebits.size();
	for (int i = 0; i < solSize; i++) {
		for (int j = 0; j < freebits.size(); j++) {
			eachsol[freebits[j]] = (i >> j) & 0x1;
		}
		for (int k = lead.size() - 1; k >= 0; k--) {
			//compute eachsol[lead[k]] use row= k
			eachsol[lead[k]] = eqSys.ma[k][eqSys.c - 1];
			for (int j = lead[k] + 1; j < eqSys.c - 1; j++) {
				if (eqSys.ma[k][j] == 1) {
					eachsol[lead[k]] = eachsol[lead[k]] ^ eachsol[j];
				}
			}
		}
		solNum++;
		sol.push_back(eachsol);
	}

	delete[]isFree;
	freebits.clear();
	lead.clear();
	eachsol.clear();
}
