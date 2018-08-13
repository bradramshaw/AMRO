#include "stdafx.h"
#include "DataExtractor.h"
#include "FindFermi.h"
using namespace std;

//int FindFermi::func(Ipp64f * params, Ipp64f * argkz, Ipp64f * argCos, Ipp64f * argSin, Ipp64f * r, int length, Ipp64f * temp, Ipp64f * out)
int FindFermi::func(Ipp64f *params, Ipp64f * argkz, Ipp64f * argCos, Ipp64f * argSin, Ipp64f *r, int length, Ipp64f *temp, Ipp64f *out) {

	ippsMulC_64f(argCos, 3.747665940, temp, length);
	ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[1 * length]); // cos cos
	ippsMulC_64f(argSin, 3.747665940, temp, length);
	ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[2 * length]); // cos sin
	ippsMulC_64f(argCos, 3.747665940 / 2, temp, length);
	ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[3 * length]); // cos cos/2
	ippsMulC_64f(argSin, 3.747665940 / 2, temp, length);
	ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[4 * length]); // cos sin/2
	ippsMulC_64f(argCos, 3.747665940 * 2, temp, length);
	ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[5 * length]); // cos 2 cos
	ippsMulC_64f(argSin, 3.747665940 * 2, temp, length);
	ippsMul_64f_I(r, temp, length);
	vdCos(length, temp, &temp[6 * length]); // cos 2 sin
	vdCos(length, argkz, &temp[7 * length]); // cos kz/2
	ippsMulC_64f(argkz, 2, temp, length);
	vdCos(length, argkz, &temp[8 * length]); // cos kz

	ippsAdd_64f(&temp[5 * length], &temp[6 * length], temp, length);// param 5
	ippsMulC_64f(temp, -35164.83516*params[5 - 1], out, length);

	ippsMul_64f(&temp[1 * length], &temp[2 * length], temp, length);// param 4
	ippsMulC_64f_I(-35164.83516 * 2 * params[4 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);

	ippsAdd_64f(&temp[1 * length], &temp[2 * length], temp, length);// param 3
	ippsMulC_64f_I(-35164.83516 * params[3 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);
	ippsAddC_64f_I(35164.83516 / 2 * params[2 - 1], out, length);// param 2
	ippsSub_64f(&temp[2 * length], &temp[1 * length], temp, length);// param 6
	ippsSqr_64f_I(temp, length); //square
	ippsMul_64f_I(&temp[3 * length], temp, length); // mult by cos cos/2
	ippsMul_64f_I(&temp[4 * length], temp, length); // mult by cos sin/2
	ippsMul_64f_I(&temp[7 * length], temp, length); // mult by cos  kz/2
	ippsMulC_64f_I(-35164.83516 * params[6 - 1], temp, length);
	ippsMulC_64f_I(-35164.83516 *params[7 - 1], &temp[8 * length], length);
	ippsAdd_64f_I(temp, out, length);
	ippsAdd_64f_I(&temp[8 * length], out, length);

	return 0;
}


 FindFermi:: FindFermi(std::string name, Ipp64f * param)
{
	DataExtractor extractor(name);
	//Ipp64f params[9] = { 0.074, 475, 525, -60, 16, 1000, 0.5, 17, 8 };
	UpdatePar(param);
	Ipp64f *start = extractor.getDataArray(); 
	nPoints = floor((extractor.getNumberOfLines()) / 2);
	tol = 1E-10;
	subMaxR = new Ipp64f[1];
	temp1 = new Ipp64f[20*nPoints];
	argCos = new Ipp64f[nPoints];
	argSin = new Ipp64f[nPoints];
	argkz = new Ipp64f[nPoints];
	theta = new Ipp64f[nPoints];
	kz = new Ipp64f[nPoints];
	funcval = new Ipp64f[nPoints];
	dfunc1 = new Ipp64f[nPoints];
	dfunc2 = new Ipp64f[nPoints];
	dfunc = new Ipp64f[nPoints];
	r = new Ipp64f[nPoints];
	rtemp = new Ipp64f[nPoints];
	kx = new Ipp64f[nPoints];
	ky = new Ipp64f[nPoints];
	//Ipp64f *startpoint = new Ipp64f[3 * nPoints];
	for (int i = 0; i < nPoints; i++)
	{
		/*cout << starts[2 * i] << " " << starts[2 * i + 1] << endl;*/
		theta[i] = start[2 * i];
		kz[i] = start[2 * i + 1];
		r[i] = 0.3;
	}
	vdSin(nPoints, theta, argSin);
	vdCos(nPoints, theta, argCos);
	ippsMulC_64f(kz, 6.6, argkz, nPoints);
	cout << "parameter: ";
	for (int i = 0; i < 10; ++i) {
		cout << params[i] << ", ";
	}
	cout << endl;
}

int FindFermi::UpdatePar(Ipp64f * param)
{
	for (int i = 0; i < 10; ++i) {
		params[i] = param[i];
	}

	return 0;
}

int FindFermi::PrintPar()
{
	cout << "parameter: ";
	for (int i = 0; i < 10; ++i) {
		cout << params[i] <<", ";
	}
	cout << endl;
	return 0;
}



int FindFermi::ReturnStart(Ipp64f * startpoint)
{
	for (int i = 0; i < 50; i++)
	{
		/*vdSin(nPoints, theta, argSin);
		vdCos(nPoints, theta, argCos);
		ippsMulC_64f_I(3.74767, argSin, nPoints);
		ippsMulC_64f_I(3.74767, argCos, nPoints);
		ippsMulC_64f(kz, 3.3, argkz, nPoints);
		func(argkz, argCos, argSin, r, nPoints, temp1, temp2, temp3, funcval);

		ippsAddC_64f(r, 1E-5, rtemp, nPoints);
		func(argkz, argCos, argSin, rtemp, nPoints, temp1, temp2, temp3, dfunc1);
		ippsAddC_64f(r, -1E-5, rtemp, nPoints);
		func(argkz, argCos, argSin, rtemp, nPoints, temp1, temp2, temp3, dfunc2);
		ippsMulC_64f_I(-1, dfunc2, nPoints);
		ippsAdd_64f(dfunc1, dfunc2, dfunc, nPoints);
		ippsDivC_64f_I(2E-5, dfunc, nPoints);

		ippsDiv_64f_I(dfunc, funcval, nPoints);
		ippsSub_64f_I(funcval, r, nPoints);*/

		//caxis lattice parameter =3.3  k_z * c
		func(params, argkz, argCos, argSin, r, nPoints, temp1, funcval);
		//func(nPoints, temp1,funcval);
		//cout << argCos[1] << endl;
		ippsAbs_64f(funcval, temp1, nPoints);
		ippsMax_64f(temp1, nPoints, subMaxR);
		//cout << i << endl;
		//cout << funcval[0] << endl;
		//cout << subMaxR[0] << endl;
		if (subMaxR[0] <= tol)
		{

			break;
		}
		ippsAddC_64f(r, 1E-6, rtemp, nPoints);
		func(params, argkz, argCos, argSin, rtemp, nPoints, temp1, dfunc1);
		//func(nPoints, temp1,dfunc1);
		ippsAddC_64f(r, -1E-6, rtemp, nPoints);
		func(params, argkz, argCos, argSin, rtemp, nPoints, temp1, dfunc2);
		//func(nPoints, temp1, dfunc2);
		ippsMulC_64f_I(-1, dfunc2, nPoints);
		ippsAdd_64f(dfunc1, dfunc2, dfunc, nPoints);
		ippsDivC_64f_I(2E-6, dfunc, nPoints);

		ippsDiv_64f_I(dfunc, funcval, nPoints);
		ippsSub_64f_I(funcval, r, nPoints);

	}

	//for (int i = 0; i < nPoints; i++)
	//{
	//	cout << r[i] << endl;
	//}

	ippsMul_64f(r, argCos, kx, nPoints);
	ippsMul_64f(r, argSin, ky, nPoints);
	for (int i = 0; i < nPoints; i++) {
		startpoint[i * 3] = kx[i];
		startpoint[i * 3 + 1] = ky[i];
		startpoint[i * 3 + 2] = kz[i];
	}
	ofstream fout;
	fout.open("FindFermi.dat");
	fout.precision(15);

	for (int i = 0; i < nPoints; i++) {

		//fout << kx[i] << "\t" << ky[i] << "\t" << kz[i] << endl;
		fout << theta[i]<<"\t"<<r[i]<< "\t" <<kz[i]<< endl;
		//cout << kx[i] << "\t" << ky[i] << "\t" << kz[i] << endl;
	}

	fout.close();
	//while (true);
	return 0;
}

int FindFermi::ReturnNumPoint()
{
	cout << "nPoints = " << &nPoints <<"  "<<nPoints<< endl;
	return nPoints;
}


FindFermi::~FindFermi()
{
	delete subMaxR;
	delete temp1;
	delete argCos;
	delete argSin;
	delete theta;
	delete argkz;
	delete kz;
	delete funcval;
	delete dfunc1;
	delete dfunc2;
	delete r;
	delete rtemp;
	delete kx;
	delete ky;

}
