#pragma once
class FindFermi
{
private:
	int func(Ipp64f *params, Ipp64f * argkz, Ipp64f * argCos, Ipp64f * argSin, Ipp64f *r, int length, Ipp64f *temp, Ipp64f *out);
	//int nPoints;
	Ipp64f *subMaxR ;
	Ipp64f *temp1;
	Ipp64f *argCos ;
	Ipp64f *argSin ;
	Ipp64f *argkz;
	Ipp64f *theta ;
	Ipp64f *kz ;
	Ipp64f *funcval;
	Ipp64f *dfunc1 ;
	Ipp64f *dfunc2;
	Ipp64f *dfunc;
	Ipp64f *r;
	Ipp64f *rtemp ;
	Ipp64f *kx ;
	Ipp64f *ky ;
	//Ipp64f *startpoint ;
	Ipp64f tol;
	Ipp64f params[9];
public:
	
	int nPoints;
	FindFermi(std::string name, Ipp64f * param);
	int UpdatePar(Ipp64f *param);
	int PrintPar();
	int ReturnStart(Ipp64f *startpoint);
	int ReturnNumPoint();
	
	virtual ~FindFermi();
};

