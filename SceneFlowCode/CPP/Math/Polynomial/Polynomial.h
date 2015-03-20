#ifndef POLYNOMIAL_INCLUDED
#define POLYNOMIAL_INCLUDED

#include <vector>

template<int Degree>
class Polynomial{
public:
	double coefficients[Degree+1];

	Polynomial(void);
	template<int Degree2>
	Polynomial(const Polynomial<Degree2>& P);
	double operator()(const double& t) const;
	double integral(const double& tMin,const double& tMax) const;

	int operator == (const Polynomial& p) const;
	int operator != (const Polynomial& p) const;
	int isZero(void) const;
	void setZero(void);

	template<int Degree2>
	Polynomial& operator  = (const Polynomial<Degree2> &p);
	Polynomial& operator += (const Polynomial& p);
	Polynomial& operator -= (const Polynomial& p);
	Polynomial  operator -  (void) const;
	Polynomial  operator +  (const Polynomial& p) const;
	Polynomial  operator -  (const Polynomial& p) const;
	template<int Degree2>
	Polynomial<Degree+Degree2>  operator *  (const Polynomial<Degree2>& p) const;

	Polynomial& operator += (const double& s);
	Polynomial& operator -= (const double& s);
	Polynomial& operator *= (const double& s);
	Polynomial& operator /= (const double& s);
	Polynomial  operator +  (const double& s) const;
	Polynomial  operator -  (const double& s) const;
	Polynomial  operator *  (const double& s) const;
	Polynomial  operator /  (const double& s) const;

	Polynomial scale(const double& s) const;
	Polynomial shift(const double& t) const;

	Polynomial<Degree-1> derivative(void) const;
	Polynomial<Degree+1> integral(void) const;

	void printnl(void) const;

	Polynomial& addScaled(const Polynomial& p,const double& scale);

	static void Negate(const Polynomial& in,Polynomial& out);
	static void Subtract(const Polynomial& p1,const Polynomial& p2,Polynomial& q);
	static void Scale(const Polynomial& p,const double& w,Polynomial& q);
	static void AddScaled(const Polynomial& p1,const double& w1,const Polynomial& p2,const double& w2,Polynomial& q);
	static void AddScaled(const Polynomial& p1,const Polynomial& p2,const double& w2,Polynomial& q);
	static void AddScaled(const Polynomial& p1,const double& w1,const Polynomial& p2,Polynomial& q);

	void getSolutions(const double& c,std::vector<double>& roots,const double& EPS) const;

private:

	int Factor(double a1,double a0,double roots[1][2],const double& EPS);
	int Factor(double a2,double a1,double a0,double roots[2][2],const double& EPS);
	int Factor(double a3,double a2,double a1,double a0,double roots[3][2],const double& EPS);
	int Factor(double a4,double a3,double a2,double a1,double a0,double roots[4][2],const double& EPS);
};

#include "Polynomial.inl"
#endif // POLYNOMIAL_INCLUDED
