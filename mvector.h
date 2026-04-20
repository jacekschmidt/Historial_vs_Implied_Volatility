#ifndef MVECTOR_H // the 'include guard'
#define MVECTOR_H // see C++ Primer Sec. 2.9.2

#include <vector>
#include <iostream>

// Class that represents a mathematical vector
class MVector
{
public:
	// constructors
	MVector() {}
	explicit MVector(int n) : v(n) {}
	MVector(int n, double x) : v(n, x) {}
	MVector(std::initializer_list<double> l) : v(l) {}

	// access element (lvalue) (see example sheet 5, q5.6)
	double &operator[](int index) 
	{ 
		return v[index];
	}

	// access element (rvalue) (see example sheet 5, q5.7)
	double operator[](int index) const {
		return v[index]; 
	}

	int size() const { return v.size(); } // number of elements

	// add data to your vector
	void push_back(double x){v.push_back(x);}


private:
	std::vector<double> v;
};

//1.2.1

inline MVector operator*(const double& lhs, const MVector& rhs)
{
	MVector temp = rhs;
	for (int i=0; i<temp.size(); i++) temp[i]*=lhs;
	return temp;
}

inline MVector operator*(const MVector& rhs, const double& lhs)
{
	MVector temp = rhs;
	for (int i=0; i<temp.size(); i++) temp[i]*=lhs;
	return temp;
}

inline MVector operator/(const MVector& rhs, const double& lhs)
{
	MVector temp = rhs;
	for (int i=0; i<temp.size(); i++) temp[i]/=lhs;
	return temp;
}

inline MVector operator+(const MVector& lhs, const MVector& rhs)
{
	MVector temp = rhs;
	for (int i=0; i<temp.size(); i++) temp[i]=lhs[i]+rhs[i];
	return temp;
}

inline MVector operator-(const MVector& lhs, const MVector& rhs)
{
	MVector temp = rhs;
	for (int i=0; i<temp.size(); i++) temp[i]=lhs[i]-rhs[i];
	return temp;
}

std::ostream& operator<<(std::ostream& os, const MVector& v)
{
	int n = v.size();
	os << "("<<v[0];
	for(int i=1; i<n;i++) os<<","<<v[i];
	os << ")";
	return os;
}


#endif
