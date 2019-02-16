/*
A class to represent an extension field 
*/

#ifndef EXTENSION_FIELD
#define EXTENSION_FIELD

#include <boost/operators.hpp>
#include <iostream>
#include <ostream>
#include <vector>

#include "polynomial.hpp"

using namespace std;

// extension field element (EFE)
template<class PFE> class EFE:
	boost::field_operators< EFE<PFE>,
	boost::equality_comparable< EFE<PFE>
	> >
{
public: 
	//! an element of the extension field is a polynomial with coefficients in the prime field PFE
	polynomial<PFE> el;
	//! the degree of the extension
	static unsigned m; 
	//! primitive polynomial defining the extension field
	static polynomial<PFE> prim_poly;
	//! constructors
	EFE(vector<PFE>& v): el(polynomial<PFE>(v)) {};
	EFE(polynomial<PFE>& el_): el(el_) {};
	EFE(unsigned x){ // set the coeff with exponent zero to x
		vector<PFE> tmpv(m, PFE(0) );
		tmpv[0] = PFE(x);
		el = polynomial<PFE>(tmpv);
	}
	EFE(){};

	EFE& operator += (const EFE& x){
		el += x.el;
		return *this;
	}
	
	EFE& operator -= (const EFE& x){
		el -= x.el;
		return *this;
	}

	EFE& operator *= (const EFE& x){
		// multiply
		polynomial<PFE> ptmp = el*x.el;
		// modulo the primitive polynomial
	
		int exp;
		while( (exp = ptmp.degree() - prim_poly.degree()) >= 0 ){
			// difference in exponents
			//unsigned exp = ptmp.degree() - prim_poly.degree();
			// primitive polynomial; to be muliplied..	
			polynomial<PFE> primtmp = prim_poly;
			PFE co_mo = ptmp.poly[ptmp.degree()];  
			co_mo /= prim_poly.poly[prim_poly.degree()];
			primtmp.multiply_monomial( co_mo , exp  );
			ptmp -= primtmp;
		}

		ptmp.poly.resize(ptmp.degree()+1);
		el = ptmp;
		return *this;
	}

	EFE& operator /= (const EFE& x){
		*this = (x.inverse() * (*this)); 
		return *this;
	}
	
	//! compute the multiplicative inverse via the extended Euclidean algorithm 
	EFE inverse() const {
		
		polynomial<PFE> rem1 = prim_poly;
		polynomial<PFE> rem2 = el;
		polynomial<PFE> aux1 = polynomial<PFE>(0);
		polynomial<PFE> aux2 = polynomial<PFE>(1);


		while(! rem2.iszero() ){
			// res.first = quotient(rem1/rem2)
			// res.second = remainder(rem1/rem2)
			
			pair<polynomial<PFE>, polynomial<PFE> > res = divide<PFE>(rem1,rem2);

			polynomial<PFE> aux_new = polynomial<PFE>(0) - res.first*aux2 + aux1;

			// prepare for the next step
			rem1 = rem2;
			rem2 = res.second;
			aux1 = aux2;
			aux2 = aux_new;
			
		}
	
		aux1.multiply_monomial(pow(rem1.poly[0],-1), 0);

		return aux1;
	}


	
	unsigned order() const {
		// multiplicative neutral element = 1
		//PFE one = PFE(1);
		vector<PFE> vecone = vector<PFE>(1,PFE(1));
		EFE one = EFE(vecone); 
		EFE tmp = *this;  //EFE(el);		
		unsigned ord = 1;
		while(!(tmp == one)){
			ord++;
			tmp *= *this;
			//cout << tmp << endl;
		}
		return ord;
	}
	/*
	// not a good solution
	vector<unsigned> tovec() const {
		vector<unsigned> vec(m,0);
		for(unsigned i=0;i<el.poly.size();++i){
			vec[i] = el.poly[i].element;
		}
		return vec;
	}
	*/

	vector<PFE> tovec() const {
		if(el.poly.size()!=m){
			vector<PFE> tmp(m,PFE(0));
			for(unsigned i=0;i<el.poly.size();++i){
				tmp[i] = el.poly[i];
			}
			return tmp;
		}else
			return el.poly;
	}

	bool iszero() const {
		vector<PFE> veczero = vector<PFE>(1,PFE(0));
		EFE zero = EFE(veczero);
		return (zero == *this);
	}
	
	bool isempty() const {
		return (el.poly.size() == 0);
	}

	bool operator ==(const EFE& x) const {
		
		for(unsigned i = 0; i < min(el.poly.size(),x.el.poly.size()); ++i)
			if(!(x.el.poly[i] == el.poly[i]) ) return false;
		
		// need to verify if the coefficients of the larger vector that are necessarily zero in the smaller vector are also zero in the larger vector 
		if(el.poly.size() > x.el.poly.size()){
			for(unsigned i = x.el.poly.size() ; i<el.poly.size(); ++i) 
				if( !(el.poly[i].iszero()) ) return false;
		} else if(el.poly.size() < x.el.poly.size()){
			for(unsigned i = el.poly.size(); i<x.el.poly.size(); ++i) 
				if( !(x.el.poly[i].iszero()) ) return false;
		}

		return true;
	}
};	

 template<class PFE>
    ostream &operator<<(ostream &stream, const EFE<PFE> x)
    {
	  	stream <<  x.el;
      	return stream; // must return stream
	};

template<class PFE>
 	EFE<PFE> pow(const EFE<PFE>& a, int exp){
		
		if(a.iszero()) return a;
		if(exp == 0) {
			vector<PFE> vecone = vector<PFE>(1,PFE(1));
			return EFE<PFE>(vecone); 
		}
		
		EFE<PFE> tmp;
		
		if(exp < 0){
			EFE<PFE> am = a.inverse(); 
			tmp = am;
			for(unsigned i =1; i < abs(exp); ++i) tmp *= am;
		} else { // exp > 0
			tmp = a;
			for(unsigned i =1; i < abs(exp); ++i) tmp *= a;
		}

		return tmp;
	};



#endif
