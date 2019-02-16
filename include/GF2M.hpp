/*
A class to represent the extension field GF(2^M) 
*/

#ifndef EXTENSION_FIELD2
#define EXTENSION_FIELD2

#include <boost/operators.hpp>
#include <iostream>
#include <ostream>
#include <vector>
#include <utility>
#include "./helpers.hpp"

using namespace std;



// divide a by b; computer a/b and remainder(a/b)
// return pair with fist: quotient, second: remainder
template<class uint>
pair<uint,uint> divide(const uint a, const uint b){

	uint quot = 0; // quotient
	uint rem = a; // remainder 
	
	if(b > a){ // then a is not divisible by b, and quot = 0
		return pair<uint,uint>(quot,rem);
	}
	
	while((rem >= b) && (rem != 0)){
		unsigned exp = degree(rem) - degree(b);
		quot |= (1 << exp);
		uint btmp = (b << exp);
		rem ^= btmp;
	}
	return pair<uint,uint>(quot,rem);
};

/* 
GF(2^m)
Parameters are:
- m: the degree of the extension
- prim_poly: primitive polynomial defining the extension field
*/
template<class uint_, uint_ m_, uint_ prim_poly_>
class GF2M:
	boost::field_operators< GF2M<uint_,m_,prim_poly_>,
	boost::equality_comparable< GF2M<uint_,m_,prim_poly_>
	> >
{
public: 
	//! an element of the extension field is a polynomial, represented by an unsigned integer
	typedef uint_ uint;
	static const unsigned m = m_;
	static const uint prim_poly = prim_poly_;
	uint el; 
	bool empty;
	//! constructors
	GF2M(uint x){el = x % 2;empty = false;};
	//GF2M(unsigned x){el = x ;empty = false;};
	GF2M(uint x,unsigned v){el = x ;empty = false;}; // generates instance with el = x
	
	GF2M(){ el = 0; empty = true;}; // -1 marks empty 


	GF2M& operator += (const GF2M& x){
		el = el ^ x.el;
		return *this;
	}
	
	GF2M& operator -= (const GF2M& x){
		el = el ^ x.el;
		return *this;
	}

	GF2M& operator *= (const GF2M& x){
		uint a = el;
		uint b = x.el;
		el = 0; // product of this and x
		while (a && b) {
				if (b & 1) // if b is odd, add the corresponding a to p 
					el ^= a;

				if (a &  ( ((uint)1) << (m-1)) ) // GF modulo: if a >= 2^m, it overflows when shifted left, so reduce 
					a = (a << 1) ^ prim_poly; // XOR with the primitive polynomial
				else
					a <<= 1; // a*2 /* equivalent to a*2 
				b >>= 1; // b / 2	
		}
		return *this;
	}
		
	GF2M& operator /= (const GF2M& x){
		*this = (x.inverse() * (*this)); 
		return *this;
	}
	
	//! compute the multiplicative inverse via the extended Euclidean algorithm 
	GF2M inverse() const {
		
		uint rem1 = prim_poly;
		uint rem2 = el;
		uint aux1 = 0;
		uint aux2 = 1;

		while( rem2 != 0 ){
			// res.first = quotient(rem1/rem2)
			// res.second = remainder(rem1/rem2)
			
			pair< uint, uint > res = divide(rem1,rem2);
			
			GF2M aux_new = GF2M(res.first,0)*GF2M(aux2,0) + GF2M(aux1,0);

			// prepare for the next step
			rem1 = rem2;
			rem2 = res.second;
			aux1 = aux2;
			aux2 = aux_new.el;
			
		}
		//aux1.multiply_monomial(pow(rem1.poly[0],-1), 0);
		return GF2M(aux1,0);
	}


	unsigned order() const {
		// multiplicative neutral element = 1
		//PFE one = PFE(1);
		GF2M one = GF2M(1); 
		GF2M tmp = *this;  //EFE(el);		
		unsigned ord = 1;
		while(!(tmp == one)){
			ord++;
			tmp *= *this;
		}
		return ord;
	}

	bool iszero() const {
		return (el == 0);
	}
	
	bool isempty() const {
		return empty;
	}

	bool operator ==(const GF2M& x) const {
		return (el == x.el);
	}
	
	bool operator < (const GF2M& x) const {
		return (el < x.el);	
	}


};	

template<class uint, uint m, uint prim_poly> 
ostream &operator<<(ostream& stream, const GF2M<uint,m,prim_poly>& x)
    {
	  	stream << bitset<32>(x.el);
      	return stream; // must return stream
	};



template<class uint, uint m, uint prim_poly> 
GF2M<uint,m,prim_poly> pow(const GF2M<uint,m,prim_poly>& a, int exp){
	
	if(a.iszero()) return a;
	if(exp == 0) {return GF2M<uint,m,prim_poly>(1);}
	
	GF2M<uint,m,prim_poly> tmp;
	
	if(exp < 0){
		GF2M<uint,m,prim_poly> am = a.inverse(); 
		tmp = am;
		for(unsigned i =1; i < abs(exp); ++i) tmp *= am;
	} else { // exp > 0
		tmp = a;
		for(unsigned i =1; i < abs(exp); ++i) tmp *= a;
	}

	return tmp;
};



#endif
