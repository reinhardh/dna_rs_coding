/*
A class to represent polynomials
*/


#ifndef POLYNOMIAL
#define POLYNOMIAL

#include <vector>

template<class coeff_t>
class polynomial:
	boost::ring_operators< polynomial<coeff_t>,
    boost::equality_comparable< polynomial<coeff_t>
	> >
{
public:
	//! contains the coefficients of the polynomial
	std::vector<coeff_t> poly;
	//! constructor
	polynomial(const std::vector<coeff_t>& p): poly(p) {};
	polynomial(coeff_t c){
		poly = std::vector<coeff_t>(1,c);
	};
	polynomial(){};

	//! Computes the degree of the polynomial, returns -1 if the polynomial is equal to 0.
	unsigned degree() const 
	{
		unsigned deg = 0;
		for(unsigned i=0;i!=poly.size();++i)
			if( !( poly[i].iszero()) )  deg = i;
		return deg;	
	}

	polynomial& operator += (const polynomial& x){
		for(unsigned i = 0;i != std::min(poly.size(),x.poly.size()); i++)
			poly[i] += x.poly[i];
		// account for differnt length
		// if(poly.size() > x.poly.size()) // nothing to do 
		if(poly.size() < x.poly.size()){
			unsigned i = poly.size();
            poly.resize(x.poly.size()); 
			for(; i<x.poly.size(); ++i)
				poly[i] = x.poly[i];
		}
		return *this;
	}
	
	polynomial& operator -= (const polynomial& x){
		for(unsigned i = 0;i != std::min(poly.size(),x.poly.size()); ++i)
			poly[i] -= x.poly[i];
		// account for differnt length
		// if(poly.size() > x.poly.size()) // nothing to do 
		if(poly.size() < x.poly.size()){
			unsigned i = poly.size();
            poly.resize(x.poly.size()); 
			for(; i<x.poly.size(); ++i)
				poly[i] = coeff_t(0) - x.poly[i];
		}
		return *this;
	}
	

	polynomial& operator *= (const polynomial& x){
		std::vector<coeff_t> tmp(x.poly.size()+poly.size()-1,0);
		for(unsigned i=0;i!=poly.size();++i){
			//multiply the ith element of poly with x.poly
			for(unsigned j=0;j!=x.poly.size();++j){
				tmp[i+j] += poly[i] * x.poly[j];
			}
		}
		
		// ``shorten'' the vector to its degree.. 
		unsigned deg = 0;
		for(unsigned i=0;i!=tmp.size();++i)
			if( !( tmp[i].iszero() ) ) deg = i;
		deg += 1;
		
		tmp.resize(deg);
		poly = tmp;

		return (*this);
	}




	//! multipy the polynomial with a mononmial with coefficient c and exponent exp
	void multiply_monomial(coeff_t c,unsigned exp){
		std::vector<coeff_t> tmp(poly.size() + exp);
		for(unsigned i=0;i!=poly.size();++i)
			tmp[i+exp] = c * poly[i];
		
		// ``shorten'' the vector to its degree.. 
		unsigned deg = 0;
		for(unsigned i=0;i!=tmp.size();++i)
			if( !( tmp[i].iszero() ) ) deg = i;
		deg += 1;
		tmp.resize(deg);
		poly = tmp;
	}
	
	void derive(){
		assert(poly.size()>0);
		poly[0] = coeff_t(0);
		for(unsigned i=0;i < poly.size()-1; ++i ){
			poly[i] = coeff_t(i+1)*poly[i+1];
		}
		poly.resize(poly.size()-1);
	}

	// evaluate the polynomial at x 
	coeff_t evaluate(coeff_t& x){
		coeff_t sum = coeff_t(0);
		for(unsigned i=0; i<poly.size();++i){
			sum += poly[i]*pow(x,i);
		}
		return sum;
	}

    bool iszero() const {
        polynomial pzero = polynomial(0);
        return (pzero == *this);
    }	

	bool operator ==(const polynomial& x){
        for(unsigned i = 0; i < std::min(poly.size(),x.poly.size()); ++i)
            if(!(x.poly[i] == poly[i]) ) return false;

        // need to verify if the coefficients of the larger vector that are necessarily zero in the smaller vector are also zero in the larger vector 
        if(poly.size() > x.poly.size()){
            for(unsigned i = x.poly.size() ; i<poly.size(); ++i)
                if( !(poly[i].iszero()) ) return false;
        } else if(poly.size() < x.poly.size()){
            for(unsigned i = poly.size(); i<x.poly.size(); ++i)
                if( !(x.poly[i].iszero()) ) return false;
        }

        return true;
    }

};

 template<class coeff_t>
    std::ostream &operator<<(std::ostream &stream, polynomial<coeff_t>  x)
    {
      	
		if(x.poly.size() == 0){
			stream << "empty polynomial";
			return stream;
		}
		 
		for(int i = x.poly.size()-1; i > 0  ; i--){
        	//stream <<  x.poly[i] << "x" << i << " + ";
        	stream <<  x.poly[i] << " ";
		}
	  	
        stream <<  x.poly[0];

      	return stream; // must return stream
    };



// divide a by b; computer a/b and remainder(a/b)
// return pair with fist: quotient, second: remainder
template<class co_t>
std::pair<polynomial<co_t>,polynomial<co_t> > divide(const polynomial<co_t>& a, const polynomial<co_t>& b){

	polynomial<co_t> quot; // quotient
	polynomial<co_t> rem = a; // remainder 
	
	if(b.degree() > a.degree()){ // then a is not divisible by b, and quot = 0
		quot.poly.resize(1);
		quot.poly[0] = co_t(0);
		return std::pair<polynomial<co_t>,polynomial<co_t> >(quot,rem);
	}
	
	quot.poly.resize(a.degree()-b.degree() + 1); 
	
	while((rem.degree() >= b.degree()) && !(rem.iszero())    ){
		unsigned exp = rem.degree() - b.degree();
		quot.poly[exp] = rem.poly[rem.degree()] / b.poly[b.degree()] ; //a_max / b_max
		polynomial<co_t> btmp = b;
		btmp.multiply_monomial(quot.poly[exp],exp); //
		rem -= btmp;
	}
	rem.poly.resize(rem.degree()+1);
	return std::pair<polynomial<co_t>,polynomial<co_t> >(quot,rem);

};





#endif
