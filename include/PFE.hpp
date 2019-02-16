/*
A class to represent a prime field 
*/


#ifndef PRIME_FIELD
#define PRIME_FIELD

#include <boost/operators.hpp>
#include <iostream>
#include <ostream>
#include <boost/bimap.hpp>

using namespace std;




template< class MapType >
void print_map(const MapType & map,
               const std::string & separator,
               std::ostream & os )
{
    typedef typename MapType::const_iterator const_iterator;

    for( const_iterator i = map.begin(), iend = map.end(); i != iend; ++i )
    {
        os << i->first << separator << i->second << std::endl;
    }
}

// prime field element (PFE)
template<class el_t = int> class PFE:
	boost::field_operators< PFE<el_t>,
	boost::equality_comparable< PFE<el_t>
	> >
{
public: 
	el_t element;
	static el_t prime;

	typedef boost::bimap< el_t, el_t> bim; 
	//! maps exponents to elements 
	static bim exp_el;
	
	//!Constructs a prime f
	PFE(el_t el): element(el % prime) {};
	PFE(){element = 0;};

	static void set_prime(const el_t& a){prime = a;}

	//! for initializying the map exponent to element; takes as parameter a primitive element of GF(prime)
	
	static void initialize_exp_el(el_t priel){
		
		int el = 1;
		for(unsigned i = 0; i < prime-1; ++i){
			// insert (exponent, corresponding element)
			
			exp_el.insert(typename bim::value_type(i,el ));
			el = (el*priel) % prime;
		//cout << i << "  " <<  ((int) pow((double) priel,(int) i)) % prime << endl;
		}
		
	}
	
	bool isempty() const {return (element==prime);}
	void setempty() {element = prime;};



	PFE& operator += (const PFE& x){
		element = (element+x.element) % prime;
		return *this;
	}
	
	PFE& operator -= (const PFE& x){
		element = (element-x.element) % prime;
		if(element < 0) element += prime;
		return *this;
	}

	PFE& operator *= (const PFE& x){
		
		element = element*x.element % prime;
		return *this;
	}


	
	PFE& operator /= (const PFE& x){
		//print_map(exp_el.left , " " , cout);
		el_t ea = exp_el.right.at(element);	
		el_t eb = exp_el.right.at(x.element);
		el_t eres = (prime - 1 - eb + ea) % (prime-1);
		element = exp_el.left.at(eres);
		return *this;
	}

	//operator const el_t (){
	//	return element;
	//}
	//operator el_t (){return element;}
	operator int (){return (int) element;};


	PFE inverse() const {
		return pow(*this, -1);
	}

	bool operator ==(const PFE& x) const {
		return (x.element == element);
	}
	bool operator < (const PFE& x) const {
		return (element < x.element);	
	}
	
	bool iszero() const { return (element == 0);  };

};	

 template<class el_t>
    ostream &operator<<(ostream &stream, PFE<el_t> x)
    {
      stream << x.element;
      return stream; // must return stream
    };


 template<class el_t>
    PFE<el_t> pow(PFE<el_t> a, int exp){
		if( a.iszero() ) return a;
		if( exp == 0 ) return PFE<el_t>(1);
		el_t ea = a.exp_el.right.at(a.element);	
		int nexp = exp*ea;
		while( nexp < 0) {
			
			nexp += a.prime-1;
		}
		nexp = nexp % (a.prime-1) ; 
		return PFE<el_t>( a.exp_el.left.at(nexp ) );
	};	


// prime field element (GF2)
class GF2:
	boost::field_operators< GF2 ,
	boost::equality_comparable< GF2 > > 
{
private:
	int prime = 2;
public: 
	int element;
	
	//!Constructs a prime f
	GF2(){element = 0;};
	GF2(int el): element(el % 2) {};

	//! for initializying the map exponent to element; takes as parameter a primitive element of GF(prime)
		
	bool isempty() const {return (element==prime);}
	void setempty() {element = prime;};

	GF2& operator += (const GF2& x){
		element = (element+x.element) % prime;
		return *this;
	}
	
	GF2& operator -= (const GF2& x){
		element = (element-x.element) % prime;
		if(element < 0) element += prime;
		return *this;
	}

	GF2& operator *= (const GF2& x){
		element = element*x.element % prime;
		return *this;
	}
	              
	GF2& operator /= (const GF2& x){
		return *this;
	}

	operator int (){return (int) element;};

	GF2 inverse() const {
		return *this;
		//return pow(*this, -1);
	}

	bool operator ==(const GF2& x) const {
		return (x.element == element);
	}
	bool operator < (const GF2& x) const {
		return (element < x.element);	
	}
	
	bool iszero() const { return (element == 0);  };

};

ostream &operator<<(ostream &stream, GF2 x)
    {
      stream << x.element;
      return stream; // must return stream
    };

GF2 pow(GF2 a, int exp){
	return a; // 0^0 = 0, 0^1 = 0
};	


#endif
