/*
discrete fourier transform for extension fields
*/


#ifndef DFT
#define DFT 

#include <boost/operators.hpp>
#include <iostream>
#include <ostream>
#include <vector>
#include <cassert>
#include "polynomial.hpp"

using namespace std;

// class to compute the powers of an element a in an extention field, 
// via a lookup table 
template<class co_t>
class fieldpower{

public:
unsigned n; // max power  
co_t nm;

private:
co_t a; // a is an element of order n
co_t am; // the inverse of a; also has order n
vector<co_t> powers_a; // the lookup table for the powers of a
vector<co_t> powers_am; // the lookup table for the powers of a^(-1)

public: 
	fieldpower(){};
	fieldpower(const unsigned n_, const co_t& a_){
		a = a_;
		n = n_;
		
		nm = co_t(n);
		nm = pow( (co_t) nm,-1  );

		// initialize powers_a
		powers_a.resize( n );
		powers_a[0] = co_t(1);
		powers_a[1] =	a;
		for(int i =2;i < powers_a.size();++i)
			powers_a[i] = powers_a[i-1]*a;

		// initialize powers_am
		am = pow( (co_t) a,-1);
		powers_am.resize( n );
		powers_am[0] = co_t(1);
		powers_am[1] =	am;
		for(int i =2;i < powers_am.size();++i)
			powers_am[i] = powers_am[i-1]*am;
		
		cout << "fieldpower: lookup tables generated.. " << endl;
	}

	co_t copow(int powe) const {
		if(powe >= 0){
			return powers_a[ powe % n ];
		} else
			return powers_am[(-powe) % n ];
	}

};


template<class co_t>
void dft_primitive( const vector<co_t>& xv, vector<co_t>& Xv, fieldpower<co_t>& fp, int prefac){
	//assert(xv.size() == Xv.size());
	for(int j=0; j < Xv.size();++j){
		// compute Xv[j]
		Xv[j] = co_t(0);
		for(int i=0;i < xv.size(); ++i){
			//Xv[j] += xv[i]*pow( (coeff_t)  a, (int)  i*j);
			if(prefac<0)
				Xv[j] += fp.nm*xv[i]*fp.copow(-i*j);
			else
				Xv[j] += xv[i]*fp.copow(i*j);
		}
	}
};

// prefac = 1: normal fft, prefac =-1: inverse fft
template<class co_t>
void fft( const vector<co_t>& xv, vector<co_t>& Xv, const fieldpower<co_t>& fp, int prefac, unsigned Q, unsigned P){
	assert(Q*P == xv.size());
	Xv = vector<co_t>(xv.size());

	// compute y_p^(r)
	vector<vector< co_t> > Y(Q,vector< co_t >(P, co_t(0) ));
	for(unsigned r=0;r<Q;++r){
		// for each r, a Fourier transform..
		for(unsigned p=0;p<P;++p){
			for(unsigned q=0;q<Q;++q){
				Y[r][p] += xv[P*q+p]*fp.copow(prefac*P*q*r);
			}
			Y[r][p] *= fp.copow(prefac*p*r);
		}
	}
	// compute A_{Qs+r}
	for(unsigned r=0;r<Q;r++){
	for(unsigned s=0;s<P;s++){
		Xv[Q*s+r] = co_t(0);
		for(unsigned p=0;p<P;++p){
			Xv[Q*s+r] += Y[r][p]*fp.copow(prefac*Q*s*p); 
		}
		if(prefac == -1) Xv[Q*s+r] *= fp.nm;
	}
	}
}


////////////////////////////////////////////////
// FFT implementaion, uses fieldpower 
template<class co_t>
class DFT_FFT{
private:
fieldpower<co_t> fp;
unsigned n; // the length of the dft
unsigned P;
unsigned Q;

public: 
	DFT_FFT(const unsigned n_, const co_t& a_,unsigned P_, unsigned Q_){
		P = P_;
		Q = Q_;
		fp = fieldpower<co_t>(n_,a_);
	}
	DFT_FFT(){}; 

	void dft( const vector<co_t>& xv, vector<co_t>& Xv) const {
		fft<co_t>(xv,Xv,fp,1,P,Q);
	};

	void idft(vector<co_t>& xv, const vector<co_t>& Xv) const {
		fft<co_t>(Xv,xv,fp,-1,P,Q);
	}

};

////////////////////////////////////////////////



// DFT implementation of certain lenth that uses a precomputed lookup table to compute the powers of a
template<class co_t>
class DFT_LA{
private:
co_t a; // the Fourier kernal
co_t am;
unsigned n; // the length of the dft
co_t nm;
vector<co_t> powers_a; // the lookup table for the powers of a
vector<co_t> powers_am; // the lookup table for the powers of a^(-1)

public: 
	DFT_LA(const unsigned n_, const co_t& a_){
		a = a_;
		n = n_;
		
		nm = co_t(n);
		nm = pow( (co_t) nm,-1  );

		// initialize powers_a
		powers_a.resize( ((n-1)*(n-1)+1) );
		powers_a[0] = co_t(1);
		powers_a[1] =	a;
		for(int i =2;i < powers_a.size();++i)
			powers_a[i] = powers_a[i-1]*a;

		// initialize powers_am
		am = pow( (co_t) a,-1);
		powers_am.resize( ((n-1)*(n-1)+1) );
		powers_am[0] = co_t(1);
		powers_am[1] =	am;
		for(int i =2;i < powers_am.size();++i)
			powers_am[i] = powers_am[i-1]*am;
		
		cout << "lookup tables generated.. " << endl;
	}

	void dft( const vector<co_t>& xv, vector<co_t>& Xv) const {
		Xv = vector<co_t>( xv.size() );
		for(int j=0; j < Xv.size();++j){
			// compute Xv[j]
			Xv[j] = co_t(0);
			for(int i=0;i < xv.size(); ++i){
				//Xv[j] += xv[i]*pow( (coeff_t)  a, (int)  i*j);
				Xv[j] += xv[i]*powers_a[i*j];
			}
		}
	};

	void idft(vector<co_t>& xv, const vector<co_t>& Xv) const {
		xv = vector<co_t>( Xv.size() );
		for(int j=0; j < xv.size();++j){
			// compute Xv[j]
			//cout << "adf"<< endl;
			xv[j] = co_t(0);
			//cout << j << endl;
			for(int i=0;i < Xv.size(); ++i){
				//cout << "--" << i << endl;
				//cout << i << "-- Xv: " <<  Xv[i] << " powers " << powers_am[i*j] << endl;
				//cout << "-- xv: " << xv[j] << endl; 
				xv[j] += nm*Xv[i]*powers_am[i*j];
			}
		}
	}

};

// primitive DFT 
template<class co_t>
class DFT_PRIM{
private:
co_t a; // the Fourier kernal
co_t am;
unsigned n; // the length of the dft
co_t nm;

public: 
	DFT_PRIM(){};
	DFT_PRIM(const unsigned n_, const co_t& a_){
		a = a_;
		n = n_;
		
		nm = co_t(n);
		nm = pow( (co_t) nm,-1  );
	}

	void dft( vector<co_t>& xv, vector<co_t>& Xv) const {

		//assert(xv.size() == Xv.size());
		for(int j=0; j < Xv.size();++j){
			// compute Xv[j]
			Xv[j] = co_t(0);
			for(int i=0;i < xv.size(); ++i){
				//Xv[j] += xv[i]*pow( (coeff_t)  a, (int)  i*j);
				Xv[j] += xv[i]*pow( (co_t) a, i*j);
			}
		}
	};

	void idft(vector<co_t>& xv, vector<co_t>& Xv) const {
		for(int j=0; j < Xv.size();++j){
			// compute Xv[j]
			xv[j] = co_t(0);
			for(int i=0;i < xv.size(); ++i)
				xv[j] += nm*Xv[i]*pow( (co_t) a,i*j);
		}
	}

};


// to do: implement the faster FFT
template<class coeff_t>
void dft( vector<coeff_t>& xv, vector<coeff_t>& Xv, coeff_t a ){

	
	cout << "startcomplookup" << endl; 
	unsigned nn = Xv.size();
	//compute lookup table
	vector<coeff_t> powers( ((nn-1)*(nn-1)+1) );
	powers[0] = coeff_t(1);
	powers[1] =	a;
	for(int i =2;i < powers.size();++i)
		powers[i] = powers[i-1]*a;
	cout << "finishcomplookup" << endl; 

	
	assert(xv.size() == Xv.size());
	for(int j=0; j < Xv.size();++j){
		// compute Xv[j]
		Xv[j] = coeff_t(0);
		for(int i=0;i < xv.size(); ++i){
			//Xv[j] += xv[i]*pow( (coeff_t)  a, (int)  i*j);
			Xv[j] += xv[i]*powers[i*j];
		}
	}
};


// nm: n^-1 ,where n is the lenght of the dft
// am: a^-1, where a is an element of order n
/*
template <class coeff_t>
void idft(vector<coeff_t>& xv, vector<coeff_t>& Xv,coeff_t am, coeff_t nm ){
	assert(xv.size() == Xv.size());
	for(int j=0; j < Xv.size();++j){
		// compute Xv[j]
		xv[j] = coeff_t(0);
		for(int i=0;i < xv.size(); ++i)
			xv[j] += nm*Xv[i]*pow( (coeff_t) am, (int)  i*j);
	}
	
}
*/


template <class coeff_t>
void idft(vector<coeff_t>& xv, vector<coeff_t>& Xv,coeff_t a ){
	assert(xv.size() == Xv.size());

	unsigned nn = xv.size();
	coeff_t n = coeff_t(nn);
	cout << "n: " << n << endl;
	n = pow(n,-1);
	
	
	//compute lookup table
	cout << "startcomplookup" << endl; 
	coeff_t am = pow( (coeff_t) a,-1);
	vector<coeff_t> powers( ((nn-1)*(nn-1)+1) );
	powers[0] = coeff_t(1);
	powers[1] =	am;
	for(int i =2;i < powers.size();++i)
		powers[i] = powers[i-1]*am;

	cout << "finishlookuptable.. " << endl;

	for(int j=0; j < Xv.size();++j){
		// compute Xv[j]
		xv[j] = coeff_t(0);
		for(int i=0;i < xv.size(); ++i){
			
			//xv[j] += n*Xv[i]*pow( (coeff_t) am, (int)  i*j);
			xv[j] += n*Xv[i]*powers[i*j];
		}
	}
	
}

#endif
