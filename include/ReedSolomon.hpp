/*
Reed Solomon coding 
*/


#ifndef ReedSolomon
#define ReedSolomon

#include <boost/operators.hpp>
#include <iostream>
#include <ostream>
#include <vector>
#include <cassert>
#include "polynomial.hpp"


using namespace std;



template<class GFE,class Dft>
class RScode{
	public:
	RScode(){};
	RScode(unsigned n_, unsigned k_,const GFE& primitive_el,Dft dft_){
		n = n_; k = k_; a = primitive_el; dftd = dft_;
		n_u = n_; 
		k_u = k_;
	}
	RScode(unsigned n_, unsigned k_,const GFE& primitive_el,Dft dft_,unsigned n_u_){
		n = n_; k = k_; a = primitive_el; dftd = dft_;
		n_u = n_u_; 
		k_u = k + n_u - n;
	}
	unsigned n;
	unsigned n_u;
	unsigned k_u;
	unsigned k;
	GFE a; // primitive element
	Dft dftd; // class providing Fourier transform
	typedef GFE Symbol; // Symbols of the code
	
	
	
	// GF: coefficients of the RS code, GF(p^m)
	// infvec: vector with k entries, containing the information
	vector<GFE>& RSencode(const vector<GFE>& infvec, vector<GFE>& c);
	
	// endode a shortened RS code 
	// infvec: contains the information, length k
	// n: the length of cw
	// n_u: length of the underlying RS code 
	vector<GFE>& RS_shortened_encode(const vector<GFE>& infvec, vector<GFE>& c);

	// GF: coefficients of the RS code, GF(p^m)
	// infvec: vector with k entries, containing the information
	vector<GFE>& RS_systematic_encode(const vector<GFE>& infvec, vector<GFE>& c);

	// decoding when the first n-k entries of the Fourier transform of c are equal zero
	// c: received vector 
	// crec: recovered vector
	// retured is a pair, where pair.fist is the number of erasurs and pair.second the number of errors
	pair<unsigned,unsigned> RSdecode(vector<GFE>& crec, const vector<GFE>& c); 

	// decode a shortened RS code 
	// infvec: contains the information, length k
	// n: the length of cw
	// n_u: length of the underlying RS code 
	pair<unsigned,unsigned> RS_shortened_decode(vector<GFE>& infvec, vector<GFE>& c);

	// decode when information is encoded in the spectrum of c, that is C
	pair<unsigned,unsigned> RS_decode_spec(vector<GFE>& infvec, const vector<GFE>& c);
};


template<class GFE, class Dft > 
vector<GFE>& RScode<GFE,Dft>::RSencode(const vector<GFE>& infvec, vector<GFE>& c){
	assert(infvec.size() == k);

	// encode
	c.resize(n); // c is the codeword
	vector<GFE> C(n); // Fourier transform of codeword
	// first n-k entries are zero
	for(unsigned i=0;i<n-k;++i) C[i] = GFE(0); 
	// next k entries contain information 
	for(unsigned i=0;i<k;++i) {
		C[i+n-k] = infvec[i]; 
	}
	dftd.idft(c,C);
	return c;
}



template<class GFE,class Dft>
//vector<GFE>& RS_shortened_encode(const vector<GFE>& infvec, vector<GFE>& c, const GFE& a, unsigned k, unsigned n, unsigned n_u){
vector<GFE>& RScode<GFE,Dft>::RS_shortened_encode(const vector<GFE>& infvec, vector<GFE>& c){
	//assert();
	vector<GFE> extrazeros(n_u-n, GFE(0));
	vector<GFE> infvec_u = infvec;
	// append zeros to the information vector
	infvec_u.insert(infvec_u.end(), extrazeros.begin(), extrazeros.end());


	// systematic encoding
	RS_systematic_encode(infvec_u,c);//a, k+n_u-n,n_u);
	// the last nu-n entries are zero, so erase them
	c.resize(n);
	return c;
}

// GF: coefficients of the RS code, GF(p^m)
// infvec: vector with k entries, containing the information
template<class GFE,class Dft>
vector<GFE>& RScode<GFE,Dft>::RS_systematic_encode(const vector<GFE>& infvec, vector<GFE>& c){
	assert(infvec.size() == k_u);

	// compute the generator polynomial 
	polynomial<GFE> gp = polynomial<GFE>(1);
	GFE ai = GFE(1); 
	for(unsigned i=1;i <= n_u-k_u; ++i){
		vector<GFE> prov(2);
		prov[0] = GFE(0) - ai;
		prov[1] = GFE(1);
		gp *= polynomial<GFE>(prov);
		ai *= a;
	}
	// generate message polynomial
	polynomial<GFE> mp(infvec);

	// construct the polynomial x^(n-k)
	vector<GFE> xnmkv(n_u-k_u+1,GFE(0));
	xnmkv[n_u-k_u] = GFE(1);
	polynomial<GFE> xnmk(xnmkv);

	polynomial<GFE> sp = mp*xnmk; // the codeword polynomial
	
	// obtain s_r(x) = p(x)*x^(n-k) mod g(x) as res.second 
	pair<polynomial<GFE>, polynomial<GFE> > res = divide<GFE>( sp ,gp);

	// c(x) = p(x)*x^(n-k) - s_r(x)
	sp -= res.second; 
	
	// make sure the cw has length n

	
	c = sp.poly;
	unsigned oldsize = c.size();
	c.resize(n_u);
	for(unsigned i=oldsize;i<n_u;++i) c[i]=GFE(0);
	return c; 
}



// decoding when the first n-k entries of the Fourier transform of c are equal zero
// c: received vector 
// crec: recovered vector
// retured is a pair, where pair.fist is the number of erasurs and pair.second the number of errors
template<class GFE, class Dft>
std::pair<unsigned,unsigned> RScode<GFE,Dft>::RSdecode(vector<GFE>& crec, const vector<GFE>& c){
	
	pair<unsigned, unsigned> erctr;
	
	crec = c; // this will be the recoverd vector

	GFE am = a.inverse();
	assert(n_u == crec.size());
	unsigned d = n_u-k_u+1;
	
	// the positions of the erasures
	vector<unsigned> er_pos(n_u);
	unsigned j=0;
	for(unsigned i=0;i<n_u;++i) 
		if(crec[i].isempty()) {
			er_pos[j] = i;
			j++;
			crec[i] = GFE(0); // set to zero
		}
	er_pos.resize(j);
	
	erctr.first = j;
	// compute the erasure locator polynomial
	polynomial<GFE> elp = polynomial<GFE>(1);
	for(unsigned i=0;i<er_pos.size(); ++i){
		vector<GFE> prov(2);
		/////// here I lose time when computing the powers.. 
		prov[0] = GFE(0) - pow(am,er_pos[i]);
		prov[1] = GFE(1);
		//prov[0] = GFE(1);
		//prov[1] = GFE(0) - pow(a,er_pos[i]);
		elp *= polynomial<GFE>(prov);
	}
	

	// obtain Crec from the received vector crec
	vector<GFE> Crec(n_u);
	dftd.dft(crec,Crec); 

	// compute syndrome from received C
	vector<GFE> syndv(n_u-k_u);
	for(unsigned i=0; i< n_u-k_u  ;++i)
		syndv[i] = Crec[i];
	polynomial<GFE> synd(syndv);

	synd *= elp; 
	synd.poly.resize(n_u-k_u); // mod x^{n-k}

	//// solving the key equation by Euclid's method 

	// construct polynomial x^(d-1)
	vector<GFE> xdm1(d,GFE(0)); 
	xdm1[d-1] = GFE(1);

	polynomial<GFE> rem1 = xdm1;
	polynomial<GFE> rem2 = synd;
	polynomial<GFE> aux1 = polynomial<GFE>(0);
	polynomial<GFE> aux2 = polynomial<GFE>(1);
	//polynomial<GFE> aux2 = elp;//polynomial<GFE>(1);
	
	
	while(! (rem2.degree() < (d-1+j)/2  )  ){
		pair<polynomial<GFE>, polynomial<GFE> > res = divide<GFE>(rem1,rem2);
		polynomial<GFE> aux_new = polynomial<GFE>(0) - res.first*aux2 + aux1;
		// prepare for the next step
		rem1 = rem2;
		rem2 = res.second;
		aux1 = aux2;
		aux2 = aux_new;
	}

	// aux2 is the error locator
	// rem2 is the errata evaluator polynomial
	
	erctr.second = aux2.degree();
	
	//cout << "erasures: " << erctr.first << " errors: " << erctr.second << endl;


	// obtain the errata locator as the product of the error evaluator and the errata evaluator
	elp *= aux2;
	

	//// Forney's algorithm to find the error values

	polynomial<GFE> elpder = elp;
	elpder.derive(); // derivative of the error locator polynomal

	// compute the error values
	
	//vector<GFE> err_val(n);
	GFE ami = GFE(1);	
	GFE ai = GFE(1);
	for(unsigned i=0;i<n_u; ++i){
		if( elp.evaluate(ami) == GFE(0)  ){
			// e_j = a^i \Gamma(a^-i) / \Lambda'(a^-i)
			//err_val[i] = GFE(0) -  ai*( rem2.evaluate(ami)/ elpder.evaluate(ami) );
			//cout << "errata at: " << i << endl;
			// correct the error
			//cout << "corrected at " << i << endl;
			GFE elpderval = elpder.evaluate(ami); 
			if(! elpderval.iszero()){ 
				crec[i] += ai*( rem2.evaluate(ami)/ elpderval);
			}else{ // otherwise something did go wrong, we cannot divide by zero..
				erctr.second = n_u; // this marks an decoding error
			}
		}
		ami *= am;
		ai *= a;
	}
	return erctr;
}


//
template<class GFE, class Dft>
pair<unsigned,unsigned> RScode<GFE,Dft>::RS_shortened_decode(vector<GFE>& infvec, vector<GFE>& c){
	
	vector<GFE> extrazeros(n_u-n, GFE(0));
	// append the zeros to the received vector
	c.insert(c.end(), extrazeros.begin(), extrazeros.end());
	
	vector<GFE> crec(n_u);
	pair<unsigned,unsigned> erctr = RSdecode(crec,c);
	
	// copy the information; recall that c is endoded systematically.. 	
	infvec = vector<GFE>(crec.begin()+n-k, crec.begin()+n );
	infvec.resize(k);
	return erctr;
}


//
template<class GFE, class Dft>
pair<unsigned,unsigned> RScode<GFE,Dft>::RS_decode_spec(vector<GFE>& infvec, const vector<GFE>& c){
	
	pair<unsigned,unsigned> erctr;

	assert(n == c.size());
	vector<GFE> crec; // the recovered codeword
	cout << "start RS decode.. " << endl;
	erctr = RSdecode(crec,c);
	cout << "..end RS decode.. " << endl;

	vector<GFE> Crec(n);
	dftd.dft(crec,Crec); // obtain C from the recovered codeword

	// recover the information from C
	infvec.resize(k);
	for(unsigned i=0;i<k;++i) infvec[i] = Crec[i+n-k]; 
	
	return erctr;
}

#endif
