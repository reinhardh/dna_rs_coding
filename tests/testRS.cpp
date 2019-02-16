/*
Test RS code over binary extension field
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/program_options.hpp>
#include "../include/PFE.hpp"
#include "../include/EFE.hpp"
#include "../include/DFT.hpp"
#include "../include/helpers.hpp"
#include "../include/GF2M.hpp"
//#include "../include/encodedecode.hpp"
#include "../include/ReedSolomon.hpp"
#include <string>
#include <fstream>
#include <streambuf>
#include <random>
#include <chrono>




using namespace std;


// define static variables 

typedef GF2 pfe;
template<> polynomial<pfe> EFE<pfe>::prim_poly = polynomial<pfe>();
template<> unsigned EFE<pfe>::m = 30;








int main(int ac, char* av[])
{

// randomness
random_device r;
default_random_engine gen(r());
uniform_int_distribution<int> randbit(0, 1);
uniform_int_distribution<int> randint(0, 16383);



// parameters
const unsigned M = 50000;
//const unsigned Q = 127;
//const unsigned P = 129;
const unsigned Q = 127;
const unsigned P = 337;
const unsigned n = P*Q; // 2^14 - 1
//const unsigned k = 13926; // k = 0.85*n
const unsigned k = 36379; // k = 0.85*n



///// GF2M definitions
/*
const unsigned prim_poly = 16553;
const unsigned m = 14;
typedef GF2M<unsigned,m,prim_poly> gf;
// efe a = efe(34); // 32 + 2 = 100010 = x + x^5
gf fa = gf(66,0); // 64 + 2 = 1000010 = x + x^5
*/

const uint_least64_t prim_poly = 4399239010919;
const unsigned m = 42;
typedef GF2M<uint_least64_t,m,prim_poly> gf;
// efe a = efe(34); // 32 + 2 = 100010 = x + x^5
gf fa = gf(1935755,0); // 64 + 2 = 1000010 = x + x^5


cout << "primitive polynomial: " << bitset<64>(prim_poly) << endl;
cout << "                      " << bitset<64>( ((uint_least64_t)1 << m) ) << endl; 
cout << "order of a: " << fa.order() << endl;



DFT_FFT<gf> dftgf(n,fa,P,Q); // Fourier transform for the outer code 
// Does the DFT work? 
vector<gf> ct_ = vector<gf>(n);
for(gf& el:ct_) 
	el = gf(randint(gen),0); // random vector
ct_[3] = gf(0);
vector<gf> C_, hatc_;
cout << ct_.size() << " " << ct_[0] << endl;
dftgf.dft(ct_,C_);
dftgf.idft(hatc_,C_);
if(hatc_ == ct_)
	cout << "DFT works!" << endl;
else
	cout << "DFT failed." << endl;


RScode<gf,DFT_FFT<gf>> outercode(n,k,fa,dftgf);


///// EFE version

/*
typedef EFE<pfe> gf;
typedef EFE<pfe> efe;
EFE<pfe>::m = 14;
// initialize the primitive polynomial
// x^14 + x^7 + x^5 + x^3 + 1
pfe ppvv[m+1] = {1,0,0,1,0,1,0,1,0,0,0,0,0,0,1}; // for m=14
vector<pfe> ppv(ppvv, ppvv+m+1);
EFE<pfe>::prim_poly = polynomial<pfe>(ppv);
// initialize the element of order n, the Fourier kernel  
pfe aa[m] = {0,1,0,0,0,0,1,0,0,0,0,0,0,0}; // element of order 2^14-1
vector<pfe> avv(aa, aa+m); // 
efe a = efe(avv); // this is an element of order n
DFT_FFT<gf> dftefe(n,a,P,Q); // Fourier transform for the outer code 



RScode<gf,DFT_FFT<efe>> outercode(n,k,a,dftefe);
*/



///// test the outer code

// generate information

vector<gf> infvec(outercode.k); // information vector; 
for(gf& el: infvec)
	el = gf(randbit(gen));

// encode

vector<gf> c; // outer codeword 
outercode.RSencode(infvec,c);

// disturb

float errorprob = 0.01;
float erasureprob = 0.1;
bernoulli_distribution bern(errorprob);
bernoulli_distribution eras(erasureprob);

int errctr = 0;
int eractr = 0;
for(gf& s : c){
	if( bern(gen) ){	
		// generated random efe
		//vector<pfe> rsymbol(m);
		//for(pfe& el : rsymbol)
		//	el = pfe( randbit(gen) );
		//s = efe(rsymbol);
		s = gf(randint(gen));
		//cout << rsymbol << endl;
		errctr++;
	} else if( eras(gen) ){
		s = gf();
		eractr++;
	}
}

cout << "inserted " << errctr << " errors and " << eractr << " erasures" << endl;
cout << "should be able to correct " << (n - k) / 2 << " errors" << endl;

// decoding

auto t1 = chrono::high_resolution_clock::now();
vector<gf> infvecrec;
pair<unsigned,unsigned> erctroc = outercode.RS_decode_spec(infvecrec,c);
auto t2 = chrono::high_resolution_clock::now();
cout << "decoding took "
     << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
     << " milliseconds\n";

cout << "outer code: " << erctroc.first << " errasures, " << erctroc.second << " errors corrected"<< endl;


if( infvec ==  infvecrec)
	cout << "I did correct all errors!" << endl; 
else	
	cout << "I couldn't correct all errors!" << endl; 


}
