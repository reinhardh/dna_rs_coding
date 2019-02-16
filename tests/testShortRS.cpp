/*
Test inner shortened RS code over binary extension field
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/program_options.hpp>
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




int main(int ac, char* av[])
{


// GF(2^6) with primitive polynomail 91
const unsigned m = 6;
const unsigned prim_poly = 91;
typedef GF2M<m,prim_poly> gf;

// randomness
random_device r;
default_random_engine gen(r());
uniform_int_distribution<int> randbit(0, 1);
uniform_int_distribution<int> randint(0, 63);


// parameters
const unsigned N = 20;
const unsigned N_u = 63; // the underlying length of the shortened inner code
const unsigned K = 18;
const unsigned Qi = 9;
const unsigned Pi = 7; //P*Q = 2^6 - 1


///// GF2M definitions

gf fa = gf(2,0); // 64 + 2 = 1000010 = x + x^5

cout << "primitive polynomial: " << bitset<32>(prim_poly) << endl;
cout << "                      " <<  bitset<32>( (1 << m) ) << endl; 
cout << "order of a: " << fa.order() << endl;

DFT_FFT<gf> dftgf(N_u,fa,Pi,Qi); // Fourier transform for the outer code 
//DFT_PRIM<gf> dftgf(N_u,fa);


//RScode<gf,DFT_PRIM<gf>> innercode(N,K,fa,dftgf,N_u);
RScode<gf, DFT_FFT<gf> > innercode(N,K,fa,dftgf,N_u);

///// test the inner code

// generate information

vector<gf> infvec(innercode.k); // information vector; 
for(gf& el: infvec)
	el = gf(randbit(gen),0);



// encode

vector<gf> c; // outer codeword 
innercode.RS_shortened_encode(infvec,c);

// disturb

uniform_int_distribution<int> randcoef(0, c.size()-1);

// insert one random error
// c[randcoef(gen)] = gf(randbit(gen),0);

// insert one random error
c[randcoef(gen)] = gf();
c[randcoef(gen)] = gf();


// decoding

auto t1 = chrono::high_resolution_clock::now();
vector<gf> infvecrec;
pair<unsigned,unsigned> erctroc = innercode.RS_shortened_decode(infvecrec,c);
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
