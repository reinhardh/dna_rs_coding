/*
Program to translate text to DNA and vice versa
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/program_options.hpp>
#include "../include/GF2M.hpp"
#include "../include/DFT.hpp"
#include "../include/helpers.hpp"
#include "../include/encodedecode.hpp"
#include "../include/ReedSolomon.hpp"
#include <string>
#include <fstream>
#include <streambuf>

using namespace std;

namespace po = boost::program_options;

int main(int ac, char* av[])
{


// command line options

int numblocks = 1;
string infile;
string outfile;
string gtruthfile;


unsigned N = 34;
unsigned K = 32;

unsigned n = 16383;
unsigned k = 10977; // k = 0.67*n

unsigned l = 4; // length of the index
unsigned nuss = 12; // number of symbols of outer code per segment


int opt;
po::options_description desc("Allowed options");
desc.add_options()
    ("help", "produce help message")
    ("encode", "encode")
    ("decode", "decode")
    ("disturb", "draw uniformly at random from the input lines, add errors to each line")
	("input",po::value<string>(&infile)->default_value(""),"inputfile")	
	("output",po::value<string>(&outfile)->default_value(""),"outputfile")	
	("groundtruth",po::value<string>(&gtruthfile)->default_value(""),"groundtruthfile")	
	
	("numblocks",po::value<int>(&numblocks)->default_value(1),"numblocks")	
	("n",po::value<unsigned>(&n)->default_value(16383),"n")	
	("k",po::value<unsigned>(&k)->default_value(10977),"k")	
	("N",po::value<unsigned>(&N)->default_value(34),"N")	
	("K",po::value<unsigned>(&K)->default_value(32),"K")	
	("l",po::value<unsigned>(&l)->default_value(4),"l")	
	("nuss",po::value<unsigned>(&nuss)->default_value(12),"nuss")	
;

po::variables_map vm;
po::store(po::parse_command_line(ac, av, desc), vm);
po::notify(vm);    

if (vm.count("help")) {
    cout << desc << "\n";
    return 1;
}


///// parameters

typedef unsigned uint;
// inner code 
const unsigned mi = 6;
const unsigned prim_poly = 91;
typedef GF2M<uint,mi,prim_poly> GFI; // GF(2^6) with primitive polynomail 91
const unsigned N_u = 63; // the underlying length of the shortened inner code


GFI fai = GFI(2,0); // 64 + 2 = 1000010 = x + x^5
DFT_FFT<GFI> dftgfi(N_u,fai,7,9); // 7*9 = 63
typedef RScode<GFI, DFT_FFT<GFI> > Innercode;
Innercode innercode(N,K,fai,dftgfi,N_u);

// outer code
const unsigned Q = 127;
const unsigned P = 129;
const unsigned n_u = P*Q; // 2^14 - 1


const unsigned prim_poly_o = 16553;
const unsigned mo = 14;
typedef GF2M<uint,mo,prim_poly_o> GFO;

GFO fao = GFO(66,0); // 64 + 2 = 1000010 = x + x^5
DFT_FFT<GFO> dftgfo(n_u,fao,P,Q); // Fourier transform for the outer code 

typedef RScode<GFO,DFT_FFT<GFO> > Outercode;
Outercode outercode(n,k,fao,dftgfo,n_u);



// check all the parameter choices
if (K*mi != nuss*mo+l*mi) {
	cerr << "Must have: K*mi = nuss*mo + l*mi, where mi=6, mo=14, but got " << K*mi << ", " << nuss*mo+l*mi << endl;	
	return 1;
}

if ( K > N) {
	cerr << "Must have: K <= N"  << endl;	
	return 1;
}

if ( k > n) {
	cerr << "Must have: k <= n"  << endl;	
	return 1;
}

if ( n > n_u) {
	cerr << "Must have: n <= n_u"  << endl;	
	return 1;
}

if ( N > N_u) {
	cerr << "Must have: N <= N_u"  << endl;	
	return 1;
}

if ( l*mi < log( numblocks*n )/log( 2 ) ) {
	cerr << "Index has lenght " << l*mi << " bit, but require " << log( numblocks*n )/log( 2 ) << " bits" << endl;	
	return 1;
}



// tell the user the parameter choices

cout << "--------------------------------" << endl;
cout << "redundancy outer code: " << float(k)/float(n) << endl;
cout << "redundancy inner code: " << float(K)/float(N) << endl;
cout << "--------------------------------" << endl;



// encoder/decoder 
EnDecode< Innercode , Outercode > endecode(innercode,outercode,l,nuss);


/////////////////////// encode 
if (vm.count("encode")) {
	
	cout << "start encoding.." << endl;

	if(infile == "" || outfile ==""){
		cout << "in/outfile not specified " << endl; 
		return 0;
	}
	cout << "infile:  " << infile << endl;
	cout << "outfile: " << outfile << endl;

	// read data
	std::ifstream t(infile.c_str());
	std::string str((std::istreambuf_iterator<char>(t)),std::istreambuf_iterator<char>());

	// encode - determines the number of blocks required automatically
	vector<string> urn(n*numblocks);
	endecode.encode(str, urn);
	numblocks = endecode.numblocks; 
	

	if(numblocks*k*mo*nuss < str.size()*8){
		cerr << "trying to store " << str.size()*8 << " bits, but can only store " << numblocks*k*mo*nuss << "many" << endl;	
		return 1;
	}


	cout << "encoded " << str.size() << " Bytes to " << numblocks << " blocks, resulting in "
	<< n*numblocks << " DNA segments of length " << N << " each." << endl;
   
	ofstream out;
	out.open(outfile.c_str());
	for(unsigned i=0;i<urn.size();++i) out << urn[i] << endl;
	out.close();
	return 0;
}


/////////////////////// decode 
if (vm.count("decode")) {
	
	if(infile == "" || outfile =="" || numblocks==0 ) {
		cout << "in/outfile/numblocks not specified " << endl; 
		return 0;
	}
	cout << "infile:  " << infile << endl;
	cout << "outfile: " << outfile << endl;
	cout << "numblocks: "<<numblocks << endl; 

	vector<string> drawnseg;

	
	string sLine = "";
	ifstream in;
	in.open(infile.c_str());
	while (!in.eof()){
		getline(in, sLine);
		drawnseg.push_back(sLine);
	}
	drawnseg.resize(drawnseg.size()-1); // erase the last, empty line

	string recstr;
	cout << "start decode.." << endl;	
	endecode.numblocks = numblocks;
	endecode.decode(recstr, drawnseg, gtruthfile);
	
	ofstream out;
	out.open(outfile.c_str());
	out << recstr;
	out.close();
	cout << "Wrote file of length " << recstr.size() << " Bytes" << endl;
}

////////////////////////// disturb
if (vm.count("disturb")) {

	// take M random draws uniformly at random and disturb those

	const unsigned M = n*numblocks*6;
	const float substprob = 0.0001; // substitution error probability
	cout << "disturb:" << endl;
	cout << "\tDraw " << M << " many times" << endl;
	cout << "\tsubstitution error probability " << substprob << endl;


	if(infile == "" || outfile =="") {
		cout << "in/outfile not specified " << endl; 
		return 0;
	}
	cout << "infile:  " << infile << endl;
	cout << "outfile: " << outfile << endl;
	
	vector<string> urn;
	
	
	string sLine = "";
	ifstream in;
	in.open(infile.c_str());
	while (!in.eof()){
		getline(in, sLine);
		urn.push_back(sLine);
	}
	urn.resize(urn.size()-1); // erase the last, empty line
	
	boost::mt19937 rng; 	
	boost::uniform_int<> unif(0,urn.size()-1); // distribution that maps to 0,..,urn.size()-1
	boost::uniform_int<> unif_N(0,N-1); // uniform distribution over {0,..,N-1}
	boost::uniform_int<> unif_4(0,4-1); // uniform distribution over {0,1,2,3}
	boost::bernoulli_distribution<> bern(substprob);
	//boost::bernoulli_distribution<> faircoin(0.5);
	
	char nucl[] = "ACGT";

	vector<string> drawnseg(M);
	
	ofstream out;
	out.open(outfile.c_str());
	
	// draw M times
	string tmpstr;
	for(unsigned i=0;i<M;++i){
		unsigned randind = unif(rng);
		tmpstr = urn[randind];
		
		// introduce errors
		// introduce on error per inner cw
		for(unsigned j=0;j<1;++j){ 			
			tmpstr[unif_N(rng)] = nucl[unif_4(rng)];
		}
	
		unsigned ctr = 0;
		for(unsigned j=0;j<tmpstr.size();++j)
			if(bern(rng)) {
				tmpstr[j] = nucl[unif_4(rng)];
				ctr++;
			}
			//cout << ctr << endl;
		//if(faircoin(rng)) flipvecdir(tmpstr); // flip every second 

		// write to file
		out << tmpstr;
		if(i!= M-1) out << endl;
	}
	out.close();

}

}
