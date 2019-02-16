

/*
Example: GF(2^14) 
with primitive polynomial 

x^14 + x^7 + x^5 + x^3 + 1
*/

#include <iostream>
#include "../include/PFE.hpp"
#include "../include/EFE.hpp"
#include "../include/GF2M.hpp"
#include "../include/helpers.hpp"


using namespace std;

//template<pfe> pfe EFE<pfe>::m=2;

// define static variables 

typedef int el_t;
template<> int PFE<el_t>::prime = 2;
//template<el_t> int PFE<el_t>::prime = 11;
template<> boost::bimap<el_t,el_t> PFE<el_t>::exp_el = boost::bimap<el_t,el_t>();

//typedef PFE<el_t> pfe;
typedef GF2 pfe;

template<> polynomial<pfe> EFE<pfe>::prim_poly = polynomial<pfe>();





int main(){


const unsigned m = 14;

const unsigned prim_poly = 16553;
typedef GF2M<unsigned,m,prim_poly> gf;

typedef GF2M<unsigned,6,91> gfi;

cout << "GF2M with m = " << gf::m << endl;
cout << "GF2M with m = " << gfi::m << endl;

gf g1 = gf(9537,0);
cout << "g1: " << g1 << endl;

gf g2 = gf(66,0);
cout << "g2: " << g2 << endl;

cout << "g1*g2: " << g1*g2 << endl;
cout << "g1 / g2 = " << g1 / g2 << endl;
cout << "g2 / g1 = " << g2 / g1 << endl;
cout << "g2 / 0 = " << g2 / gf(0) << endl;
cout << "0 / g2 = " << gf(0) / g2 << endl;




cout << "g1: " << g1 << endl;
cout << "g1^0 : " << pow(g1,0) << endl;
cout << "g1^1: " << pow(g1,1) << endl;
cout << "g1^2: " << pow(g1,2) << endl;
cout << "g1^3: " << pow(g1,3) << endl;
cout << "g1^-1: " << pow(g1,-1) << endl;
cout << "g1^-2: " << pow(g1,-2) << endl;
cout << "g1^-3: " << pow(g1,-3) << endl;





cout << "primitive polynomial: " << bitset<32>(prim_poly) << endl;
cout << "                      " <<  bitset<32>( (1 << m) ) << endl; 

cout << "g1^-1: " << g1.inverse()   << endl;
cout << "g1ˆ-1*g1 " << g1.inverse()*g1 << endl;

cout << "g2^-1: " << g2.inverse()   << endl;
cout << "g2ˆ-1*g2 " << g2.inverse()*g2 << endl;



cout << "test divide: " << endl;
pair<unsigned,unsigned> res = divide(64+32,32+4+1);
cout << "original:  " <<  bitset<32>(64+32) << "\t" << bitset<32>(32+4+1) << endl; 
cout << "result:    " <<  bitset<32>(res.first) << "\t" << bitset<32>(res.second) << endl; 
if( ( gf(res.first)*gf(32+4+1,0) + gf(res.second) ) == gf(64+32,0) )
	cout << "correct!" << endl;


gf g3;
if(g3.isempty())
	cout << "is empty: " << bitset<64>(g3.el) << endl;

gf g4(1);
if(!g4.isempty())
	cout << "is not empty" << endl;




// x^24 + x^16 + x^15 + x^14 + x^13 + x^10 + x^9 + x^7 + x^5 + x^3 + 1
/*
typedef GF2M<24,16901801> gf24;
for(unsigned i = 0; i < 100000; i++){
	gf24 el = gf24(i,0);
	unsigned ord = el.order();
	//if(ord < 16777216)
	cout << i << " " << el.order() << endl;
}
*/

// test the helper function 


unsigned fromv[4] = {4,3,0,1};
vector<unsigned> from(fromv, fromv+4);

vector<unsigned> to;


cout << "from: " << from << endl;
gfm2gfn<unsigned,4,2>(from,to);
cout << "to:   " << to << endl;
gfm2gfn<unsigned,2,4>(to,from);
cout << "from: " << from << endl;



typedef GF2M<uint,8,0> GFSTR;
const unsigned prim_poly_o = 16553;
const unsigned mo = 16;
typedef GF2M<uint,mo,prim_poly_o> GFO;


string str ("This is text we want to store on DNA");

// convert the string to symbols in the outer code
vector<GFSTR> strnu; // string in numbers
for(uint c: str)
	strnu.push_back(GFSTR(c,0));
vector<GFO> datao;


GFM2GFN<GFSTR,GFO>(strnu,datao);

//convert do DNA sequence
// test DNAMapGF

DNAmapGF<GFO> dnamap;
string dnamol;
dnamap.cw2frag(datao,dnamol);

cout << "DNA molecule" << endl;
cout << dnamol << endl;

vector<GFO> datao2;
cout << "and back.. " << endl;
dnamap.frag2cw(datao2,dnamol);
cout << "done back.. " << endl;

//vector<GFO> datao2 = datao;

// convert the string back
vector<GFSTR> strnu2; // string in numbers
GFM2GFN<GFO,GFSTR>( datao2 ,strnu2);
string str2;
for(GFSTR ch : strnu2)
	str2 += (char)ch.el;

if(str == str2){
	cout << "Here is the original character sequence:" << endl;
	cout << str2 << endl;
} else {
	cout << ">" << str << "<" << endl;
	cout << ">" << str2 << "<" << endl;
}


	typedef GF2M<uint,14,0> GFOO;
	typedef GF2M<uint,42,0> GFINT; // 56/14 = 4



	vector<GFOO> filesize(3,GFOO(0,0)); // length is 3
	vector<GFINT> filesizeint(1,GFINT( 4179 ,0));
	GFM2GFN<GFINT,GFOO>( filesizeint ,filesize);
	filesize.resize(3); // make sure length is 3
	cout << "filesizeint: " << filesizeint << endl;
	cout << "filesize " << filesize << endl;


	vector<GFINT> filesizeint2(0,0);
	vector<GFOO> filesizegfo = filesize;
	cout << filesizegfo << endl;
	cout << GFOO::m << endl;
	cout << GFINT::m << endl;
	
	cout << "filesizegfo " << filesizegfo << endl;
	GFM2GFN<GFOO,GFINT>(filesizegfo,filesizeint2);
	
	cout << "filesizeint " << filesizeint2 << endl;

	unsigned filesizeend = filesizeint2[0].el;
	cout << "Determined filesize as: " << filesizeend << endl;





}
