

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

	// initialize
	//pfe::initialize_exp_el(1);

	typedef EFE<pfe> efe;
	const unsigned m = 14;
	// initialize the primitive polynomial
	// pfe ppvv[m+1] = {1,0,0,0,0,0,0,1,0,1,0,1,0,0,1}; // for m=14
	// x^14 + x^7 + x^5 + x^3 + 1
	pfe ppvv[m+1] = {1,0,0,1,0,1,0,1,0,0,0,0,0,0,1}; // for m=14
	vector<pfe> ppv(ppvv, ppvv+m+1);
	cout << "shoud be 1 ? " << ppv[5] << endl;
	EFE<pfe>::prim_poly = polynomial<pfe>(ppv); 

	pfe v1v[m] = {1,0,0,0,0,0,1,0,1,0,1,0,0,1};
	vector<pfe> v1(v1v, v1v+m);
	efe e1 = efe(v1);

	pfe v2v[m] = {0,1,0,0,0,0,1,0,0,0,0,0,0,0};
	vector<pfe> v2(v2v, v2v+m);
	efe e2 = efe(v2);

	pfe v3v[m] = {1,0,0,0,0,0,1,0,1,0,1,0,0,0};
	vector<pfe> v3(v3v, v3v+m);
	efe e3 = efe(v3);

	pfe v4v[m] = {1,1,0,1,0,0,1,0,1,0,1,0,0,0};
	vector<pfe> v4(v4v, v4v+m);
	efe e4 = efe(v4);

	pfe v5v[m] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	vector<pfe> v5(v5v, v5v+m);
	efe e5 = efe(v5);



	
	cout << "e1: " << e1 << endl;
	cout << "e1^0: " << pow(e1,0) << endl;
	cout << "e1^1: " << pow(e1,1) << endl;
	cout << "e1^2: " << pow(e1,2) << endl;
	cout << "e1^3: " << pow(e1,3) << endl;
	cout << "e1^-1: " << pow(e1,-1) << endl;
	cout << "e1^-2: " << pow(e1,-2) << endl;
	cout << "e1^-3: " << pow(e1,-3) << endl;


	
	cout << "primitive polynomial: " << EFE<pfe>::prim_poly << endl;
	cout << "e1: " << e1 << endl;
	cout << "e2: " << e2 << endl;

	cout << "e1 * e2 = " << e1 * e2 << endl;
	cout << "e1 / e2 = " << e1 / e2 << endl;
	cout << "e2 / e1 = " << e2 / e1 << endl;	
	cout << "e2 / 0 = " << e2 / e5 << endl;
	cout << "0 / 1 = " << e5 / e2 << endl;


	
	cout << "------" << endl;
	cout << e1 << " * " << e3 << endl;
	cout <<  " = " << e1*e3 << endl;
	cout << "------" << endl;
	cout << "the same? "<< (e1 == e2) << endl;

	cout << "------" << endl;
	cout << "order(e1): " << e1.order() << endl;
	cout << "order(e2): " << e2.order() << endl;
	cout << "order(e3): " << e3.order() << endl;
	cout << "order(e4): " << e4.order() << endl;


	cout << "e4: " << e4 << endl;
	
	efe e4i = e4.inverse();
	
	cout << "size: " << e4i.el.poly.size() << endl;


	cout << "e4^-1: " << e4i   << endl;
	cout << "e4: " << e4   << endl;
	
	efe prod = e4i*e4;
	
	cout << "prod size: " << prod.el.poly.size() << endl;

	cout << "prod  : " << prod  << endl;

	cout << "=================" << endl;

	cout << "e3: " << e3   << endl;
	cout << "e3^-1: " << e3.inverse()   << endl;
	cout << "e3*e3^-1  : " << e3*e3.inverse()  << endl;

	cout << "=================" << endl;

	cout << "e2: " << e3   << endl;
	cout << "e2*e2^-1  : " << e2*e2.inverse()  << endl;

	cout << "e3/e3: " << e3 / e3 << endl;


///// 

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



}
