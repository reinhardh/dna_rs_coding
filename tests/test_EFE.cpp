

/*
Example: GF(11^3) 
with primitive polynomial 10 x^2 + 5 x + 7
*/

#include <iostream>
#include "../include/PFE.hpp"
#include "../include/EFE.hpp"


using namespace std;


//template<pfe> pfe EFE<pfe>::m=2;

// define static variables 

typedef int el_t;
template<> int PFE<el_t>::prime = 11;
//template<el_t> int PFE<el_t>::prime = 11;
template<> boost::bimap<el_t,el_t> PFE<el_t>::exp_el = boost::bimap<el_t,el_t>();

typedef PFE<el_t> pfe;

template<> polynomial<pfe> EFE<pfe>::prim_poly = polynomial<pfe>();





int main(){

	// initialize
	pfe::initialize_exp_el(7);

	typedef EFE<pfe> efe;

	// initialize the primitive polynomial
	vector<pfe> ppv(4);
	ppv[0] = 7;
	ppv[1] = 4;
	ppv[2] = 0;
	ppv[3] = 1;
	EFE<pfe>::prim_poly = polynomial<pfe>(ppv); 

	
	vector<pfe> v1(3);
	v1[0] = 0;
	v1[1] = 5;
	v1[2] = 10;
	efe e1 = efe(v1);

	vector<pfe> v2(3);
	v2[0] = 1;
	v2[1] = 0;
	v2[2] = 5;
	efe e2 = efe(v2);

	vector<pfe> v3(3);
	v3[0] = 0;
	v3[1] = 3;
	v3[2] = 7;
	efe e3 = efe(v3);

	vector<pfe> v4(3);
	v4[0] = 6;
	v4[1] = 5;
	v4[2] = 4;
	efe e4 = efe(v4);


	
	cout << "e1: " << e1 << endl;
	cout << "e1^1: " << pow(e1,1) << endl;
	cout << "e1^2: " << pow(e1,2) << endl;
	cout << "e1^3: " << pow(e1,3) << endl;


	
	cout << "primitive polynomial: " << EFE<pfe>::prim_poly << endl;
	cout << "e1: " << e1 << endl;
	cout << "e2: " << e2 << endl;

	cout << "e1 * e2 = " << e1 * e2 << endl;

	
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

	//pfe a1 = pfe(4);
	//pfe a2 = pfe(10);

	//cout << a2/a1 << endl;


}
