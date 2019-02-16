


#ifndef HELPERSDNA 
#define HELPERSDNA

#include <boost/operators.hpp>
#include <iostream>
#include <ostream>
#include <boost/bimap.hpp>
#include <string>
#include <map>
#include<vector>

using namespace std;


template<class uint>
unsigned degree(uint v){ 
	uint r = 0;
	while (v >>= 1) {
		r++;
	}
	return r;
}


template<class cf_t> 
ostream &operator<<(ostream &stream, const vector<cf_t>& x) 
{ 
	for(unsigned i = 0; i<x.size(); ++i) stream << x[i] << "  ";
	return stream; // must return stream 
}; 



template<class vectt>
unsigned hamming_distance(vectt& v1,vectt& v2){
	assert(v1.size()=v2.size());
	unsigned hd = 0;
	for(unsigned i=0;i<v1.size();++i)
		if(v1[i] != v2[i]) hd++;
	return hd;
}

template<class vectt>
void flipvecdir(vectt& v){
	vectt tmp = v;
	for(unsigned i=0;i<v.size();++i){
		v[i] = tmp[v.size()-1-i];
	}
}


// convert decimal number to a number with another base
// nu: 		an unsigned integer
// base: 	base of the new representation
// size: 	size of the vector required for the new representation
template<class base_t>
vector<base_t> tobase(uint64_t nu, uint64_t base, unsigned size){
	
	vector<base_t> nn(size,(base_t) 0 ); 
	if(nu == 0) return nn;
	
	uint64_t div = nu;
	unsigned i = 0;
	while(div != 0){
		//nn[i] = base_t(div % base); // rest
		assert( i < nn.size() );
		nn[i] = base_t(div % base);
		div = div/base;	
		i++;
	}
	
	return nn;
};


// convert vector of numbers in some basis to decimal number
template<class base_t>
uint64_t frombase(const vector<base_t>& nn,uint64_t base){
	unsigned deg = degree(base);
	assert(deg*nn.size() < 64); // otherwise we have an overflow below
	uint64_t nu = 0;
	uint64_t basepow = 1;
	for(unsigned i=0;i<nn.size();++i){
		uint64_t adf = nn[i];
		nu += ((uint64_t)adf) * basepow;//((unsigned) (adf)) * ((int) pow((double) base,(int) i));
		basepow *= base;
	}
	return nu;
};


template<class GFM,class GFN>
void GFM2GFN(vector<GFM>& from, vector<GFN>& to){
	// choose a and b as smallest numbers such that a*m == b*n
	unsigned a = 1;
	unsigned b = 1;
	while (a*GFM::m != b*GFN::m){
		if(a*GFM::m < b*GFN::m)
			a += 1;
		else
			b += 1;	
	}
	while( (from.size() % a) != 0 ) from.push_back(GFM(0));
    to.resize(0);
    // every a entries of `from' are mapped to b entries of `to'
	for(unsigned i=0;i<from.size()/a; ++i){
		// the unsigned char type below is important, otherwise errors occur when converting..
		vector<uint> ind;
		for(typename vector<GFM>::iterator it = from.begin()+i*a; it != from.begin()+i*a+a; ++it){
			ind.push_back(it->el);
		}		
		unsigned mm = GFM::m;
		assert(mm < 64);
		uint64_t powmm = 1;
		powmm = powmm << mm;
		uint64_t pownm = 1;
		unsigned nm = GFN::m;
		assert(nm < 64);
		pownm = pownm << nm;
		int64_t	indec = frombase( ind , powmm );     
        vector<uint> tmpvec = tobase<uint>(indec, pownm ,b);
		for(uint j : tmpvec)
			to.push_back(GFN(j,0));
		//to.insert(to.begin()+i*b,tmpvec.begin(),tmpvec.end());
	}
};


template<class uint,unsigned m, unsigned n>
void gfm2gfn(vector<uint>& from, vector<uint>& to){
	// choose a and b as smallest numbers such that a*m == b*n
	unsigned a = 1;
	unsigned b = 1;
	while (a*m != b*n){
		if(a*m < b*n)
			a += 1;
		else
			b += 1;	
	}
	while((from.size() % a) != 0 ) from.push_back(0);
   
    to.resize(0);
    // every a entries of `from' are maped to three entries of `to'
	for(unsigned i=0;i<from.size()/a; ++i){
        // the unsigned char type below is important, otherwise errors occur when converting..
		unsigned indec = frombase(vector<unsigned char>(from.begin()+i*a, from.begin()+i*a+a), 1<<m );       
        vector<uint> tmpvec = tobase<uint>(indec,1<<n,b);
        to.insert(to.begin()+i*b,tmpvec.begin(),tmpvec.end());
	}
};




// converts each 2 characters in data_char to 3 elements of GF(47), which are appended at 
template<class pfe>
void char2pfe(string& data_char, vector<pfe>& data_b47){
    if((data_char.size() % 2) != 0 ) data_char.push_back('\n');
   
    data_b47.resize(0);
	// convert vector of char to vector with base 47;
    // every 2 char's are coded to three pfe's 
    for(unsigned i=0;i<data_char.size()/2; ++i){
        // the unsigned char type below is important, otherwise errors occur when converting..
		unsigned indec = frombase(vector<unsigned char>(data_char.begin()+i*2, data_char.begin()+i*2+2),256);
        
        vector<pfe> tmpvec = tobase<pfe>(indec,47,3);
        data_b47.insert(data_b47.begin()+i*3,tmpvec.begin(),tmpvec.end());
	}
};


template<class pfe>
void pfe2char(string& data_char, vector<pfe>& data_b47){
    
    assert(data_b47.size() % 3  == 0);

    data_char.resize(2*(data_char.size()/3) );
    // convert vector of char to vector with base 47;
    // every 2 char's are coded to three pfe's 
    for(unsigned i=0;i<data_b47.size()/3; ++i){
        unsigned indec = frombase(vector<pfe>(data_b47.begin()+i*3, data_b47.begin()+i*3+3),47);
        
        vector<char> tmpvec = tobase<char>(indec,256,2);
        data_char.insert(data_char.begin()+i*2,tmpvec.begin(),tmpvec.end());
    }

};




/*
maps each element of GF2M to a string with letters {A,C,G,T}
*/

template<class GF>
class DNAmapGF{
public: 
	DNAmapGF(){assert(GF::m % 2 == 0);};	
	// codeword to DNA fragment
	void cw2frag(const vector<GF>& cw, string& frag ){
		frag.resize(0);
		for(GF gf: cw){
			for(unsigned i=0; i<GF::m; i+=2){
				switch( (gf.el >> i) & 3 ) { // shift by two, and then check the bits..
    				case 0 : frag.append("A");
							 break;
					case 1 : frag.append("C"); 
							 break;
					case 2 : frag.append("G");
							 break;
					case 3 : frag.append("T");
							 break;
				}
			}
		}
	}
	
	// DNA fragment to codeword
	void frag2cw(vector<GF>& cw, const string& frag ){
		assert( (frag.size() % (GF::m/2)) == 0 );
		cw.resize( frag.size() / (GF::m/2) );
		for(unsigned i=0;i<cw.size();++i){ 
			typename GF::uint el = 0;
			for(unsigned j=0; j<GF::m/2; j++){
				switch(frag[i*GF::m/2 + j]){	
    				case 'A' :  // add zero..  
							 break;
					case 'C' : el += ((typename GF::uint) 1) << 2*j ;
							 break;
					case 'G' : el += ((typename GF::uint) 2) << 2*j ;
							 break;
					case 'T' : el += ((typename GF::uint) 3) << 2*j ;
							 break;
				}
			}
			cw[i] = GF(el,0);
		}
	}
};






/*
maps each element of GF(47) to a string with letters {A,C,G,T}
*/

template<class PFE>
class DNAmap{
private: 
	//typedef boost::bimap< PFE , std::string > mapdna_type;
	//mapdna_type mapdna;
	
	map<PFE,string> pfetostr;
	map<string,PFE> strtopfe;

public: 
	DNAmap(){
		//initialize map
		unsigned prime = 47;
		
		char nucl[] = "ACGT"; // nucleotides	

		vector<string> allpos(4*4*4);
		for(unsigned i=0;i<allpos.size();++i){
			char cur[] = "AAA";
			cur[0] = nucl[i % 4];
			cur[1] = nucl[(i/4) % 4];
			cur[2] = nucl[(i/16) % 4];
			allpos[i] = string(cur);
		}
		
		unsigned j = 0;
		for(unsigned i=0;i<prime;++i){
			while( allpos[j][1] == allpos[j][2]) j++;
			//mapdna.insert( typename mapdna_type::value_type(PFE(i), allpos[j] ));
			pfetostr[PFE(i)] = allpos[j];
			strtopfe[allpos[j]] = PFE(i);
			
			j++;
		}
		//print_map(mapdna.right, " ", cout);
	}
	
	// codeword to DNA fragment
	void cw2frag(const vector<PFE>& cw, string& frag ){
		frag.resize(0);
		for(unsigned i=0;i<cw.size();++i){
			frag.append(pfetostr[cw[i]]); //mapdna.left.at(cw[i]); 
		}
	}
	
	// DNA fragment to codeword
	void frag2cw(vector<PFE>& cw, const string& frag ){
		assert( (frag.size() % 3) == 0   );
		cw.resize(frag.size()/3);
		for(unsigned i=0;i<cw.size();++i){ 
		
			//char charartmp[] = "AAA"; // nucleotides	
			//charartmp[0] = frag[i*3];
			//charartmp[1] = frag[i*3+1];
			//charartmp[2] = frag[i*3+2];
			//string tmp(charartmp);
			string tmp = string(frag.begin()+i*3, frag.begin()+i*3+3  );
				
			//print_map(mapdna.right, " ", cout);
			
			//cout <<  mapdna.right.at(tmp )  << endl; 
			//cw[i] = PFE( mapdna.right.at(tmp ) ); 
			cw[i] = strtopfe[tmp];
			//cout << cw[i] << endl;
		}
	}
};

#endif
