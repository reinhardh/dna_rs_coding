
#include <vector>
#include "./helpers.hpp"


#include <algorithm>
#include <numeric>
#include<boost/foreach.hpp>
#include <random>


#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>

using namespace boost::accumulators;
using namespace std;





/*
Innercode: over GF2M
*/

template<class Innercode, class Outercode>
class EnDecode {
	//private:
	public:
	void outercode_stats(const vector< vector< typename Innercode::Symbol > >& segm_ordered);
	Innercode innercode;
	Outercode outercode;
	typedef typename Outercode::Symbol GFO;
	typedef typename Innercode::Symbol GFI;
	typedef typename Outercode::Symbol::uint uint; 
	DNAmapGF<GFI> dnamap;
	vector< vector<GFI> > groundtruth;
	const unsigned l; // length of index
	//const uint gfim = GFI::m;
	//const unsigned len = 4;
	//const uint nbitsindex = 24;
	//const uint nbitsindex = 24;
	typedef GF2M<uint,24,0> GFUINT; // l*GFI::m = number bits of index
	typedef GF2M<uint,42,0> GFINT; // 56/14 = 4

	unsigned m;
	unsigned nuss; // number of symbols of outer code per molecule
    unsigned numblocks; 

	EnDecode(){};
	EnDecode(const Innercode& innercode_, const Outercode& outercode_ ,unsigned l_,unsigned nuss_):
	l(l_)
	{
		innercode = innercode_;
		outercode = outercode_;
		nuss = nuss_;
		assert(innercode.n>l);
		m = innercode.k - l;
	
	};

	void encode(string& str, vector<string>& urn);

	void decode(string& str, const vector<string>& drawnseg, const string& gtruthfile);
};




////////////////////////////
template<class Innercode, class Outercode>
void EnDecode<Innercode,Outercode>::encode(string& str, vector<string>& urn){
	
	typedef GF2M<uint,8,0> GFSTR;

	// convert the string to symbols in the outer code

	vector<GFSTR> strnu; // string in numbers
	for(unsigned char c: str){
		unsigned nu = (unsigned) c;
		assert(nu < 256);
		strnu.push_back(GFSTR( nu , 0));
	}
	vector<GFO> datao;

	GFM2GFN<GFSTR,GFO>(strnu,datao); // convert string

	// add the filelength at the beginning:
	assert(GFO::m >= 10); // then the filesize can be 2^(3*10)
	vector<GFO> filesize(3,GFO(0,0)); // length is 3
	vector<GFINT> filesizeint(1,GFINT( str.size() ,0));
	GFM2GFN<GFINT,GFO>( filesizeint ,filesize);
	filesize.resize(3); // make sure length is 3
	cout << "Filesize in Byte: " << str.size() << endl;
	datao.insert(datao.begin() , filesize.begin() , filesize.end() );

	// determine blocklength
	numblocks = 1;
	while( datao.size() > nuss * outercode.k * numblocks ) numblocks++;
	// datao has to grow such that nuss * Outercode::k * numblocks	

	while( datao.size() > nuss * outercode.k * numblocks ) datao.push_back(GFO(0));

	urn.resize(numblocks*outercode.n);

	// add pseudorandom numbers 
	mt19937 rng( 5489U ) ; // constructed with seed 5489U
	uniform_int_distribution<uint> randgfe(0, 1<<GFO::m - 1 );
	
    for(unsigned b=0;b<numblocks;++b){ // run over blocks 
		cout << "encode block " << b << endl;	
		// encode outer codes, add pseudorandom numbers	
		vector< vector<GFO> > c(nuss); // outer codeword for given block
    	for(unsigned nb=0;nb<nuss;++nb){ 
			// outer information vector
        	vector<GFO> infvec( datao.begin()+(b*nuss+nb)*outercode.k ,
		                   		datao.begin()+(b*nuss+nb+1)*outercode.k );
			// encode the information, obtain outer codeword
			outercode.RS_shortened_encode(infvec,c[nb]); 
			// add pseudorandom element to each element of the codeword
			
			for(GFO& gfo: c[nb])
				gfo += GFO(randgfe(rng),0);
				// gfo += GFO( rng() % (1<<GFO::m) ,0);
 		}

		// add index, encode inner code
        for(unsigned i=0;i<outercode.n;++i){ // run over elements in a block 
			
			vector<GFO> infmol(nuss);
            for(unsigned nb=0;nb<nuss;++nb){
				infmol[nb] = c[nb][i];
			}
			
			// tranform from GFO to GFI
			vector<GFI> infmolin;
			GFM2GFN<GFO,GFI>(infmol,infmolin);
			//cout << "infmolin: " << infmol << endl;
            // compute index
			vector<GFUINT> index(1, GFUINT( i+outercode.n*b ,0) );
			vector<GFI> indexgfi;
			GFM2GFN<GFUINT,GFI>(index,indexgfi);
		
			// add pseudorandom sequence to index 
			mt19937 rng( 5489U ) ; // constructed with seed 5489U
			uniform_int_distribution<uint> randgfe(0, 1<<GFI::m - 1 );
			for(GFI& gfi: indexgfi)
				gfi += GFI(randgfe(rng),0);
				// gfi += GFI( rng() % (1<<GFI::m) ,0);

			// construct information of inner codeword: [information | index | information]
			//vector<GFI> infveci; 	
			//infveci.insert(infveci.end(), infmolin.begin(), infmolin.begin() + infmolin.size()/2);
			//infveci.insert(infveci.end(), indexgfi.begin(), indexgfi.end());
			//infveci.insert(infveci.end(), infmolin.begin() + infmolin.size()/2, infmolin.end());
			
			// construct information of inner codeword: [information | index ]
			vector<GFI> infveci;

			// index first
			infveci.insert(infveci.end(), indexgfi.begin(), indexgfi.end());
			infveci.insert(infveci.end(), infmolin.begin(), infmolin.end());

			assert(infveci.size() == innercode.k);
            // encode the inner code, write it to urn[i] 
            vector<GFI> icw;  
			innercode.RS_shortened_encode(infveci,icw);
			dnamap.cw2frag(icw,urn[i+outercode.n*b]);
			flipvecdir(urn[i+outercode.n*b]);
        }
		cout << "encoded block " << b << " of " << numblocks << endl;
    } 	
}


template<class T>
bool pairCompare(const T & x, const T & y) {
  return x.second < y.second; 
}


////////////////////////////
template<class Innercode, class Outercode>
void EnDecode<Innercode,Outercode>::decode(string& str, const vector<string>& drawnseg, const string& gtruthfile){
	/*
	Inner decoding, for each read:
		- inner decode, obtain inner information, extract index
		- for each index we keep a map which counts how many times an information has been
		  observed with that index. In that map, we only keep the information that has a
		  minimal number of errors, out ot all molecules we have seen with that index
	After going though all reads, choose the one for a given index, that has been observed
	most frequently.
	*/


   	////
	// read grouundtruth if available for making tests
	if(gtruthfile.size() > 0){
		vector<string> drawnseg_gtruth;
		string infile = "../../dna_coding_messy/File1_ODNA.txt";
		string sLine = "";
		ifstream in;
		in.open(gtruthfile.c_str());
		while (!in.eof()){
			getline(in, sLine);
			drawnseg_gtruth.push_back(sLine);
		}
		drawnseg_gtruth.resize(drawnseg_gtruth.size()-1); // erase the last, empty line
		// run over all reads
		for(const string& read: drawnseg_gtruth){
			vector<GFI> rececw;		// read: inner codeword
			string flipedread = read;
			flipvecdir(flipedread);
			dnamap.frag2cw(rececw,flipedread); 	// map read to GFI vector
		
			// decode molecule
			pair<unsigned,unsigned> erctrictmp; // for error counts
			vector<GFI> irec(innercode.k); // recovered information from inner codeword
			erctrictmp = innercode.RS_shortened_decode(irec,rececw);

		 	groundtruth.push_back(irec);
		}
		cout << "Read groundtruth for evaluation, " << drawnseg_gtruth.size() << " sequences" << endl;
	}

	//ofstream cand_file("candidates.dat");
	////

	// to keep track    
	pair<unsigned,unsigned> erctric(0,0);
	vector<int> frag_freqs(numblocks*outercode.n, 0); // to count the frequencies of each DNA sequence
	unsigned ic_decerr = 0; // counts the inner decoding errors
	unsigned errctr = 0;
	unsigned numfrag = 0;

	typedef typename Innercode::Symbol GFI;
	typedef map<vector<GFI>, unsigned > seqctrmap; // information -> how many times it appeared
	vector< seqctrmap > segm_ordered_map(numblocks*outercode.n, seqctrmap() );
	// the minimal number of errors as reported by inner decoder for all molecules with given index
	vector<unsigned> curerr(numblocks*outercode.n, innercode.n ); // error of the element in segm_ordered_map;

	// run over all reads
	for(const string& read: drawnseg){
		vector<GFI> rececw;		// read: inner codeword
		string flipedread = read;
		flipvecdir(flipedread);
		dnamap.frag2cw(rececw,flipedread); 	// map read to GFI vector
		
		// decode molecule
		pair<unsigned,unsigned> erctrictmp; // for error counts
		vector<GFI> irec(innercode.k); // recovered information from inner codeword
		erctrictmp = innercode.RS_shortened_decode(irec,rececw);
		
		
		if(erctrictmp.second == innercode.n_u){ // decoding error
			ic_decerr++;
		} else {
			erctric.first += erctrictmp.first;
			erctric.second += erctrictmp.second;
			// extract index
			//// information part of the inner codeworkd is [information | index | information]
			//unsigned len = (innercode.k-l)/2;
			//vector<GFI> indexgfi(irec.begin()+len, irec.begin()+len+l);
			// information part of the inner codeworkd is [information | index]
			//vector<GFI> indexgfi(irec.end()-l-1, irec.end());
			vector<GFI> indexgfi(irec.begin(), irec.begin() + l);

			// remove pseudorandom sequence from index 
			mt19937 rng( 5489U ) ; // constructed with seed 5489U
			uniform_int_distribution<uint> randgfe(0, 1<<GFI::m - 1 );
			for(GFI& gfi: indexgfi)
				gfi -= GFI( randgfe(rng) ,0);
				//gfi -= GFI( rng() % (1<<GFI::m)  ,0);

			vector<GFUINT> index;
			GFM2GFN<GFI,GFUINT>(indexgfi,index);
			//unsigned maxindex = 1<<(GFI::m*l);
			//unsigned ind = ((int)index[0].el - maxindex/3); // % maxindex;
			//cout << ind << " " << maxindex << endl;
			unsigned ind = (unsigned)index[0].el; // index
			if(ind >= numblocks*outercode.n){
				// if the index is not in the appropriate range, there must be an error
				errctr ++;
			} else { // index can be right
				if(erctrictmp.second < curerr[ind]){ 
					// found information with smaller error than previously seen, thus discard map
					curerr[ind] = erctrictmp.second;
					segm_ordered_map[ind].empty();
					segm_ordered_map[ind].insert(make_pair(irec,1));
				} else if (erctrictmp.second == curerr[ind]){
					typename seqctrmap::iterator cur = segm_ordered_map[ind].find(irec);
					if(cur != segm_ordered_map[ind].end()){ // segment already appeared  
						cur->second++;  
					} else {
						segm_ordered_map[ind].insert(make_pair(irec,1));
					}	
				} // else: (erctrictmp.second > curerr[ind]), nothing to do
			
				// if groundtruth given, compare to groundtruth whether there is a decoding error
				/*
				if(! groundtruth.empty() ){ 
					if(groundtruth[ind] != irec) {
						errctr++;
					} else {
						// this fragment has been decoded without error
						frag_freqs[ind]++;
					}
				}
				*/
				/////////////////////////
				/*
				// write to candidate output, as a test:
				assert(irec.size() == innercode.k);
				// encode the inner code, write it to urn[i] 
				vector<GFI> icw; 
				innercode.RS_shortened_encode(irec,icw);
				string str;
				dnamap.cw2frag(icw,str);
				flipvecdir(str);
				cand_file << str << endl;
				*/
				/////////////////////////

			}	
		
		} 
		numfrag++;
	} // run over all reads

	//cand_file.close();


	// take the one (out of all segments with smallest error) with largest frequency as the estimate,
	// store it in segm_ordered
	vector< vector<GFI> > segm_ordered(numblocks*outercode.n, vector<GFI>());
	for(unsigned ind=0;ind<segm_ordered_map.size();++ind){
		if(! segm_ordered_map[ind].empty() ){
			typename seqctrmap::iterator max = max_element(segm_ordered_map[ind].begin(), segm_ordered_map[ind].end(), pairCompare<typename seqctrmap::value_type> );
			segm_ordered[ind] = max->first;
		}
	}


	
	//cout << "read error prob: " << (float) errctr / (float) numfrag << endl;
	cout << "number of reads: " << numfrag << endl;
	//cout << "inner code: " << (float)erctric.first / (float) numfrag << " erasures on average corrected" << endl;
	cout << "inner code: " << (float)erctric.second/ (float) numfrag << " errors on average corrected per sequence" << endl;
	//cout << "inner code, dec. errors: " << ic_decerr << endl;

	int ct = 0;
	for(unsigned i=0;i<curerr.size();++i){
		//cout << i << " " << curerr[i] << endl;
		if( curerr[i] > (float)(innercode.n - innercode.k)/2.0 ){
		//if( curerr[i] >= (innercode.n - innercode.k) ){
			segm_ordered_map[i] = seqctrmap();
			ct++;
			segm_ordered[i] = vector<GFI>();
		}
	}
	// cout << "errasures: " << ct << endl;
	
	if(gtruthfile.size() > 0){
		// error statistics based on groundtruth
		outercode_stats(segm_ordered);
	}

	///// decode the blocks
	vector<GFO> inf;

	// for removing pseudorandom numbers 
	mt19937 rng( 5489U ) ; // constructed with seed 5489U
	uniform_int_distribution<uint> randgfe(0, 1<<GFO::m - 1 );

	for(unsigned b=0;b<numblocks;++b){ // run over blocks 		
		// obtain list of codewords
		vector< vector<GFO> > cw(nuss, vector<GFO>(outercode.n) ); // outer codewords for given block
		unsigned erasure_ctr = 0;
		for(unsigned i =0;i<outercode.n;++i){ // run over elements in a block 
			if( ! segm_ordered[b*outercode.n+i].empty()){

				// extract information vector
				const vector<GFI>& inf = segm_ordered[b*outercode.n + i];
				//unsigned len = (innercode.k-l)/2;
				vector<GFI> infi;
				// information is [information | index]
				//infi.insert(infi.end(), inf.begin(), inf.begin() + len);
				//infi.insert(infi.end(), inf.begin()+len+l+1, inf.end());
				//infi.insert(infi.end(), inf.begin(), inf.end()-l);
				infi.insert(infi.end(), inf.begin()+l, inf.end());
				vector<GFO> info;
				GFM2GFN<GFI,GFO>(infi,info); // information vector in 

				for(unsigned nb=0;nb<nuss;++nb)
					cw[nb][i] = info[nb];
	
			} else { // else we have an erasure
				for(unsigned nb=0;nb<nuss;++nb)
					cw[nb][i] = GFO();
				erasure_ctr++;
			}
		}

		cout << erasure_ctr << " many erasures in block " << b << " of length " << outercode.n <<  endl;
	
		///// decode outer codword
		cout << "remove random sequence" << endl;
		for(unsigned nb=0;nb<nuss;++nb){
			// add pseudorandom element to each element
			for(GFO& gfo: cw[nb])
				if( ! gfo.isempty() )
					gfo -= GFO(randgfe(rng),0);
					//gfo -= GFO( rng() % (1<<GFO::m) ,0);
				else // errasure
					//rng(); 
					randgfe(rng); // draw to increase state of rng
		}

		cout << "decoding block " << b << endl;
		vector< vector<GFO> > infrec(nuss);
		for(unsigned nb=0;nb<nuss;++nb){
			//pair<unsigned,unsigned> erctroc = outercode.RS_decode_spec(infrec[nb],cw[nb]);
			pair<unsigned,unsigned> erctroc = outercode.RS_shortened_decode(infrec[nb],cw[nb]);
			cout << "\tpart " << nb << "/" << nuss << endl;
			if(nb == nuss-1){
				cout << "\touter code: " << erctroc.first << " errasures, " << erctroc.second << " errors corrected"<< endl;
			}
		}
		
		// add information to string
		for( const vector<GFO>& cn : infrec )
			inf.insert( inf.end(), cn.begin(), cn.end() );

	} // for(unsigned b=0


	vector<GFINT> filesizeint(0,0);
	vector<GFO> filesizegfo(inf.begin(),inf.begin()+3);
	GFM2GFN<GFO,GFINT>(filesizegfo,filesizeint);
	unsigned filesize = filesizeint[0].el;
	inf = vector<GFO>(inf.begin()+3, inf.end());

	typedef GF2M<uint,8,0> GFSTR;
	// convert the string to symbols in the outer code
	vector<GFSTR> strnu; // string in numbers
	GFM2GFN<GFO,GFSTR>(inf,strnu);
	
	str.resize(0);
	for(GFSTR ch : strnu)
		str += (char)ch.el;
	str.resize(filesize);

}







///// determine error probabilities outer code, based on groundtruth 
template<class Innercode, class Outercode>
void EnDecode<Innercode,Outercode>::outercode_stats(const vector< vector<typename Innercode::Symbol> >& segm_ordered){
	unsigned total_erasure_ctr = 0;
	unsigned total_error_ctr = 0;
	unsigned total_blockerror_ctr = 0;
	//// determine the number of errors and erasures in each block
	for(unsigned b=0;b<numblocks;++b){
		unsigned block_erasure_ctr = 0;
		unsigned block_error_ctr = 0;
		for(unsigned i =0;i<outercode.n;++i){ // run over elements in a block 
			if( ! segm_ordered[b*outercode.n+i].empty()){
            	//compare to groundtruth
				if(! groundtruth.empty() ){ 
					if(groundtruth[b*outercode.n+i] != segm_ordered[b*outercode.n+i]) block_error_ctr++;
				}
			} // else we have an erasure
            else{
                block_erasure_ctr++;
            }	
		}
		total_erasure_ctr += block_erasure_ctr;
		total_error_ctr += block_error_ctr;

		cout << b <<"-> errors: " << block_error_ctr << "	erasures: " << block_erasure_ctr <<"	total:" << block_error_ctr*2+block_erasure_ctr << endl;
		
		if(block_error_ctr*2+block_erasure_ctr > outercode.n-outercode.k){
			total_blockerror_ctr++;
		}
	}

	cout << "num block errors: " << total_blockerror_ctr << endl;
	cout << "OC: symbol err. prob.: " << (float) total_error_ctr / (float) (outercode.n*numblocks) << endl;
	cout << "OC: symbol era. prob.: " << (float) total_erasure_ctr / (float) (outercode.n*numblocks) << endl;

	cout << "OC: symbol err. prob.: " << total_error_ctr  << endl;
	cout << "OC: symbol era. prob.: " << total_erasure_ctr << endl;
}
