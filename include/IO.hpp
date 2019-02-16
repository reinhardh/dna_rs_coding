

#include "../include/helpers.hpp"



/*
Read the DNA sequences from a .fastq file which has the format:
something
DNA
something
something
something
DNA
something
...
*/

class fastq_iterator
  : public std::iterator<std::input_iterator_tag, std::string>
{
private:
  std::istream* is;
  std::string line;
  void fetchline() {
    if (is && !std::getline(*is, line)) is = 0;
  }

public:
  fastq_iterator(std::istream& is) : is(&is) {
	fetchline();
    fetchline(); // skip one line
  }
  fastq_iterator() : is(0) { }
  std::string operator*() {
    return line;
  }
  fastq_iterator& operator++() {
    fetchline(); // skip three lines 
    fetchline();
    fetchline();
    fetchline();
    return *this;
  }
  bool operator==(const fastq_iterator& rhs) const {
    return (!rhs.is) == (!is);
  }
  bool operator!=(const fastq_iterator& rhs) const {
    return !(rhs == *this);
  }
};






void fliplett(string& frag){
	for(unsigned i=0;i<frag.size();++i){
		switch(frag[i]){
			case 'T': frag[i] = 'A';
				break;
			case 'A': frag[i] = 'T';
				break;
			case 'G': frag[i] = 'C';
				break;
			case 'C': frag[i] = 'G';
				break;
		}
	}
}


void readflip(string infile, vector<string>& drawnseq, bool rev){
        ifstream in;
        in.open(infile.c_str());
        string sLine = "";
        vector<int> hist(200,0);
        unsigned totalctr = 0;

        int ctr = 0;
        while (!in.eof()){
            getline(in, sLine);
            ctr++;
            if(ctr % 4 == 2){
                hist[sLine.size()]++;
                totalctr++;

                if(rev){
                    flipvecdir(sLine);
                    fliplett(sLine);
                }

                if((sLine.size() == 117)  ){
                    sLine.resize(117);
                    drawnseq.push_back(sLine);
                }

            }
            //if(ctr == 1000000) break;
        }

        cout << "Total # sequences: " << totalctr << endl;
        cout << "# of seq. of len. 117: " << hist[117] << endl;
}

