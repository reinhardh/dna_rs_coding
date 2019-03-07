Summary
=======

This file gives instructions on encoding and decoding information on short DNA molecules.

Description of the scheme
=========================

We store the information on M fragments of length L. We index the fragments, and protect the index with $R$ parity bits. We propose the following scheme:

- Map the information to K fragments of length L - \log(M)-R. 

- Multiply the information with a pseudorandom sequence. 

- Split the information into K blocks of length L - \log(M)-R, and add extra blocks by encoding each row separately by using a Reed-Solomon code over the extension field 2^m. This results in M-K extra molecules.

- Add a unique index in the middle of each fragment, along with parity bits protecting the index.

We are not avoiding homopolymers. Due to the randomization of the data, long homopolymers and corresponding errors are very unlikely, and the error-correcting code will deal with those. 

Parameters of the code
======================

Here are the parameters of the code, below are some examples to see how the code can be used to store information on DNA:

- l: length of index in muliples of 6 bits (default choice l=4)
- n: length of a block of the outer code (default is n=16383, must be less or equal than 16383)
- k: number of information symbols of outer code (default is k=10977, must obey k<=n)
- N: length of inner codeword (default is N=34)
- K: length of inner codeword symbols (default is K=32)
- nuss: number of symbols of outer code per segment (default is nuss=12)

Here are two constraints:
- The index needs to be sufficiently long so that each sequence has a unique index, specificaly: l*6 < log_2( numblocks*n )
- For the inner and outer code parameters to go together, must have: K*mi = nuss*mo + l*mi, where mi=6, mo=14

Example
=======

In the following, we describe a few examples of how information can be encoded to DNA segments, and decode it back in order to recover the information. For all examples, the first step is to compile the program texttodna. Towards this goal, change to the folder ./simulate and execute the following command which compiles the code:

	`make texttodna' 

The program texttodna can be used to encode data, map it to DNA segments, and to recover the data from the DNA segments as illustrated in the following examples.

Example 1
---------

In our first example, we encode information on one block of length n=12472 and we choose k=9000, which meand the outer code can correct nse substitution errors and ner erasures provided that 2*nse + ner <= 12472-9000. The rest of the parameters are the default parameters, i.e.,  N=34, K=32, nuss=12. This results in sequences of length N*ni/2 = 43*6/2 = 102.
 

1. The following command takes the text in the file *Linda_reformat.zip*, encodes it on DNA segments, and stores the segments in *linda_encoded.txt*:

	`./texttodna --n=12472 --k=9000 --encode --input=../data/Linda_reformat.zip --output=../data/linda_encoded.txt`

2. The following command takes random lines (DNA segments) from *linda_encoded.txt*, introduces errors, and saves the resulting file to *linda_drawnseg.txt*:

	`./texttodna --disturb --input=../data/linda_encoded.txt --output=../data/linda_drawnseg.txt`

3. The following command decodes the perturbed data in *linda_drawnseg.txt*, and writes the result to *linda_rec.zip*; this should recover the original text:

	`./texttodna --decode --n=12472 --k=9000 --numblocks=1 --input=../data/linda_drawnseg.txt --output=../data/linda_rec.zip`

Example 2
---------

In our second example, we encode information on two blocks of maximum length n=16383, and choose k=12000. In addition we make our sequences shorter than in the previous example, by setting nuss=9. We set K = 25. Note that his satisfies the constraint K*mi = nuss*mo + l*mi, where mi=6, mo=14. Finally, we add more redundancy in the inner code by setting N = 28. This results in sequences of length N*ni/2 = 28*6/2 = 84.

1. Encoding:
	
	`./texttodna --k=12000 --nuss=9 --K=25 --N=28 --encode --input=../data/DNA_channel_statistics.pdf --output=../data/paper_encoded.txt`

2. Decoding:

	`./texttodna --k=12000 --nuss=9 --K=25 --N=28  --decode --numblocks=3 --input=../data/paper_encoded.txt --output=../data/paper_recovered.pdf`


Installation
============

The code is written in C++, and compillation requires installation of the boost libary. 


Installation of required software on Linux (not tested)
-------------------------------------------------------
	sudo apt-get install gcc
	sudo apt-get install make
	sudo apt-get install libboost-all-dev

Licence
==========

All files are provided under the terms of the Apache License, Version 2.0, see the included file "apache_licence_20" for details.
