


>
Summary
=======

This file gives instructions on encoding and decoding information on short DNA molecules.

Description of the scheme
=========================

We store the information on M sequences of length L. The sequences are indexed and protected with error correcting codes. The scheme performs the following steps:

- Multiply the information with a pseudorandom sequence. 

- Split the information into B*k information pieces of length L - log(M)-R, where R is the redundancy on a sequence level. The information of each piece will be stored on a sequence.

- For each block b = 1,...,B, add extra sequences by encoding each row separately by using a Reed-Solomon code over the extension field 2^m. This results in n-k extra sequences per blocks and thus in B*(n-k) extra sequences.

- Add a unique index to each sequence, and inner code each sequence.

We are not avoiding homopolymers. Due to the randomization of the data, long homopolymers and corresponding errors are very unlikely, and the error-correcting code will deal with those. 

Parameters of the code
======================

Here are the parameters of the code, and below are some examples to see how the code can be used to store information on DNA:

- l: length of the index in muliples of 6 bits (default choice is l=4)
- n: length of a block of the outer code (default is n=16383, n must be less or equal than 16383)
- k: number of information symbols of the outer code (default is k=10977, must obey k<=n)
- N: length of inner codeword (default is N=34)
- K: length of inner codeword symbols (default is K=32)
- nuss: number of symbols of outer code per segment (default is nuss=12)

Here are two constraints:
- The index needs to be sufficiently long so that each sequence has a unique index, specifically: l * 6 > log( numblocks * n ), where the logarithm has base 2.
- For the inner and outer code parameters to go together, must have: K * mi = nuss * mo + l * mi, where mi=6, mo=14.

The number of bits per nucleotide is:
2*k/n * (K-l)/N

The redundancy of the outer code is n-k which means that k-n of the sequences of each block can be lost and we can still recover the information. Stated differently if no more than a fraction of 1-k/n of the sequences are lost, then the information can be recovered.

The redundancy of the inner code is N-K which means that a sequence can have (N-K)/2 substitution error and we can still recover it.

Here are a few concrete examples of the choices of the parameters that satisfy the constraints above.
Let the length of the index be l = 4, then we can have at most 2^(mi * l) = 2^24 = 16777216 sequences. That means we can for example choose n = 16383 and store the data on up to 1024 = floor(2^24/16383) many blocks, each consisting of n=16383 sequences. 
Suppose we choose the redundancy of the inner code such that N-K = 3.
Then the following are a subset of the choices that are possible:

| K  | nuss|  N | length of sequence in nucleotides  |
|:--:|:---:|:--:| :-----:|
| 11  | 3   | 14 | 42 |
| 18 | 6   | 21 | 63 |
| 25 | 9   | 28 | 84 |
| 32 | 12  | 35 | 105 |
| 39 | 15  | 42 | 126|
| 46 | 18  | 49 | 147|




Example
=======

In the following, we describe a few examples of how information can be encoded to DNA segments, and decode it back in order to recover the information. For all examples, the first step is to compile the program texttodna. Towards this goal, change to the folder ./simulate and execute the following command which compiles the code:

	`make texttodna' 

The program texttodna can be used to encode data, map it to DNA segments, and to recover the data from the DNA segments as illustrated in the following examples.

Example 1
---------

In our first example, we encode information on one block of length n=12472 and we choose k=9000, which means the outer code can correct nse substitution errors and ner erasures provided that 2*nse + ner <= 12472-9000. The rest of the parameters are the default parameters, i.e.,  N=34, K=32, nuss=12. This results in sequences of length N*ni/2 = 43*6/2 = 102.
 

1. The following command takes the text in the file *data.zip*, encodes it on DNA segments, and stores the segments in *data_encoded.txt*:

	`./texttodna --n=12472 --k=9000 --encode --input=../data/data.zip --output=../data/data_encoded.txt`

2. The following command takes random lines (DNA segments) from *data_encoded.txt*, introduces errors, and saves the resulting file to *data_drawnseg.txt*:

	`./texttodna --disturb --input=../data/data_encoded.txt --output=../data/data_drawnseg.txt`

3. The following command decodes the perturbed data in *data_drawnseg.txt*, and writes the result to *data_rec.zip*; this should recover the original text:

	`./texttodna --decode --n=12472 --k=9000 --numblocks=1 --input=../data/data_drawnseg.txt --output=../data/data_rec.zip`

Example 2
---------

In our second example, we encode information on two blocks of maximum length n=16383, and choose k=12000. In addition we make our sequences shorter than in the previous example, by setting nuss=9. We set K = 25. Note that his satisfies the constraint K*mi = nuss*mo + l*mi, where mi=6, mo=14. Finally, we add more redundancy in the inner code by setting N = 28. This results in sequences of length N*ni/2 = 28*6/2 = 84.

1. Encoding:
	
	`./texttodna --k=12000 --nuss=9 --K=25 --N=28 --encode --input=../data/DNA_channel_statistics.pdf --output=../data/paper_encoded.txt`

2. Decoding:

	`./texttodna --k=12000 --nuss=9 --K=25 --N=28  --decode --numblocks=3 --input=../data/paper_encoded.txt --output=../data/paper_recovered.pdf`


Installation
============

The code is written in C++, and compilation requires installation of the boost library. 


Installation of required software on Linux
------------------------------------------
	sudo apt-get install gcc
	sudo apt-get install make
	sudo apt-get install libboost-all-dev

Installation of required software on macOS
------------------------------------------
	/usr/bin/ruby –e “$curl –fsSL [https://raw.githubusercontent.com/Homebrew/install/master/install](https://raw.githubusercontent.com/Homebrew/install/master/install))”
	brew install boost
	sudo make –f Makefile
	
Licence
==========

All files are provided under the terms of the Apache License, Version 2.0, see the included file "apache_licence_20" for details.
