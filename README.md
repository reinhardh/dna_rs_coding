Summary
=======

This file gives instructions on encoding and decoding information on short DNA molecules.

Description of the scheme
=========================

The coding scheme relies on the following:
- The substitution error probability is small, but substitution errors are much more likely than deletions and insertions (by about a factor 10). So we only account for them with an outer code.
- If the information that is written on the nucleotides behaves like a random sequence, then long homopolymers are very unlikely; specifically the probability that there is a run of length at least $k$ in a random DNA molecule of length $n$ can be upper bounded by $n (4/4^k) = n 4^{1-k}$. 


We store the information on $M$ fragments of length $L$. We index the fragments, and protect the index with $R$ parity bits. We propose the following scheme:

- Map the information to $K$ fragments of length $L - \log(M)-R$. 

- Multiply the information with a pseudorandom sequence. 

- Split the information into $K$ blocks of length $L - \log(M)-R$, and add extra blocks by encoding each row separately by using a Reed-Solomon code over the extension field $2^m$. This results in $M-K$ extra molecules.

- Add a unique index in the middle of each fragment, along with parity bits protecting the index.

We are not avoiding homopolymers. Due to the randomization of the data, long homopolymers and corresponding errors are very unlikely, and the RS code will deal with those. 

Parameters of the code
======================


This particular scheme:
- Information is in char = 2^8 bit
- Map to (4^3)^14
- Map to sequence of 4^3 = 2^6, so 3 char to 4 symbols
- Map to k = 13926 fragments of length (4^3)^14 = 2^28, so map 7 byte to two fragments
- Outer decode
- Muliply each block with (same) pseudorandom sequence
- Add index of length 4 * 4^3 
- Inner decode (add 2 Symbols redundancy)

- l: length of index in muliples of 6 bits (default choice l=4)
- n: length of a block of the outer code (default is n=16383, must be less or equal than 16383)
- k: number of information symbols of outer code (default is k=10977, must obey k<=n)
- N: length of inner codeword (default is K=20)
- K: length of inner codeword symbols (default is K=18)
- nuss: number of symbols of outer code per segment (default is nuss=6)

Here are two constraints:
- The index needs to be sufficiently long so that each sequence has a unique number: l*6 < log_2( numblocks*n )
- For the inner and outer code parameters to go together, must have: K*mi = nuss*mo + l*mi, where mi=6, mo=14

Example
=======

In the following, we describe an example workflow that takes part of Archimedes Book, encodes it to DNA segments, and decodes it in order to recover the information. 

To start, change to the folder ./simulate and execute the following commands:

- compile texttodna: (make sure that the correct path to BOOST_LIB is provided in the Makefile)
	
	`make texttodna' 

The program texttodna can be used to encode data, map it to DNA segments, and to recover the data from the DNA segments as follows.

1. The following command takes the text in the file *archimedes_short.txt*, encodes it on DNA segments, and stores the segments in *arch_short_dna.txt*:

	`./texttodna --n=12472 --k=9000 --encode --input=../data/Linda_reformat.zip --output=../data/linda_encoded.txt`

2. The following command takes random lines (DNA segments) from *arch_short_dna.txt*, introduces errors, and saves the resulting file to *arch_short_drawnseg.txt*:

	`./texttodna --disturb --input=../data/linda_encoded.txt --output=../data/linda_drawnseg.txt`

3. The following command decodes the perturbed data in *arch_short_dna.txt*, and writes the result to *arch_short_rec.txt*; this should recover the original text:

	`./texttodna --decode --n=12472 --k=9000 --numblocks=1 --input=../data/linda_drawnseg.txt --output=../data/linda_rec.zip`



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
