Summary
=======

This file gives instructions for encoding and decoding information on short DNA molecules.

Instructions to encode the first Biohackers series from DNA
===========================================================

In August 2020, the first episode of the new Netflix series Biohackers was stored on DNA, see [https://de.biohackersnetflix.com/](https://de.biohackersnetflix.com/) for details. With the code in this repository, you can reconstruct the first episode yourself from DNA. Here are three steps for reconstruction:

1. Download the DNA from [https://de.biohackersnetflix.com/](http://d3u0hl24lyh7px.cloudfront.net/encoded_withprimers.txt), and make sure it is in your `Downloads` folder.

2. Install the code in this repository: There are a few options for installing the code, explained below. Here is a single, simple suggestion that works for most operating systems:
	
	i.  Install Docker, by downloading [Docker Desktop](https://docs.docker.com/get-docker/).
	
	ii. Download a virtual machine containing the pre-installed code from this repository, by typing the following command in the Terminal or bash application on a computer with a  Mac or linux operating system, respectively:
	
	`docker pull rheckel/dna_storage`
	
	iii. Start the virtual machine:
	
	`docker run -it -v ~/Downloads:/data rheckel/dna_storage bash`
	
	You are now runing a virtual machine in which the code from this repository is already installed.

3. In the following you'll decode the series from the DNA. For this you need to type the following commands in the Terminal application at a Mac or the bash application at linux operating system:
	
	i. Change to the folder in which the code is located:
	
	`cd dna_rs_coding/simulate/`

	ii. Type the following to start the decoding:
	
	`./texttodna --decode --n=16383 --k=12700 --N=35 --K=32 --nuss=12 --numblocks=237 --input=/data/netflix_encoded_shuffled.txt --output=../data/netflix_BH101_FullEP_rec.mp4`
	
	You can now see the code running and the video is being reconstructed. Once the code is done running, the file `netflix_BH101_FullEP_rec.mp4` appears in the folder `Downloads`, and you successfully reconstructed the video.


Description of the scheme
=========================

We store the information on M sequences of length L. The sequences are indexed and protected with error correcting codes. The scheme performs the following steps:

- Multiply the information with a pseudorandom sequence. 

- Split the information into B*k information pieces of length L - log(M)-R, where R is the redundancy on a sequence level. The information of each piece will be stored on a sequence.

- For each block b = 1,...,B, add extra sequences by encoding each row separately by using a Reed-Solomon code over the extension field 2^m. This results in n-k extra sequences per blocks and thus in B*(n-k) extra sequences.

- Add a unique index to each sequence, and inner code each sequence.

- For the inner code we decide on a symbol size mi of 6 bit, for the outer code the symbol size mo is 14 bit. These values reflect considerations on accessible DNA sequence lengths and optimal code length.

We are not avoiding homopolymers. Due to the randomization of the data, long homopolymers and corresponding errors are very unlikely, and the error-correcting code will deal with those. 

Parameters of the code
======================

Here are the parameters of the code, and below are some examples to see how the code can be used to store information on DNA:


- l: length of the index in multiples of 6 bits (default choice is l=4)
- n: length of a block of the outer code (default is n=16383, n must be less or equal than 16383 (max outer code length equals 2^mo-1))
- k: number of information symbols of the outer code (default is k=10977, must obey k<=n)
- N: length of inner codeword (default is N=34; max inner code length is 2^mi-1)
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
Then the following are a subset of the choices that are possible due to the constraints given above:

| K  | nuss|  N | length of sequence in nucleotides  |
|:--:|:---:|:--:| :-----:|
| 11  | 3   | 14 | 42 |
| 18 | 6   | 21 | 63 |
| 25 | 9   | 28 | 84 |
| 32 | 12  | 35 | 105 |
| 39 | 15  | 42 | 126|
| 46 | 18  | 49 | 147|


Demo
====

In the following, we describe a few examples of how information can be encoded to DNA segments, and decode it back in order to recover the information. For all examples, the first step is to compile the program texttodna. Towards this goal, change to the folder ./simulate and execute the following command which compiles the code:

	`make texttodna` 

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

This example runs on a standard laptop (2017 MacBook Pro) in 1-5 minutes.

Example 2
---------

In our second example, we encode information on two blocks of maximum length n=16383, and choose k=12000. In addition we make our sequences shorter than in the previous example, by setting nuss=9. We set K = 25. Note that his satisfies the constraint K*mi = nuss*mo + l*mi, where mi=6, mo=14. Finally, we add more redundancy in the inner code by setting N = 28. This results in sequences of length N*ni/2 = 28*6/2 = 84.

1. Encoding:
	
	`./texttodna --k=12000 --nuss=9 --K=25 --N=28 --encode --input=../data/DNA_channel_statistics.pdf --output=../data/paper_encoded.txt`

2. Decoding:

	`./texttodna --k=12000 --nuss=9 --K=25 --N=28  --decode --numblocks=3 --input=../data/paper_encoded.txt --output=../data/paper_recovered.pdf`

Example 3
---------

The input can be a fastq file. In that case provide a fastq file as input as follows:
1. Decoding:

	`./texttodna --k=12000 --decode --numblocks=1 --fastq_infile=../data/reads.fastq  --output=../data/paper_recovered.pdf`

Use the flag --reverse if the reads are in reverse complement form.

Installation
============

The code is written in C++, and compilation requires installation of the boost library. The package has been tested on the following operating systems:

	- Ubuntu 16.04
	- Mac OS Mojave 10.14
    
The installation (compile time) is less than 2 minutes. The easiest way to run the code for Windows and Mac users is to install and run the code in a Docker container. 


Running the code in a Docker container - recommendened for Windows and Mac Users
--------------------------------------------------------------------------------
1. Install Docker desktop
2. Start an Ubuntu virtual machine, by runing the following commands in a terminal:
	
    For OSX:
	
	`docker run -it -v ~/Documents:/data ubuntu:16.04 bash`
	
    For Windows:

	`docker run -it -v c:/Users:/data ubuntu:16.04 bash`

This generates a docker container which has access to the folder ~/Documents/data on your local machine. Make sure the docker container has sufficient memory when decoding larger files - for the example here little memory is needed.

3. Install the required software in the Docker container:

        apt-get update
        apt-get install g++
        apt-get install make
        apt-get install git
        apt-get install libboost-all-dev

4. Download the code:

	`git clone https://github.com/reinhardh/dna_rs_coding/`

5. Compile the code:

        cd dna_rs_coding
        cd simulate
        make texttodna

Now the examples above can be run.

To work with your own data in the docker image:

    - copy files to an apropriate folder (e.g., cp data.zip /dna_rs_coding/data)
    - move to the folder /dna_rs_coding/simulate (e.g., by executing the command cd /dna_rs_coding/simulate) and run the code as in the examples to encode and decode files
    - copy results back from virtual machine to local folder on computer (e.g., with command: cp data_encoded.txt c:/Users/robert/)
    
Please note that once the virtual machine is ended (by command exit), all data generated within the machine is deleted, and only the data copied to local folders remains.


Installation of required software on Linux
------------------------------------------
    sudo apt-get install gcc
    sudo apt-get install make
    sudo apt-get install git
    sudo apt-get install libboost-all-dev
    
Follow steps 4-5 above to run and compile the code.

Installation of required software on macOS
------------------------------------------
	/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
	brew install boost
	sudo make –f Makefile

Follow steps 4-5 above to run and compile the code.
	
License
==========

All files are provided under the terms of the Apache License, Version 2.0, see the included file "apache_licence_20" for details.
