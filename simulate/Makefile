
BOOST_LIB = /Users/reinhardheckel/include/boost_1_55_0/stage/lib

texttodna: texttodna.cpp
	g++ -O3 texttodna.cpp -o texttodna -std=c++11 -lboost_program_options -L$(BOOST_LIB)
	#g++ -O3 texttodna.cpp -o texttodna -std=c++11 -lboost_program_options -L/Users/reinhardheckel/include/boost_1_55_0/stage/lib

tests: test1

test1:
	./texttodna --n=12472 --k=9000 --encode --input=../data/Linda_reformat.zip --output=../data/linda_encoded.txt
	./texttodna --disturb --input=../data/linda_encoded.txt --output=../data/linda_drawnseg.txt
	./texttodna --decode --n=12472 --k=9000 --numblocks=1 --input=../data/linda_drawnseg.txt --output=../data/linda_rec.zip
	diff ../data/Linda_reformat.zip ../data/linda_rec.zip
