
texttodna: texttodna.cpp
	g++ -O3 texttodna.cpp -o texttodna -std=c++11 -lboost_program_options
	#g++ -O3 texttodna.cpp -o texttodna -std=c++11 -lboost_program_options -L/Users/reinhardheckel/include/boost_1_55_0/stage/lib

test: test1

test1:
	./texttodna --n=12472 --k=9000 --encode --input=../data/data.zip --output=../data/data_encoded.txt
	./texttodna --disturb --input=../data/data_encoded.txt --output=../data/data_drawnseg.txt
	./texttodna --decode --n=12472 --k=9000 --numblocks=1 --input=../data/data_drawnseg.txt --output=../data/data_rec.zip
	diff ../data/data.zip ../data/data_rec.zip
