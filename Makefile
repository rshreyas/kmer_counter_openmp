CXX=g++
CPPFLAGS= -std=gnu++11 -g -O3 -fopenmp
LDLIBS =
OBJECTS = fasta.o
FASTA_READER_PATH = ./fasta_reader/

all: kmer_count

fasta.o: $(FASTA_READER_PATH)fasta.c
	g++ -c -O3 -fomit-frame-pointer $(FASTA_READER_PATH)fasta.c

kmer_count: $(OBJECTS)
	$(CXX) $(CPPFLAGS) kmer_counter.cpp -o kmer_count $(LDLIBS) $(OBJECTS)

clean:
	rm -f *.o *~ kmer_count
