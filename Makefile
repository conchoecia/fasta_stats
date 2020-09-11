CPPF=-g -O3 -W -Wall -Wno-unused-parameter -Wpointer-arith -Wshadow -Wundef -std=c++11
FSD=src

all: bin/fasta_stats

bin/fasta_stats: ${FSD}/fasta_stats.cc ${FSD}/itoa.cc ${FSD}/open_compressed.h
	mkdir -p bin
	g++ ${CPPF} -o bin/fasta_stats ${FSD}/*.cc
