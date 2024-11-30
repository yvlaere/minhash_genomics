# minhash_genomics
Tool for similarity estimation of genomic data using the minhash algorithm

# Methodology
1. Load fasta files
2. Split sequences into kmers
3. Hash kmers
   - using murmurhash3 (https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp)
   - hashing only the lexicographically smaller of a kmer and its reverse complement
4. Store the smallest s kmer hashes in a bottom sketch
5. Calculate the jaccard index and mutation distance using the bottom sketches

# Inspiration 
mash paper:
https://doi.org/10.1186/s13059-016-0997-x
