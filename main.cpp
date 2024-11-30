#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <set>
#include <math.h>
using namespace std;

string fasta1_path = "test_sets/GCA_000717535.1_ASM71753v1_genomic.fna";
string fasta2_path = "test_sets/GCA_000835235.1_ASM83523v1_genomic.fna";

// read fasta file
string parse_fasta(const string& file1) {
    // open the fasta file
    ifstream file(file1);
    
    // check if file is opened successfully
    if (!file) {
        cerr << "Error: Could not open file " << file1 << endl;
    }

    // parse the fasta file
    string fasta;
    string line;
    while (getline(file, line)) {
        // add all nucleotide lines to each other (change later for files with multiple sequences)
        if (line[0] != '>') {
            fasta += line;
        }
    }

    // close file
    file.close();

    return fasta;
}

// create array, split fasta into kmers
vector<string> get_kmers(const string& fasta, size_t K) {

    // calculate sequence length
    const size_t N = fasta.size();
    
    // check if the kmer size is compatible with the fasta sequence
    if (K < 0 || N < K) {
        cerr << "Invalid value of K for the given sequence length N." << endl;
    }

    // declare kmer array
    vector<string> kmers;

    // split sequence by kmers
    for (size_t i = 0; i <= (N - K); ++i) {
        kmers.push_back(fasta.substr(i, K));
    }
    
    return kmers;
}

char comp(char nt) {
    switch(nt) {
        case 'a':
        case 'A':
            return 'T';
        case 't':
        case 'T':
            return 'A';
        case 'c':
        case 'C':
            return 'G';
        case 'g':
        case 'G':
            return 'C';
        default:
            cerr << "Invalid nucleotide " << nt << "." << endl;
            return 'N';
    }
}

string rev_comp(const string& kmer, size_t K) {
    
    // initialize the reverse complement kmer
    string rc_kmer = kmer;
    
    // reverse the kmer
    for (size_t i = 0; i < K ; ++i) {
        rc_kmer[i] = comp(kmer[K - i - 1]);
    }

    return rc_kmer;
}

// murmurhash3, implementation adapted from https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp
// it does a bunch of random bitwise operations and mathematical operations to create a hash from a key
// #include <cstring>  contains std::hash, which could have also been used
inline __attribute__((always_inline)) uint32_t getblock32 (const uint32_t * p, int i)
{
  return p[i];
}

inline __attribute__((always_inline)) uint32_t rotl32 (uint32_t x, int8_t r)
{
  return (x << r) | (x >> (32 - r));
}

inline __attribute__((always_inline)) uint32_t fmix32 (uint32_t h)
{
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;

  return h;
}

uint32_t MurmurHash3_x86_32(const char * key, size_t len, uint32_t seed)
{
    const uint8_t * data = (const uint8_t*)key;

    const int nblocks = len / 4;

    uint32_t h1 = seed;

    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;

    //----------
    // body

    const uint32_t * blocks = (const uint32_t *)(data + nblocks*4);

    for(int i = -nblocks; i; i++)
    {
        uint32_t k1 = getblock32(blocks,i);

        k1 *= c1;
        k1 = rotl32(k1,15);
        k1 *= c2;
        
        h1 ^= k1;
        h1 = rotl32(h1,13); 
        h1 = h1*5+0xe6546b64;
    }

    //----------
    // tail

    const uint8_t * tail = (const uint8_t*)(data + nblocks*4);

    uint32_t k1 = 0;

    switch(len & 3)
    {
    case 3: k1 ^= tail[2] << 16;
    case 2: k1 ^= tail[1] << 8;
    case 1: k1 ^= tail[0];
            k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1;
    };

    //----------
    // finalization

    h1 ^= len;

    h1 = fmix32(h1);

    return h1;
} 

void create_bottom_sketch(vector<string>& kmers, vector<uint32_t>& bottom_sketch, size_t sketch_size, size_t K) {
    
    // hash kmers
    uint32_t seed = 0;
    for (string& kmer : kmers) {
        
        //cout << "kmer" << endl;
        //cout << kmer << endl;

        // initialize hash
        uint32_t kmer_hash;

        // get the reverse complement kmer
        string rc_kmer = rev_comp(kmer, K);

        // hash the lexicographically smaller of the two
        // the lexicographically smallest of a kmer and its refcomp is always the same, so this saves a lot of time
        if (kmer <= rc_kmer) {

            // the variable kmer is a string object, to get a pointer to its content, use .c_str()
            const char * pkmer = kmer.c_str();
            kmer_hash = MurmurHash3_x86_32(pkmer, K, seed);
        }
        else {

            // the variable rc_kmer is a string object, to get a pointer to its content, use .c_str()
            const char * prc_kmer = rc_kmer.c_str();
            kmer_hash = MurmurHash3_x86_32(prc_kmer, K, seed);
        }

        // only store the kmer if it is smaller than the largest hash in the sketch
        // bottom_sketch contains the sorted smallest hashes, smallest hash is the first value
        if (kmer_hash < bottom_sketch.back()) {

            // find the insertion point for the new hash
            auto pos = lower_bound(bottom_sketch.begin(), bottom_sketch.end(), kmer_hash);

            // only insert if the value is not already present
            if (pos == bottom_sketch.end() || *pos != kmer_hash) {

                // insert the new hash at the correct position
                bottom_sketch.insert(pos, kmer_hash);

                // keep the vector size fixed by removing the largest hash
                bottom_sketch.pop_back();
            }
        }
    }
}

int main() {

    // get fasta sequences from file
    cout << "Reading fasta files" << endl;
    string fasta1 = parse_fasta(fasta1_path);
    string fasta2 = parse_fasta(fasta2_path);

    // get kmers from sequences
    cout << "Extracting kmers" << endl;
    const size_t Kmer_size = 21;
    vector<string> kmers1 = get_kmers(fasta1, Kmer_size);
    vector<string> kmers2 = get_kmers(fasta2, Kmer_size);

    // initialize bottom sketch
    const size_t sketch_size = 1000;
    vector<uint32_t> bottom_sketch1(sketch_size, numeric_limits<uint32_t>::max());
    vector<uint32_t> bottom_sketch2(sketch_size, numeric_limits<uint32_t>::max());

    // create bottom sketch
    cout << "Creating sketches" << endl;
    create_bottom_sketch(kmers1, bottom_sketch1, sketch_size, Kmer_size);
    create_bottom_sketch(kmers2, bottom_sketch2, sketch_size, Kmer_size);

    // calculate Jaccard index

    // convert vectors to sets
    std::set<uint32_t> set1(bottom_sketch1.begin(), bottom_sketch1.end());
    std::set<uint32_t> set2(bottom_sketch2.begin(), bottom_sketch2.end());

    // find the intersection (common elements)
    std::set<uint32_t> intersection;
    std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(),
                          std::inserter(intersection, intersection.begin()));

    cout << intersection.size() << endl;

    // find the union (all unique elements from both sets)
    std::set<uint32_t> union_set;
    std::set_union(set1.begin(), set1.end(), set2.begin(), set2.end(),
                   std::inserter(union_set, union_set.begin()));

    cout << union_set.size() << endl;

    // jaccard index = size of intersection / size of union
    double jaccard_index = static_cast<double>(intersection.size()) / union_set.size();

    cout << "Jaccard Index: " << jaccard_index << endl;

    // convert the jaccard index to a mutation rate (ANI = 1 - mutation rate)
    double d = (-1.0/Kmer_size)*log((2*jaccard_index)/(1 + jaccard_index));

    cout << "Mutation rate: " << d << endl;

    return 0;
}

