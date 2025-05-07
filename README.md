## Overview

This repository provides a learning implementation of similarity estimation between genomic sequences using MinHash, inspired by Mash (https://github.com/marbl/Mash) (https://doi.org/10.1186/s13059-016-0997-x) .
## Prerequisites
- C++ compiler (e.g. g++, clang)

## Methodology and Algorithms

### 1. K-mer Extraction

Each DNA sequence is split into overlapping k-mers. For a sequence of length _L_, there are _(L − k + 1)_ k-mers. Using k-mers greatly reduces large datasets to a manageable feature set. The reverse complement of each k-mer is created and only the lexicographically smallest of the k-mer and its reverse complement is kept, ensuring forward and reverse sequences share the same k-mers.

### 2. Hashing with MurmurHash3

Each k-mer is hashed to a 32-bit integer using MurmurHash3, a non-cryptographic but fast and well-known hash function. The code for this is adapted from https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp

### 3. Bottom-s Sketch Construction

To avoid storing all k-mers, a bottom-s MinHash sketch (also called “signature”) is created. This makes the comparison of multiple datasets much faster:

1. Initialize an empty array of size _s_.
    
2. Only insert a new hash if its smaller than the largest hash in the sketch. The sketch is sorted for fast evaluation.
    
The result is the _s_ smallest hash values. Two sequences’ sketches can be compared rapidly in _O(s)_ time, independent of original length.

### 4. Jaccard Similarity Estimation

Given two bottom-s sketches $S(A)$ and $S(B)$​, the Jaccard index can be estimated as

$$J(A,B) = \frac{|A \cap B|}{|A \cup B|} \approx \frac{|S(A\cup B) \cap S(A) \cap S(B)|}{|S(A \cup B)|}$$

Empirically, the fraction of matching hash values in the sketches approximates the true Jaccard index.

### 5. Mutation Distance Approximation

Assuming a simple Poisson model of substitutions, Mash shows that an estimated Jaccard $\hat{J}$ yields an average nucleotide identity (ANI)–derived distance:

$$d=−\frac{1}{k}ln⁡(\frac{2\hat{J}}{1 + \hat{J}})$$

This converts the estimated jaccard index into an approximate mutation rate per basepair, which $\approx 1 - ANI$.

## References

1. Mash: Ondov, B.D., Treangen, T.J., Melsted, P. _et al._ Mash: fast genome and metagenome distance estimation using MinHash. _Genome Biol_ **17**, 132 (2016). https://doi.org/10.1186/s13059-016-0997-x
    
2. Appleby, A. MurmurHash3 (https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp)
