#include <chrono>
#include <cstdlib>
#include <seqan/seq_io.h>

int main(int argc, char const * argv[])
{
    using namespace seqan;

    // k-mer size
    uint8_t k = static_cast<uint8_t>(atoi(argv[2]));

    // Read input
    CharString id;
    DnaString seq;
    SeqFileIn seqFileIn;
    if (!open(seqFileIn, argv[1]))
        throw "Unable to open file.\n";
    readRecord(id, seq, seqFileIn);
    close(seqFileIn);

    // Shape for hashing over Dna alphabet
    Shape<Dna, SimpleShape> kmerShape;
    resize(kmerShape, k);
    auto start = std::chrono::high_resolution_clock::now();
    auto it = begin(seq);

    // Initialise
    hashInit(kmerShape, it);

    uint64_t possible = length(seq) - k + 1;
    std::vector<uint64_t> kmerHashes;
    kmerHashes.reserve(possible);

    // The hash value is the arithmitic coding with base 4, i.e. A=0, C=1, G=2, T=3 and hash(TGCA) = 3*4^3 + 2*4^2 + 1*4^1 + 0*4^0
    // Internally the Shape remembers the left most character and the current hash. When we hash the next k-mer
    // (shift one to the right), we only need to remove the hash value of the leftmost character and add the value of
    // the new character. 
    for (uint64_t i = 0; i < possible; ++i)
    {
        kmerHashes.push_back(hashNext(kmerShape, it));
        ++it;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    std::cerr << "The text of length " << length(seq) << " contains " << kmerHashes.size() << ' ' << (int)k << "-mers. Run time: " << duration << " ms.\n";
}
