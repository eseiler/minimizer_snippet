#include <chrono>
#include <cstdlib>
#include <seqan/seq_io.h>

int main(int argc, char const * argv[])
{
    using namespace seqan;

    uint8_t k = static_cast<uint8_t>(atoi(argv[2]));
    CharString id;
    DnaString seq;
    SeqFileIn seqFileIn;
    if (!open(seqFileIn, argv[1]))
        throw "Unable to open file.\n";
    readRecord(id, seq, seqFileIn);
    close(seqFileIn);

    Shape<Dna, SimpleShape> kmerShape;
    resize(kmerShape, k);
    auto start = std::chrono::high_resolution_clock::now();
    auto it = begin(seq);
    hashInit(kmerShape, it);

    uint64_t possible = length(seq) - k + 1;
    std::vector<uint64_t> kmerHashes;
    kmerHashes.reserve(possible);

    for (uint64_t i = 0; i < possible; ++i)
    {
        kmerHashes.push_back(hashNext(kmerShape, it));
        ++it;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    std::cerr << "The text of length " << length(seq) << " contains " << kmerHashes.size() << ' ' << (int)k << "-mers. Run time: " << duration << " ms.\n";
}
