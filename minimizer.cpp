#include <chrono>
#include <seqan/seq_io.h>
#include <seqan/modifier.h>

#include "minimizer.h"

int main(int argc, char const * argv[])
{
    using namespace seqan;
    uint8_t k = static_cast<uint8_t>(atoi(argv[2]));
    uint8_t w = static_cast<uint8_t>(atoi(argv[3]));
    Minimizer minimizer;
    minimizer.resize(k, w);

    CharString id;
    DnaString seq;
    SeqFileIn seqFileIn;
    if (!open(seqFileIn, argv[1]))
        throw "Unable to open file.\n";
    readRecord(id, seq, seqFileIn);
    close(seqFileIn);

    auto start = std::chrono::high_resolution_clock::now();
    auto hashvalues = minimizer.getHash(seq);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    auto uniqueBegins = minimizer.minBegin;
    uniqueBegins.erase(unique(uniqueBegins.begin(), uniqueBegins.end()), uniqueBegins.end());
    auto uniqueEnds = minimizer.minEnd;
    uniqueEnds.erase(unique(uniqueEnds.begin(), uniqueEnds.end()), uniqueEnds.end());

    std::cerr << "The text of length " << length(seq) << " contains " << uniqueBegins.size() << " distinct minimizers(" << (int)k << ',' << (int)w <<"). Run time: " << duration << " ms.\n";

    // for (uint64_t i = 0; i < 30/*uniqueBegins.size()*/; ++i)
    // {
    //     std::cerr << hashvalues[i] << '\t' << uniqueBegins[i] << '\t' << uniqueEnds[i] << '\n';
    // }
}
