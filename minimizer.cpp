#include <seqan/seq_io.h>
#include <seqan/modifier.h>

#include "minimizer.h"

int main(int argc, char const * argv[])
{
    using namespace seqan;
    Minimizer minimizer;
    minimizer.resize(19, 25);

    CharString id;
    DnaString seq;
    SeqFileIn seqFileIn;
    if (!open(seqFileIn, argv[1]))
        throw "Unable to open file.\n";
    readRecord(id, seq, seqFileIn);
    close(seqFileIn);

    auto hashvalues = minimizer.getHash(seq);

    for (uint64_t i = 0; i < 30/*minimizer.minEnd.size()*/; ++i)
    {
        std::cerr << hashvalues[i] << '\t' << minimizer.minBegin[i] << '\t' << minimizer.minEnd[i] << '\n';
    }
}
