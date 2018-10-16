#include <chrono>
#include <seqan/seq_io.h>
#include <seqan/modifier.h>

#include "minimizer.h"

struct separators : std::numpunct<char> {
   char do_thousands_sep() const { return ' '; }

   std::string do_grouping() const { return "\3"; }
};

int main(int argc, char const * argv[])
{
    using namespace seqan;

    std::cerr.imbue(std::locale(std::locale(), new separators));
    // k-mer size
    uint8_t k = static_cast<uint8_t>(atoi(argv[2]));
    // window size
    uint8_t w = static_cast<uint8_t>(atoi(argv[3]));

    // Read Input
    CharString id;
    DnaString seq;
    SeqFileIn seqFileIn;
    if (!open(seqFileIn, argv[1]))
        throw "Unable to open file.\n";

    uint64_t distinctMinimizers{0};
    double duration{0.0};
    uint64_t textLength{0};

    while (!atEnd(seqFileIn))
    {
        Minimizer minimizer;
        minimizer.resize(k, w);
        readRecord(id, seq, seqFileIn);
        auto start = std::chrono::high_resolution_clock::now();
        auto volatile hashvalue = minimizer.getHash(seq);
        auto end = std::chrono::high_resolution_clock::now();
        duration += std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
        auto uniqueBegins = minimizer.minBegin;
        uniqueBegins.erase(unique(uniqueBegins.begin(), uniqueBegins.end()), uniqueBegins.end());
        distinctMinimizers += uniqueBegins.size();
        textLength += length(seq);
    }

    close(seqFileIn);

    std::cerr << "The text of length " << textLength << " contains " << distinctMinimizers << " distinct minimizers(" << (int)k << ',' << (int)w <<"). Run time: " << duration << " ms.\n";
}
