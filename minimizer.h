using namespace seqan;

struct Minimizer
{
public:

    // Random, but static value for xor for hashes. Counteracts consecutive minimizers.
    // E.g., without it, the next minimizer after a poly-A region AAAAA would be most likely something like AAAAC.
    uint64_t const seed{0x8F3F73B5CF1C9ADE};
    // Shape for forward hashes
    Shape<Dna, SimpleShape> kmerShape;
    // Shape for hashes on reverse complement
    Shape<Dna, SimpleShape> revCompShape;
    // k-mer size
    uint8_t k{19};
    // window size
    uint8_t w{25};
    // start positions of minimizers
    std::vector<uint64_t> minBegin;
    // end positions of minimizers
    std::vector<uint64_t> minEnd;

    template<typename TIt>
    inline void hashInit(TIt it)
    {
        seqan::hashInit(kmerShape, it);
    }

    template<typename TIt>
    inline auto hashNext(TIt it)
    {
        return seqan::hashNext(kmerShape, it);
    }

    template<typename TIt>
    inline void revHashInit(TIt it)
    {
        seqan::hashInit(revCompShape, it);
    }

    template<typename TIt>
    inline auto revHashNext(TIt it)
    {
        return seqan::hashNext(revCompShape, it);
    }

    inline auto length()
    {
        return seqan::length(kmerShape);
    }

    inline void resize(uint8_t newKmerSize, uint8_t neww)
    {
        k = newKmerSize;
        w = neww;
        seqan::resize(kmerShape, k);
        seqan::resize(revCompShape, k);
    }

    std::vector<uint64_t> getHash(DnaString & text)
    {
        if (k > seqan::length(text))
            return std::vector<uint64_t> {};

        // Reverse complement without copying/modifying the original string
        typedef ModifiedString<ModifiedString<DnaString, ModComplementDna>, ModReverse> TRC;
        TRC revComp(text);

        uint64_t possible = seqan::length(text) > w ? seqan::length(text) - w + 1 : 1;
        uint8_t windowKmers = w - k + 1;

        std::vector<uint64_t> kmerHashes;
        // Stores hash, begin and end for all k-mers in the window
        std::deque<std::tuple<uint64_t, uint64_t, uint64_t>> windowValues;
        kmerHashes.reserve(possible);
        minBegin.reserve(possible);
        minEnd.reserve(possible);

        auto it = begin(text);
        auto rcit = begin(revComp);
        hashInit(it);
        revHashInit(rcit);

        // Initialisation. We need to compute all hashes for the first window.
        for (uint8_t i = 0; i < windowKmers; ++i)
        {
            // Get smallest canonical k-mer
            uint64_t kmerHash = hashNext(it) ^ seed;
            uint64_t revcHash = revHashNext(rcit) ^ seed;
            if (kmerHash <= revcHash)
            {
                uint64_t distance = std::distance(begin(text), it);
                windowValues.push_back(std::make_tuple(kmerHash, distance, distance + k - 1));
            }
            else
            {
                uint64_t distance = std::distance(begin(revComp), rcit);
                windowValues.push_back(std::make_tuple(revcHash, distance, distance + k - 1));
            }
            ++it;
            ++rcit;
        }

        auto max = *std::min_element(std::begin(windowValues), std::end(windowValues));
        kmerHashes.push_back(std::get<0>(max));
        minBegin.push_back(std::get<1>(max));
        minEnd.push_back(std::get<2>(max));

        // For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
        // that results from the window shifting
        for (uint64_t i = 1; i < possible; ++i)
        {
            windowValues.pop_front();
            uint64_t kmerHash = hashNext(it) ^ seed;
            uint64_t revcHash = revHashNext(rcit) ^ seed;
            if (kmerHash <= revcHash)
            {
                uint64_t distance = std::distance(begin(text), it);
                windowValues.push_back(std::make_tuple(kmerHash, distance, distance + k - 1));
            }
            else
            {
                uint64_t distance = std::distance(begin(revComp), rcit);
                windowValues.push_back(std::make_tuple(revcHash, distance, distance + k - 1));
            }
            ++it;
            ++rcit;

            auto max = *std::min_element(std::begin(windowValues), std::end(windowValues));
            kmerHashes.push_back(std::get<0>(max));
            minBegin.push_back(std::get<1>(max));
            minEnd.push_back(std::get<2>(max));
        }

        return kmerHashes;
    }
};
