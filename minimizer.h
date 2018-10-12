using namespace seqan;

struct Minimizer
{
public:

    uint64_t const seed{0x8F3F73B5CF1C9ADE};
    Shape<Dna, SimpleShape> kmerShape;
    Shape<Dna, SimpleShape> revCompShape;
    uint8_t k{19};
    uint8_t w{25};
    std::vector<uint64_t> minBegin;
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

        typedef ModifiedString<ModifiedString<DnaString, ModComplementDna>, ModReverse> TRC;
        TRC revComp(text);

        uint64_t possible = seqan::length(text) > w ? seqan::length(text) - w + 1 : 1;
        uint8_t windowKmers = w - k + 1;

        std::vector<uint64_t> kmerHashes;
        std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> windowValues;
        kmerHashes.reserve(possible);
        minBegin.reserve(possible);
        minEnd.reserve(possible);
        windowValues.reserve(windowKmers);

        auto it = begin(text);
        auto rcit = begin(revComp);
        hashInit(it);
        revHashInit(rcit);

        for (uint8_t i = 0; i < windowKmers; ++i)
        {
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

        for (uint64_t i = 1; i < possible; ++i)
        {
            windowValues.erase(std::begin(windowValues));
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
