using namespace seqan;

template<typename TValue>
struct Minimizer
{
public:

    // Random, but static value for xor for hashes. Counteracts consecutive minimizers.
    // E.g., without it, the next minimizer after a poly-A region AAAAA would be most likely something like AAAAC.
    uint64_t const seed{0x8F3F73B5CF1C9ADE};
    // Shape for forward hashes
    Shape<TValue, SimpleShape> kmerShape;
    // Shape for hashes on reverse complement
    Shape<TValue, SimpleShape> revCompShape;
    // k-mer size
    uint16_t k{19};
    // window size
    uint32_t w{25};
    // start positions of minimizers
    std::vector<uint64_t> minBegin;
    // end positions of minimizers
    std::vector<uint64_t> minEnd;

    std::vector<uint32_t> coverage;
    std::vector<uint64_t> coverageBegin;
    std::vector<uint64_t> coverageEnd;


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

    inline void resize(uint16_t newKmerSize, uint32_t neww)
    {
        k = newKmerSize;
        w = neww;
        seqan::resize(kmerShape, k);
        seqan::resize(revCompShape, k);
    }

    std::vector<uint64_t> getHash(String<TValue> & text)
    {
        if (k > seqan::length(text))
            return std::vector<uint64_t> {};

        // Reverse complement without copying/modifying the original string
        typedef ModifiedString<ModifiedString<String<TValue>, ModComplementDna>, ModReverse> TRC;
        TRC revComp(text);

        uint64_t possible = seqan::length(text) > w ? seqan::length(text) - w + 1 : 1;
        uint32_t windowKmers = w - k + 1;

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
        for (uint32_t i = 0; i < windowKmers; ++i)
        {
            // Get smallest canonical k-mer
            uint64_t kmerHash = hashNext(it) ^ seed;
            uint64_t revcHash = revHashNext(rcit) ^ seed;
            uint64_t distance = std::distance(begin(text), it);
            if (kmerHash <= revcHash)
            {
                windowValues.push_back(std::make_tuple(kmerHash, distance, distance + k - 1));
            }
            else
            {
                windowValues.push_back(std::make_tuple(revcHash, distance, distance + k - 1));
            }
            ++it;
            ++rcit;
        }

        auto min = std::min_element(std::begin(windowValues), std::end(windowValues));
        kmerHashes.push_back(std::get<0>(*min));
        minBegin.push_back(std::get<1>(*min));
        minEnd.push_back(std::get<2>(*min));

        // For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
        // that results from the window shifting
        for (uint64_t i = 1; i < possible; ++i)
        {
            if (min == std::begin(windowValues))
            {
                windowValues.pop_front();
                min = std::min_element(std::begin(windowValues), std::end(windowValues));
            }
            else
                windowValues.pop_front();

            uint64_t kmerHash = hashNext(it) ^ seed;
            uint64_t revcHash = revHashNext(rcit) ^ seed;
            uint64_t distance = std::distance(begin(text), it);
            if (kmerHash <= revcHash)
            {
                windowValues.push_back(std::make_tuple(kmerHash, distance, distance + k - 1));
            }
            else
            {
                windowValues.push_back(std::make_tuple(revcHash, distance, distance + k - 1));
            }
            ++it;
            ++rcit;

            if (std::get<0>(windowValues.back()) < std::get<0>(*min))
                min = std::end(windowValues) - 1;

            kmerHashes.push_back(std::get<0>(*min));
            minBegin.push_back(std::get<1>(*min));
            minEnd.push_back(std::get<2>(*min));
        }
        kmerHashes.erase(std::unique(std::begin(kmerHashes), std::end(kmerHashes)), std::end(kmerHashes));
        return kmerHashes;
    }

    inline uint32_t get_threshold(uint32_t t, uint16_t e)
    {
        (void) t;
        uint32_t destroyed{0};
        minBegin.erase(std::unique(std::begin(minBegin), std::end(minBegin)), std::end(minBegin));
        uint32_t available{static_cast<uint32_t>(minBegin.size())};
        // destroyed += e;

        for (uint16_t i = 0; i < e; ++i)
        {
            get_coverage();
            auto max = std::max_element(coverage.begin(), coverage.end());
            if (i == e - 1)
                destroyed += *max;
            else
            {
                destroyed += *max;
                std::vector<uint64_t> newBegin;
                std::vector<uint64_t> newEnd;
                newBegin.reserve(available - destroyed);
                newEnd.reserve(available - destroyed);

                auto idx = std::distance(coverage.begin(), max);
                auto cb = coverageBegin[idx];
                auto ce =  coverageEnd[idx];
                for (uint64_t i = 0; i < minBegin.size(); ++i)
                {
                    auto mb = minBegin[i];
                    auto me = minEnd[i];
                    if ((mb >= cb && mb <= ce) || (me >= cb && me <= ce))
                        continue;
                    newBegin.push_back(mb);
                    newEnd.push_back(me);
                }
                minBegin = std::move(newBegin);
                minEnd = std::move(newEnd);
                coverageBegin.clear();
                coverageEnd.clear();
                coverage.clear();
            }
        }
        return destroyed > available ? 0 : available - destroyed;
    }

    inline void get_coverage()
    {
        uint64_t bIndex{1};
        uint64_t eIndex{0};
        coverageBegin.push_back(minBegin[0]);
        coverage.push_back(1);

        minBegin.erase(std::unique(std::begin(minBegin), std::end(minBegin)), std::end(minBegin));
        minEnd.erase(std::unique(std::begin(minEnd), std::end(minEnd)), std::end(minEnd));

        while ((bIndex < minBegin.size() ) || (eIndex < minEnd.size()))
        {
            uint64_t begin = bIndex < minBegin.size() ? minBegin[bIndex] : 0xFFFFFFFFFFFFFFFFULL;
            uint64_t end   = minEnd[eIndex];
            // Overlap
            if (begin < end)
            {
                coverageEnd.push_back(begin-1);
                coverageBegin.push_back(begin);
                coverage.push_back(coverage.back()+1);
                ++bIndex;
            }
            // Flatten consecutive positions, where one kmer ends and other one starts
            if (begin == end)
            {
                coverageEnd.push_back(begin-1);
                coverageBegin.push_back(begin);
                coverage.push_back(coverage.back()+1);
                while (minBegin[bIndex] == minEnd[eIndex])
                {
                    ++bIndex;
                    ++eIndex;
                }
                --eIndex;
            }
            // Kmer ends
            if (end < begin)
            {
                coverageEnd.push_back(end);
                coverageBegin.push_back(end+1);
                coverage.push_back(coverage.back()-1);
                ++eIndex;
            }
        }
        coverageBegin.pop_back();
        coverage.pop_back();
    }
};
