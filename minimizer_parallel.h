using namespace seqan;

struct Minimizer
{
public:

    // Random, but static value for xor for hashes. Counteracts consecutive minimizers.
    // E.g., without it, the next minimizer after a poly-A region AAAAA would be most likely something like AAAAC.
    uint64_t const seed{0x8F3F73B5CF1C9ADE};
    uint8_t threads{8};
    // k-mer size
    uint8_t k{19};
    // window size
    uint8_t w{25};
    // start positions of minimizers
    std::vector<uint64_t> minBegin;
    // end positions of minimizers
    std::vector<uint64_t> minEnd;

    inline void resize(uint8_t newKmerSize, uint8_t neww)
    {
        k = newKmerSize;
        w = neww;
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
        kmerHashes.resize(possible);
        minBegin.resize(possible);
        minEnd.resize(possible);

        uint64_t offset = possible / threads;

        std::vector<std::future<void>> tasks;

        for (uint8_t taskNo = 0; taskNo < threads; ++taskNo)
        {
            tasks.emplace_back(std::async([&, taskNo] {
                std::deque<std::tuple<uint64_t, uint64_t, uint64_t>> windowValues;
                uint64_t id = offset * taskNo;
                Shape<Dna, SimpleShape> kmerShape;
                Shape<Dna, SimpleShape> revCompShape;
                seqan::resize(kmerShape, k);
                seqan::resize(revCompShape, k);

                auto it = begin(text) + offset * taskNo;
                auto rcit = begin(revComp) + offset * taskNo;

                hashInit(kmerShape, it);
                hashInit(revCompShape, rcit);

                // Initialisation. We need to compute all hashes for the first window.
                for (uint8_t i = 0; i < windowKmers; ++i)
                {
                    // Get smallest canonical k-mer
                    uint64_t kmerHash = hashNext(kmerShape, it) ^ seed;
                    uint64_t revcHash = hashNext(revCompShape, rcit) ^ seed;
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
                kmerHashes[id] = std::get<0>(max);
                minBegin[id] = std::get<1>(max);
                minEnd[id] = std::get<2>(max);
                ++id;

                uint64_t threadPossible = taskNo != threads - 1 ? offset : length(text) - taskNo  * offset - w + 1;

                // For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
                // that results from the window shifting
                for (uint64_t i = 1; i < threadPossible; ++i)
                {
                    windowValues.pop_front();
                    uint64_t kmerHash = hashNext(kmerShape, it) ^ seed;
                    uint64_t revcHash = hashNext(revCompShape, rcit) ^ seed;
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
                    kmerHashes[id] = std::get<0>(max);
                    minBegin[id] = std::get<1>(max);
                    minEnd[id] = std::get<2>(max);
                    ++id;
                }
            }));
        }
        for (auto &&task : tasks)
        {
            task.get();
        }
        return kmerHashes;
    }
};
