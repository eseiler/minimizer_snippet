##### Setup

```console
git clone https://github.com/eseiler/minimizer_snippet.git
git clone https://github.com/seqan/seqan.git
mkdir build
cd build
cmake ../minimizer_snippet -DCMAKE_CXX_FLAGS="-O3" -DCMAKE_CXX_COMPILER=g++ -DCMAKE_PREFIX_PATH=$(pwd)/../seqan/util/cmake/ -DSEQAN_INCLUDE_PATH=$(pwd)/../seqan/include/
make minimizer
./minimizer ../minimizer_snippet/test.fasta 19 25
make full
./full ../minimizer_snippet/test.fasta 25

```

##### Example Output

```console
The text of length 1000000 contains 246097 distinct minimizers(19,25). Run time: 80 ms.
The text of length 1000000 contains 999976 25-mers. Run time: 13 ms.
```
