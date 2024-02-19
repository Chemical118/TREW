TREW : Telomeric Repeat motif Estimation tool with Whole-genome sequencing data
===============================================================================

TREW is tool that can identify Telomeric repeat motif (TRM) with any whole-genome sequencing data. This tool looks for repeated sequences in a single read to find candidates for TRMs, iterating through them to finally find a TRM.

[![Build](https://github.com/Chemical118/TREW/actions/workflows/build.yaml/badge.svg)](https://github.com/Chemical118/TREW/actions/workflows/build.yaml)[![CI](https://github.com/Chemical118/TREW/actions/workflows/ci.yaml/badge.svg)](https://github.com/Chemical118/TREW/actions/workflows/ci.yaml)[![codecov](https://codecov.io/gh/Chemical118/TREW/graph/badge.svg?token=WRDCVZUAWH)](https://codecov.io/gh/Chemical118/TREW)

### Install

You can install `TREW` by downloading a binary from the release or building from the source.

[Windows (x86_64)](https://github.com/Chemical118/TREW/releases/latest/download/trew-windows-x86_64.tar.gz)  
[Linux (x86_64)](https://github.com/Chemical118/TREW/releases/latest/download/trew-linux-x86_64.tar.gz)  
[MacOS (x86_64, arm64)](https://github.com/Chemical118/TREW/releases/latest/download/trew-macos-universal.tar.gz)  

### Install from source

> We found that Intel® oneAPI DPC++/C++ Compiler has the potential to create ~20% faster programs on Intel® CPUs, especially, using multi-threading.

Windows (Visual Studio)

```sh
git clone https://github.com/Chemical118/TREW.git
cd TREW

git clone https://github.com/Microsoft/vcpkg.git
.\vcpkg\bootstrap-vcpkg.bat

mkdir build
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=vcpkg\scripts\buildsystems\vcpkg.cmake
cmake --build build --config Release

.\build\Release\trew_test.exe -- test\test.fastq.gz test\test.fastq test\test_long.fastq.gz test\test_long.fastq
.\build\Release\trew.exe -h
```

Unix

```sh
git clone https://github.com/Chemical118/TREW.git
cd TREW

git clone https://github.com/Microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.sh

mkdir build
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build build

./build/trew_test.exe -- test/test.fastq.gz test/test.fastq test/test_long.fastq.gz test/test_long.fastq
./build/trew -h
```

#### Install from source with Intel® oneAPI DPC++/C++ Compiler

Install compiler at [Intel® oneAPI DPC++/C++ Compiler website](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html#gs.4ooj5x) with default destination.

Windows
```sh
git clone https://github.com/Chemical118/TREW.git
cd TREW

git clone https://github.com/Microsoft/vcpkg.git
.\vcpkg\bootstrap-vcpkg.bat

mkdir build
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\icx.exe -DCMAKE_CXX_COMPILER=C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\icx.exe -DCMAKE_TOOLCHAIN_FILE=vcpkg\scripts\buildsystems\vcpkg.cmake
cmake --build build --config Release

.\build\Release\trew_test.exe -- test\test.fastq.gz test\test.fastq test\test_long.fastq.gz test\test_long.fastq
.\build\Release\trew.exe -h
```

Linux
```sh
git clone https://github.com/Chemical118/TREW.git
cd TREW

git clone https://github.com/Microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.sh

mkdir build
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=$HOME/intel/oneapi/compiler/latest/bin/icx -DCMAKE_CXX_COMPILER=$HOME/intel/oneapi/compiler/latest/bin/icpx -DCMAKE_TOOLCHAIN_FILE=vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build build

./build/trew_test.exe -- test/test.fastq.gz test/test.fastq test/test_long.fastq.gz test/test_long.fastq
./build/trew -h
```
### Quick Start

Short-read sequencing

```sh
trew short MIN_MER MAX_MER <short_read_data1.fastq.gz> <short_read_data2.fastq>... -t <number of threads> 
```

Long-read sequencing

```sh
trew long MIN_MER MAX_MER <long_read_data1.fastq.gz> <long_read_data2.fastq>... -t <number of threads> 
```

MIN_MER : minimum length of sequence to find telomere [MIN_MER >= 3]  
MAX_MER : maximum length of sequence to find telomere [MAX_MER <= 64]

> Note that to get the correct putative TRM, you must run the program on the same species sequencing data.

The following is a recommended command line to run TREW.

```sh
trew short 5 32 <short_read_species1_data1.fastq.gz> <short_read_species1_data2.fastq>... -t <number of threads>
trew short 5 32 <short_read_species2_data1.fastq.gz> <short_read_species2_data2.fastq>... -t <number of threads>

trew long 5 32 <long_read_species1_data1.fastq.gz> <long_read_species1_data2.fastq>... -t <number of threads>
trew long 5 32 <long_read_species2_data1.fastq.gz> <long_read_species2_data2.fastq>... -t <number of threads>
```

### Output

```
>H:<read_data1.fastq.gz>
<length of repeat>,<repeat sequence>,<number of repeat>,<number of reverse repeat>,<number of pure repeat>
<length of palindromic repeat>,<repeat sequence>,<number of repeat>,-1,<number of pure repeat>
...

>L:<read_data1.fastq.gz>
<length of repeat>,<repeat sequence>,<number of repeat>,<number of reverse repeat>,<number of pure repeat>
<length of palindromic repeat>,<repeat sequence>,<number of repeat>,-1,<number of pure repeat>
...

>H:<read_data2.fastq.gz>
<length of repeat>,<repeat sequence>,<number of repeat>,<number of reverse repeat>,<number of pure repeat>
<length of palindromic repeat>,<repeat sequence>,<number of repeat>,-1,<number of pure repeat>
...

>L:<read_data2.fastq.gz>
<length of repeat>,<repeat sequence>,<number of repeat>,<number of reverse repeat>,<number of pure repeat>
<length of palindromic repeat>,<repeat sequence>,<number of repeat>,-1,<number of pure repeat>
...

>Putative_TRM
<putative telomeric repeat motif>,<score>
...
or
NO_PUTATIVE_TRM,-1
```

`:H` means _high_ baseline search, this raw result might show homogeneous repeats.
Also, `:L` means _low_ baseline search, this raw result might show heterogeneous repeats.
The value of score in `Putative_TRM` can range from 1 to 10, but we recommend checking all putative TRMs regardless of the value of score.

### Author

Hyunwoo Ryu <wowo0118@korea.ac.kr>

*Special thanks to*  
Jiho Choi <sdatoli@korea.ac.kr>  
Kyungmo Ku <kyungmoku7141@gmail.com>
