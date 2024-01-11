TREW : Telomeric Repeat motif Estimation tool with Whole-genome sequencing
==========================================================================

TREW is tool that can identify Telomeric repeat motif (TRM) with short-read sequencing data. This tool looks for repeated sequences in a single read to find candidates for TRMs, iterating through them to finally find a TRM.

[![Build](https://github.com/Chemical118/TREW/actions/workflows/build.yaml/badge.svg)](https://github.com/Chemical118/TREW/actions/workflows/build.yaml)[![CI](https://github.com/Chemical118/TREW/actions/workflows/ci.yaml/badge.svg)](https://github.com/Chemical118/TREW/actions/workflows/ci.yaml)[![codecov](https://codecov.io/gh/Chemical118/TREW/graph/badge.svg?token=WRDCVZUAWH)](https://codecov.io/gh/Chemical118/TREW)

### Install

You can download a binary from release or build from source.

[Windows (x86_64)](https://github.com/Chemical118/TREW/releases/latest/download/trew-windows-x86_64.tar.gz)
[Linux (x86_64)](https://github.com/Chemical118/TREW/releases/latest/download/trew-linux-x86_64.tar.gz)
[MacOS (x86_64, arm64)](https://github.com/Chemical118/TREW/releases/latest/download/trew-macos-universal.tar.gz)

### Install from source

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

The following is a recommended command line to run TREW.

```sh
trew short 5 32 <short_read_data1.fastq.gz> <short_read_data2.fastq>... -t <number of threads>

trew long 5 32 <long_read_data1.fastq.gz> <long_read_data2.fastq>... -t <number of threads>
```

### Output

Short-read sequencing

```
><short_read_data1.fastq.gz>
<length of repeat>,<repeat sequence>,<number of repeat>
...

><short_read_data2.fastq.gz>
<length of repeat>,<repeat sequence>,<number of repeat>
...
```

Long-read sequencing

```
><long_read_data1.fastq.gz>
<length of repeat>,<repeat sequence>,<number of repeat>,<number of reverse repeat>,<number of pure repeat>
<length of palindromic repeat>,<repeat sequence>,<number of repeat>,_,<number of pure repeat>
...

><long_read_data2.fastq.gz>
<length of repeat>,<repeat sequence>,<number of repeat>,<number of reverse repeat>,<number of pure repeat>
<length of palindromic repeat>,<repeat sequence>,<number of repeat>,_,<number of pure repeat>
...
```

In long-read sequencing data, we count with repeat orientation to accurately recognize TRMs. If the number of repeats and number of reverse repeats are similar, this repeat is very unlikely to be a TRM.

### Author

Hyunwoo Ryu <wowo0118@korea.ac.kr>

*Special thanks to*
Jiho Choi <sdatoli@korea.ac.kr>
Kyungmo Ku <kyungmoku7141@gmail.com>
