# dedup

`dedup` is command-line tool for removing duplicate sequences from FASTA files. It can identify and remove sequences that are either identical (100% identity) or highly similar (e.g., 95% identity).


## How to Compile

To compile the tool, you need a C++ compiler that supports C++17 and has support for SSE4.2 intrinsics (most modern x86-64 CPUs do).

**On Linux or macOS:**

Use `g++` (version 8 or newer) or `clang++`.

```bash
make
```

**On Windows:**

It is recommended to use the MinGW-w64 toolchain to get a `g++` compiler.

```bash
# Make sure MinGW is in your PATH
make
```

This will produce an executable named `dedup` (or `dedup.exe` on Windows).

To run smoke tests:

```bash
make test
```

## Usage

The tool is run from the command line with the following options:

```
./dedup -i <input.fasta> -o <output.fasta> [options]
```

### Options

| Flag | Argument | Description | Default |
| :--- | :--- | :--- | :--- |
| `-i` | `FILE` | **Required.** Input FASTA file containing the sequences. | |
| `-o` | `FILE` | **Required.** Output FASTA file where the unique sequences will be written. | |
| `-p` | `INT` | Identity percentage threshold. Accepted values: `100`, `95`, `90`, `85`, `80`. | `100` |
| `-k` | `INT` | K-mer size for indexing when using a percentage threshold < 100. | `12` |
| `-t` | `INT` | Number of threads to use. | Hardware concurrency |
| `-h` | | Show the help message. | |

### Examples

**1. Remove exact duplicate sequences:**

This is the default mode. It will keep the first occurrence of each unique sequence.

```bash
./dedup -i my_sequences.fasta -o unique_sequences.fasta
```

**2. Remove sequences that are 95% identical or more:**

This will remove sequences that are highly similar to a sequence that appeared earlier in the file.

```bash
./dedup -i my_sequences.fasta -o unique_95_percent.fasta -p 95
```

**3. Use 8 threads and a k-mer size of 10 for 90% identity:**

```bash
./dedup -i large_dataset.fasta -o unique_90_percent.fasta -p 90 -k 10 -t 8
```
