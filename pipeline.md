# Pipeline

open questions:

- wget fastq files from GitHub,
- all code blocks are launched from $HOME and not from the latest
  folder?


## Google Colab

Let's explore the environment:

``` code
%%shell

date
whoami
```

We are root! Maximal clearance level.

Check version numbers:

``` code
%%shell

uname -a
bash --version
git --version
gcc --version
python --version
R --version
```

What about hardware resources?

``` code
%%shell

df -h
cat /proc/cpuinfo
cat /proc/meminfo
```

It seems that Google colab instances are virtual x86-64 machines with
2 CPU-cores, 16 GB of RAM and 80 GB of storage space. The operating
system is Ubuntu LTS 20.04. LTS stands for long-term support, the most
recent LTS at the time of writing is 22.04. Version 20.04 will be
supported until April 2025. Its `gcc` version is a bit old (9.2), and
might not allow to compile code based on very recent standards (e.g.;
`mumu` wich uses C++20 features).

I will assume that you also have [python](https://www.python.org/)
(version 3.5 or more), [R](https://cran.r-project.org/) (version 3.5
or more), and [bash](https://www.gnu.org/software/bash/) (version 4 or
more).

Git and compilation tools are already installed.

Is it possible to 

Note: it seems that in `%%shell` code blocks, `cd` moves, variable
declarations, and function declarations are limited to the current
block (no effect on downstream code blocks).

Let's create some folders:

``` code
%%shell

mkdir -p src data results
```


## install dependencies

You will need to install
[vsearch](https://github.com/torognes/vsearch),
[cutadapt](https://github.com/marcelm/cutadapt/), and
[swarm](https://github.com/torognes/swarm).

### install cutadapt

using conda

...

``` code
%%shell

apt search cutadapt
```


### install swarm

We could install `swarm` and `vsearch` using `conda`, but for
educational purposes, let's compile them ourselves. We will put their
source code in the `src` folder:

``` code
%%shell

(cd ./src/
 git clone https://github.com/torognes/swarm.git
 cd ./swarm/
 make
)
```

``` code
%%shell

./src/swarm/bin/swarm --version
```

### install vsearch

``` code
%%shell

(cd ./src/
 git clone https://github.com/torognes/vsearch.git
 cd ./vsearch/
 ./autogen.sh
 ./configure CFLAGS="-O3" CXXFLAGS="-O3"
 make
)
```

``` code
%%shell

./src/vsearch/bin/vsearch --version
```

Installing [lulu](https://github.com/tobiasgf/lulu) or
[mumu](https://github.com/frederic-mahe/mumu) is not necessary.


## dataset

A subset of the Neotropical Forest Soil dataset
([PRJNA317860](https://www.ebi.ac.uk/ena/browser/view/PRJNA317860);
[Mahé et al., 2017](https://doi.org/10.1038/s41559-017-0091)),
corresponding to the following run accessions:

```
SRR23272700
SRR23272716
SRR23272737
SRR23272741
SRR23272752
SRR23272767
SRR23272778
SRR23272788
SRR23272799
SRR23272803
SRR23272822
SRR23272833
SRR23272848
SRR23272859
SRR23272860
SRR23272861
SRR23272874
SRR23272881
SRR23272890
SRR23272901
```

and subsampled at 1%, using `vsearch`:

```bash

function subsample() {
    local -ri SEED=1
    local -r PERCENTAGE="1.0"
    local -r SUBSAMPLED_FASTQ="$(sed 's/NG-7070_// ; s/_lib.*_1976//' <<< ${1})"

    vsearch \
        --fastx_subsample "${1}" \
        --randseed "${SEED}" \
        --sample_pct "${PERCENTAGE}" \
        --quiet \
        --fastqout - | \
        gzip - > "${SUBSAMPLED_FASTQ}"
}


export -f subsample

find . -name "NG-7070_*.fastq.gz" -type f -exec bash -c 'subsample "$0"' {} \;
```

Note: in a pair of R1 and R2 fastq files, both files have the same
number of reads, so using a fix seed (not zero) guarantees that
subsamplings results for both R1 and R2 fastq files are identical
(same number of reads, same reads, in the same order)


## aim

Mention some of vsearch\'s lesser known features such as **sff to
fastq** conversion. I\'d also like to demonstrate how our tools can be
piped (`|`) together to create seamless pipelines. I also had in mind
the creation of a reference database from scratch for taxonomic
assignment.


## checks
==================================

``` code
%%shell
for FASTA in ./data/*.fas ; do
    cat ${FASTA}
done
```

``` code
%%shell
#!/usr/bin/env bash

# check software versions
for SOFT in bash R python3 swarm vsearch cutadapt ; do
    "${SOFT}" --version 2>&1 | head -n 1
    echo
done
```
