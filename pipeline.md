# Pipeline

## Google Colab

Let's explore the environment:

``` code
%%shell

git --version
gcc --version
uname -a
```

Runtime seems to be based on Ubuntu LTS 20.04. The `gcc` version is old
(9.2) and does not allow to compile `mumu` locally.


## Requirements and dependencies

You will need to install
[vsearch](https://github.com/torognes/vsearch),
[cutadapt](https://github.com/marcelm/cutadapt/),
[swarm](https://github.com/torognes/swarm), and
[lulu](https://github.com/tobiasgf/lulu).

I will assume that you also have [python](https://www.python.org/)
(version 3.5 or more), [R](https://cran.r-project.org/) (version 3.5
or more), and [bash](https://www.gnu.org/software/bash/) (version 4 or
more).


``` {.bash}
%%shell

mkdir -p src data results

cd ./src/
git clone https://github.com/torognes/swarm.git
cd ./swarm/
make
```


``` {.bash}
%%shell

./src/swarm/bin/swarm --version
```


``` {.bash}
%%shell

cd ./src/

git clone https://github.com/torognes/vsearch.git
cd ./vsearch/
./autogen.sh
./configure CFLAGS="-O3" CXXFLAGS="-O3"
make
```

``` {.bash}
%%shell

./src/vsearch/bin/vsearch --version
```


## aim

Mention some of vsearch\'s lesser known features such as **sff to
fastq** conversion. I\'d also like to demonstrate how our tools can be
piped (`|`) together to create seamless pipelines. I also had in mind
the creation of a reference database from scratch for taxonomic
assignment.


### checks
==================================

``` {.bash}
%%shell
for FASTA in ./data/*.fas ; do
    cat ${FASTA}
done
```

``` {.bash}
%%shell
#!/usr/bin/env bash

# check software versions
for SOFT in bash R python3 swarm vsearch cutadapt ; do
    "${SOFT}" --version 2>&1 | head -n 1
    echo
done
```
