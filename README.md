How to install
==============

Make sure you have matlab installed and then download the [latest release](https://github.com/kaizhang/RFECS/releases).

How to use
==========

Unpack the tarball and type `rfecs -h` to see help information.

As an example, you can go to the `example` directory, and run `rfecs input.yaml RFECS_output -p ../RFECS_core -g hg19 --chrom chr21`. After 10~20 minutes, you will see the results in a directory named `RFECS_output`.

Build from source
=================

1. Download the latest release of [stack tool](https://github.com/commercialhaskell/stack/releases).
2. `stack install`.
