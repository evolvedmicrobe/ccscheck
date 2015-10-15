### ccscheck - Evaluate the empirical accuracy and diagnose any issues with PacBio CCS sequencing.


```
    ccscheck input.ccs.bam out_folder_name any_size_ref.fasta
```



ccscheck will produce the following data files in the output folder specified when it completes.


| CSV File Name      | Data Contained                                                                                     |
| :-------          | :-----------                                                                                       |
| qv_calibration | Empirical counts of the correct and incorrect basecalls divided by insertions, deletions and SNPs. |
| snrs       | The SNR values for the emitted ZMWs.                                                               |
| variants   | A detailed description of all variants found when aligning the reads.                              |
| zmws       | A detailed description and associated statistics for each emitted ZMW.                             |
| zscores    | Z-scores for each ZMW/subread combination.                                                         |

This tool is an internal tool that we are making available for external users, but it is not officially supported by PacBio.

##Stand alone binaries
Are only available as betas right now but can be downloaded for:
* [Mac OSX El Capitan](http://www.evolvedmicrobe.com/ccscheck/ccscheck_macosx.tar.gz)
* [Ubuntu 14.04] (http://www.evolvedmicrobe.com/ccscheck/ccscheck.ubuntu14.tar.gz)


##Build requirements

* POSIX/64-bit/Little-Endian
* git
* C Compiler
* C# Compiler with CLR available.
* python


##Build instructions

```
    git clone https://github.com/evolvedmicrobe/ccscheck.git
    cd ccscheck
    ./build.sh
```
