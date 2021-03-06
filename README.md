# AmpliSep
AmpliSep is a tool for dividing amplicon-based sequencing reads by pool based on their alignment coordinates

## General notes
split .bam by pools

## User Guide
*The following information is relevant for AmpliSep v.1.0.1*
### Synopsis
amplisep.py [options]

### Installation and Requirements
* You need to have your data aligned (.bam files).
* AmpliSep uses ***pysam*** python module, so it must be preinstalled.

### To use
| Options  | Input  | Description |
| :----------------------- |:------:|:---------------|
| <pre lang="console">-f, --file</pre> | FILE | Specifies the file name of the reads to divided. The input reads<br> can be only in BAM format. |
| <pre lang="console">-o, --outdir</pre> | DIR | Specifies the output directory name. If this is not specified, the<br> output will be written to the same folder where the input file is. |
| <pre lang="console">-d, --design</pre> | FILE | Specifies the target regions file name. The input should be in<br> BED format. |

**Command line example**
```console
foo@bar:~$ python amplisep.py -f test.bam -o /home/user/result/ -d design.bed
```
## License
AmpliSep is free for non-commercial use without warranty.
