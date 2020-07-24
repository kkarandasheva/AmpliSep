# AmpliSep
AmpliSep is a tool for dividing amplicon-based sequencing reads by pool based on their alignment coordinates

## General notes
split .bam by pools

## User Guide
### Installation and Requirements
You need to have your data aligned (.bam files).


### To use
**Options**
| Options  | Input  | Description |
| :----------------------- |:------:|:---------------|
| <pre lang="console">-f, --file</pre> | FILE | Specifies the file name of the reads to divided. The input reads<br> can be only in BAM format. |
| <pre lang="console">-o, --outdir</pre> | DIR | Specifies the output directory name. If this is not specified, the<br> output will be written to the same folder where the input file is. |
| <pre lang="console">-d, --design</pre> | FILE | Specifies the target regions file name. The input should be in<br> BED format. |

**Command line example**
```console
foo@bar:~$ python amplisep.py -f test.bam -o /home/user/result/ -d design.bed
```
