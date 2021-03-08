# FAStools

FAStools is a compendium of perl scripts that I've written durign the years to manage and analyze FASTA-formatted biological sequences. This repository will be populated as long as the scripts are ported to a release verison (e.g. used in a publication), since by now they have a "my personal use" design.

The scripts are intended to be used in the linux/macOSX shell. Although some of them script may operate well in Windows, they are neither written nor intended nor maintained for MS. And they will never be.

FAStools are written in perl, but does not use perl modules. Sometimes it makes extensive use of basic linux shell commands.

As a distinctive trait, FAStools script names will always begin with the capital FAS prefix. This is useful for localizing the scripts if users store them in PATH.

Each script has its own mini-help available if run with no options. Sequences can generally be in a file (gzip is supported via zcat) or they con be piped-in.
