# SuperCharge
A tool for merging or extending CHGCAR of VASP

## Introduction

Three python scripts were included.
- BaderMerge.py
- PoscarExtend.py
- SuperCharge.py
SuperCharge.py has covered all of the functions needed.

## Goals

To extend the potentials in CHGCAR as a S1\*S2\*S3 times. The original matrix was Na\*Nb\*Nc and final matrix will be S1Na\*S2Nb\*S3Nc.
If the atoms selected is not 'All', the atoms selected will be calculated from Bader files and output same as a S1Na\*S2Nb\*S3Nc matrix.

## Rely

To run the script, please check the related packages of Python:
- Python 3.6.5 (Recommended)
- Pandas
- numpy

The necessary files must exist with the script in the same foloder:
- CHGCAR
- parameter.txt
- BvAt00xx.dat (Bader files if atoms selected is not 'All')

## Parameter format
All of the parameters was included in the file parameter.txt. It has the following format:
- (LINE#1) S1 S2 S3
- (LINE#2) all/list/range

If second line is 'list':
- (LINE#3) atom number1
- (LINE#4) atom number2
- ...

If second line is 'range':
- (LINE#3) Axis (x/y/z)
- (LINE#4) minium number
- (LINE#5) maxium number

*Example*
```
2 3 1
range
0
0.5
```

## Output

Two files will be generated from the script:
- CHGCAR.new
- POSCAR.new

## Usage
```
python SuperCharge.py
```
