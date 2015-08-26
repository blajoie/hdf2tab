<img height=40 src='http://my5C.umassmed.edu/images/3DG.png' title='3D-Genome' />
&nbsp;&nbsp;
<img height=30 src='http://my5C.umassmed.edu/images/dekkerlabbioinformatics.gif' />
&nbsp;&nbsp;
<img height=30 src='http://my5C.umassmed.edu/images/umasslogo.gif' />

# hdf2tab

convert hdf5 (h5) files into tsv (tab seperated value) text (matrix) files.

```
hdf2tab/
	hdf2tab.py - Extract c-data from HDF5 file into TXT (matrix.gz)
```

## Installation

The main requirements are numpy and h5py.

You can install the dependencies with

for req in $(cat requirements.txt); do pip install $req; done

## Full Documentation

See the [hdf2tab Wiki](https://github.com/blajoie/hdf2tab/wiki) for full documentation of Hi-C HDF5 file format.
<br>

See the [HDF5 spec Wiki](https://github.com/blajoie/hdf2tab/wiki/H5-Spec) for focumentation of the dekkerlab Hi-C HDF5 file format. (h5 dataset/attributes)

Download/Clone the [hdf2tab.py](https://github.com/blajoie/hdf2tab) HDF5 helper script (python).
<br>
numpy/scipy/h5py required. Python 2.7+

See the [Usage Wiki](https://github.com/blajoie/hdf2tab#usage</a>) for help running the hdf2tab.py scripy.

Hi-C data is pushed through the dekkerlab standard pipeline (to be available on GIT soon) and binned at multiple bin sizes (resolutions):
- 10mb
- 2.5mb
- 1mb
- 500kb
- 250kb
- 100kb
- 40kb

See the [Wiki](https://github.com/blajoie/hdf2tab/wiki) for full documentation, examples, operational details and other information.

## Communication

- [Bryan Lajoie](https://github.com/blajoie)
- [Noam Kaplan](https://github.com/NoamKaplan)
- Twitter: [@my5C](https://twitter.com/my5C)

## What does it do?

hdf2tab can read/subset a hdf5 file containing Hi-C data into a my5C fomatted tsv matrix file.

## Usage

```

$ python scripts/hdf2tab.py  --help

usage: hdf2tab.py [-h] -i IN_FILE [-v] [--info] [-o OUT_FILE]
				  [-wm {cis,seperate,all}] [-z ZOOM_COORDS [ZOOM_COORDS ...]]
				  [-xz X_ZOOM_COORDS [X_ZOOM_COORDS ...]]
				  [-yz Y_ZOOM_COORDS [Y_ZOOM_COORDS ...]]
				  [-e BED_FILE [BED_FILE ...]]
				  [-xe X_BED_FILE [X_BED_FILE ...]]
				  [-ye Y_BED_FILE [Y_BED_FILE ...]]
				  [--element_exten ELEMENT_EXTEN]
				  [-c SELECTED_CHRS [SELECTED_CHRS ...]]
				  [-xc X_SELECTED_CHRS [X_SELECTED_CHRS ...]]
				  [-yc Y_SELECTED_CHRS [Y_SELECTED_CHRS ...]] [-b BLOCKSIZE]
				  [-p PRECISION] [--maxdim MAX_DIMENSION] [--outputchrs]
				  [--outputbins] [--outputfactors] [--version]

Extract c-data from HDF5 file into TXT (matrix.gz)

optional arguments:
  -h, --help			show this help message and exit
  -i IN_FILE, --input IN_FILE
						interaction matrix hdf5 file (default: None)
  -v, --verbose		 Increase verbosity (specify multiple times for more)
						(default: None)
  --info				interaction matrix hdf5 file (default: False)
  -o OUT_FILE, --output OUT_FILE
						interaction matrix output file (default: None)
  -wm {cis,seperate,all}
						write mode (cis=cis only maps, seperate=all cis/trans
						seperate fiels, all=single matrix) (default: all)
  -z ZOOM_COORDS [ZOOM_COORDS ...], --zoom ZOOM_COORDS [ZOOM_COORDS ...]
						x/y axis zoom coordinate (default: [])
  -xz X_ZOOM_COORDS [X_ZOOM_COORDS ...], --xzoom X_ZOOM_COORDS [X_ZOOM_COORDS ...]
						x axis zoom coordinate (default: [])
  -yz Y_ZOOM_COORDS [Y_ZOOM_COORDS ...], --yzoom Y_ZOOM_COORDS [Y_ZOOM_COORDS ...]
						y axis zoom coordinate (default: [])
  -e BED_FILE [BED_FILE ...], --elements BED_FILE [BED_FILE ...]
						x/y axis elements bed file (bed3+ format) (default:
						[])
  -xe X_BED_FILE [X_BED_FILE ...], --xelements X_BED_FILE [X_BED_FILE ...]
						x axis elements bed file (bed3+ format) (default: [])
  -ye Y_BED_FILE [Y_BED_FILE ...], --yelements Y_BED_FILE [Y_BED_FILE ...]
						y axis elements bed file (bed3+ format) (default: [])
  --element_exten ELEMENT_EXTEN
						bp extension for all elements, +/- bp to start/end of
						each element (default: 0)
  -c SELECTED_CHRS [SELECTED_CHRS ...], --chrs SELECTED_CHRS [SELECTED_CHRS ...]
						x/y axis subset of chromosomes to extract (default:
						[])
  -xc X_SELECTED_CHRS [X_SELECTED_CHRS ...], --xchrs X_SELECTED_CHRS [X_SELECTED_CHRS ...]
						x axis subset of chromosomes to extract (default: [])
  -yc Y_SELECTED_CHRS [Y_SELECTED_CHRS ...], --ychrs Y_SELECTED_CHRS [Y_SELECTED_CHRS ...]
						y axis subset of chromosomes to extract (default: [])
  -b BLOCKSIZE, --blocksize BLOCKSIZE
						block size of HDF5 file (default: None)
  -p PRECISION		  output precision (# of digits) (default: 4)
  --maxdim MAX_DIMENSION
						maximum dimension of allxall matrix - else cis only
						(default: 12000)
  --outputchrs		  output the chromosome list file, no matrix output
						(default: False)
  --outputbins		  output the bin position file, no matrix output
						(default: False)
  --outputfactors	   output the balancing factor list file, no matrix
						output (default: False)
  --version			 show program's version number and exit

```
  
## Usage Examples

```

(dump info of HDF5 file)
$ python scripts/hdf2tab.py -i test/C-10000000/iced/NPC__mm9__genome__C-10000000-iced.hdf5 --info

(dump balance factors from HDF5 file)
$ python scripts/hdf2tab.py -i test/C-10000000/iced/NPC__mm9__genome__C-10000000-iced.hdf5 --outputfactors
	-> NPC__mm9__genome__C-10000000-iced.factors

(dump genome wide tsv matrix
$ python scripts/hdf2tab.py -i test/C-10000000/iced/NPC__mm9__genome__C-10000000-iced.hdf5 -v

(dump all cis maps)
$ python scripts/hdf2tab.py -i test/C-10000000/iced/NPC__mm9__genome__C-10000000-iced.hdf5 -wm cis -v

(dump cis maps for chr1, chr2 and chr3 from the HDF5 file)
$ python scripts/hdf2tab.py -i test/C-10000000/iced/NPC__mm9__genome__C-10000000-iced.hdf5 -wm cis --chrs chr1 chr2 chr3 -v 
	-> NPC__mm9__genome__C-10000000-iced__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr2.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr3.matrix.gz

(dump chr1 (y-axis) vs genome (x-axis), single matrix)
$ python scripts/hdf2tab.py -i test/C-10000000/iced/NPC__mm9__genome__C-10000000-iced.hdf5 -wm all -yc chr1 -v
	-> NPC__mm9__genome__C-10000000-iced.matrix.gz

(dump chr1 (y-axis) vs genome (x-axis), seperate chr:chr matrices)
$ python scripts/hdf2tab.py -i test/C-10000000/iced/NPC__mm9__genome__C-10000000-iced.hdf5 -wm seperate -yc chr1 -v
	-> NPC__mm9__genome__C-10000000-iced__chr10__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr11__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr12__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr13__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr14__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr15__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr16__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr17__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr18__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr19__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr1__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr2__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr3__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr4__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr5__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr6__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr7__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr8__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr9__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chrM__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chrX__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chrY__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced.matrix.gz

(chr1,chr2,chr3,chrX:20000000-80000000,chr7:10000000-20000000 chr8:30000000-40000000 on y-axis)
(chr4,chr5,chr6,chrX:20000000-80000000,chr9:50000000-60000000 chr10:70000000-80000000 on x-axis)
(output single matrix file as 'NPC-zoomA6253'.matrix.gz file)
$ python scripts/hdf2tab.py -i test/C-10000000/iced/NPC__mm9__genome__C-10000000-iced.hdf5 -wm all -yc chr1 chr2 chr3 -xc chr4 chr5 chr6 -z chrX:20000000-80000000 -yz chr7:10000000-20000000 chr8:30000000-40000000 -xz chr9:50000000-60000000 chr10:70000000-80000000 -o NPC-zoomA6253
	-> NPC-zoomA6253.matrix.gz
	
(subset HDF by bed file)	
$ python scripts/hdf2tab.py -i test/C-10000000/iced/NPC__mm9__genome__C-10000000-iced.hdf5 -e test/chromosome-bands-qB.bed -v
	-> NPC__mm9__genome__C-10000000-iced.matrix.gz
	
```

## Change Log

## Bugs and Feedback

For bugs, questions and discussions please use the [Github Issues](https://github.com/blajoie/hdf2tab/issues).

## LICENSE

Licensed under the Apache License, Version 2.0 (the 'License');
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

<http://www.apache.org/licenses/LICENSE-2.0>

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an 'AS IS' BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

