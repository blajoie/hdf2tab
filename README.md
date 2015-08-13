<img height=40 src='http://my5C.umassmed.edu/images/3DG.png' title='3D-Genome' />
&nbsp;&nbsp;
<img height=30 src='http://my5C.umassmed.edu/images/dekkerlabbioinformatics.gif' />
&nbsp;&nbsp;
<img height=30 src='http://my5C.umassmed.edu/images/umasslogo.gif' />

# hdf2tab

h5 helper scripts for encode/dekkerlab hi-c data.

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

See the [hdf2tab Wiki](https://github.com/blajoie/hdf2tab/wiki) for full documentation of encode/dekkerlab Hi-C HDF5 file format.
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

hdf2tab can read/subset a hdf5 file containing Hi-C data from the encode/dekkerlab and convert it into a my5C fomatted tsv matrix file.

## Usage

```

$ python hdf2tab/hdf2tab.py  --help

usage: hdf2tab.py [-h] -i INFILE [-v] [-o OUTFILE]
                  [-z ZOOM_COORDS [ZOOM_COORDS ...]] [-e ELEMENTSFILE]
                  [-m BLOCKMEM] [-p PRECISION] [--info] [--or] [--cis]
                  [--chrs SELECTED_CHRS [SELECTED_CHRS ...]]
                  [--elementexten ELEMENTEXTEN] [--maxdim MAX_DIMENSION]
                  [--outputchrs] [--outputbins] [--outputfactors] [--version]

Extract c-data from HDF5 file into TXT (matrix.gz)

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --input INFILE
                        interaction matrix hdf5 file (default: None)
  -v, --verbose         Increase verbosity (specify multiple times for more)
                        (default: None)
  -o OUTFILE, --output OUTFILE
                        interaction matrix output file (default: None)
  -z ZOOM_COORDS [ZOOM_COORDS ...], --zoom ZOOM_COORDS [ZOOM_COORDS ...]
                        zoom coordinate (can only select symmetrical subsets)
                        (default: None)
  -e ELEMENTSFILE, --elements ELEMENTSFILE
                        elements bed file (bed3+ format) (default: None)
  -m BLOCKMEM, --bmem BLOCKMEM
                        block size for extracting (default=hdf chunk size)
                        (default: None)
  -p PRECISION          output precision (# of digits) (default: 4)
  --info                interaction matrix hdf5 file (default: False)
  --or                  output file relative to input file path (default:
                        False)
  --cis                 extract cis maps only (default: False)
  --chrs SELECTED_CHRS [SELECTED_CHRS ...]
                        subset of chromosomes to extract, [+] = all, [-] =
                        none, zoom selected overrides --chrs (default:
                        ['default'])
  --elementexten ELEMENTEXTEN
                        bp extension for all elements, +/- bp to start/end of
                        each element (default: 0)
  --maxdim MAX_DIMENSION
                        maximum dimension of allxall matrix - else cis only
                        (default: 4000)
  --outputchrs          output the chromosome list file, no matrix output
                        (default: False)
  --outputbins          output the bin position file, no matrix output
                        (default: False)
  --outputfactors       output the balancing factor list file, no matrix
                        output (default: False)
  --version             show program's version number and exit

## Usage Examples

(dump info of HDF5 file)
$ python hdf2tab/hdf2tab.py -i test/C-10000000/iced/NPC__mm9__genome__C-10000000-iced.hdf5 --info

(dump balance factors from HDF5 file)
$ python hdf2tab/hdf2tab.py -i test/C-10000000/iced/NPC__mm9__genome__C-10000000-iced.hdf5 --outputfactors
	-> NPC__mm9__genome__C-10000000-iced.factors

(dump genome wide tsc matrix
$ python hdf2tab/hdf2tab.py -i test/C-10000000/iced/NPC__mm9__genome__C-10000000-iced.hdf5 -v

(dump all cis maps)
$ python hdf2tab/hdf2tab.py -i test/C-10000000/iced/NPC__mm9__genome__C-10000000-iced.hdf5 --cis -v

(dump cis maps for chr1, chr2 and chr3 from the HDF5 file)
$ python hdf2tab/hdf2tab.py -i test/C-10000000/iced/NPC__mm9__genome__C-10000000-iced.hdf5 --cis --chrs chr1 chr2 chr3 -v 
	-> NPC__mm9__genome__C-10000000-iced__chr1.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr2.matrix.gz
	-> NPC__mm9__genome__C-10000000-iced__chr3.matrix.gz
	
(subset HDF by bed file)	
$ python hdf2tab/hdf2tab.py -i test/C-10000000/iced/NPC__mm9__genome__C-10000000-iced.hdf5 -e test/chromosome-bands-qB.bed -v
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

