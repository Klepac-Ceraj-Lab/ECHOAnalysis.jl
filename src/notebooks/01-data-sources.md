# Data for Analysis

## Description

Raw data for this paper can be found here **link to OSF repo**

In order to perform the analyses outlined in the notebooks in this project,
**description of DataDeps**

## Metadata

Patient metadata is kept in a FileMaker Pro database,
and exported into a `csv` file that requires some effort to parse effectively.

### Export from FilemakerPro

The first step is to
export the relevant sections of the database as `.xlsx` (MS Excel) files.
First, make sure all records are being displayed
by clicking the `Show All` button (if it is not greyed out).

![Show all button](src/img/01.showall.png)

Next, export one section as an excel file.
The following demonstrates the steps
to get the metadata associated with timepoints,

1. Select `File > Export Records`
1. Select the location to save the file (eg `~/Desktop/timepoints`),
    and click `Save`
1. When the `Excel Options` dialogue box pops up, click `Continue`
1. In the `Specify Field Order for Export` dialogue box:
    1. If you've previously exported anything, clear those files
        by clicking `Clear All`
    1. Select a table to export.
        In this example, we'll select `Timpointinfo`
    1. Then click `Move All`, then `Export`

![Export](src/img/01.export.png)

![Save](src/img/01.save.png)

![Clear previous](src/img/01.clearall.png)

![Table select](src/img/01.tableselect.png)

![Move all and export](src/img/01.moveexport.png)

### Simple Scrub

The next step is to convert this file to a `.csv`.
This can be done using software like MS Excel, Apple Numbers, or OpenOffice Calc.

Once you have a `.csv` file,
run the `simple_scrub.jl` script contained in the `bin/` folder.
This script does a few simple things.
first, it renames the headers to be a bit easier to work with,
replacing spaces with underscores and `::` with two underscores.
Second, it removes empty columns and rows
(that is, columns and rows that have only missing values).

For more information, run the script with the `--help`:

```
$ julia --project=@. bin/simple_scrub.jl --help
usage: simple_scrub.jl [-d] [-v] [-q] [-l LOG] [--delim DELIM]
                       [--dry-run] [-o OUTPUT] [-h] input

positional arguments:
  input                Table to be scrubbed. By default, this will be
                       overwritten

optional arguments:
  -d, --debug          Show @debug level logging messages
  -v, --verbose        Show @info level logging messages
  -q, --quiet          Only show @error logging messages
  -l, --log LOG        Write logs to this file. Default - no file
  --delim DELIM         (default: ",")
  --dry-run            Show logging, but take no action. Most useful
                       with --verbose
  -o, --output OUTPUT  write to this file instead of overwriting
                       original
  -h, --help           show this help message and exit
```

For example, to run this script on our example above:

```sh
$ julia --project=@. bin/simple_scrub.jl ~/Desktop/timepoints.csv \
  -o data/metadata/timepoints.csv \
  -vl data/metadata/timepoints.log
```

Additional transformations of the metadata to get it into a usable form
can be found in the next notebook.

## Metagenome data

Compressed quality-scored sequencing files (`.fastq.gz`)
from the sequencing facility were concatenated
and run through the bioBakery metagenome workflow.
Each sample was sequenced on 4 lanes of Illumina {{TODO: add info about Illumina}}
generating paired-end reads.

### Concatenating raw sequencing files

Each set of 4 forward and reverse reads
were concatenated using the `catfastq.jl` script.

```
$ julia --project=@. bin/catfastq.jl --help
usage: catfastq.jl [-i IN-FOLDER] [--delete] [-o OUTPUT-FOLDER]
                   [--dry-run] [--debug] [-v] [-q] [-l LOG] [-h]

optional arguments:
  -i, --in-folder IN-FOLDER
                        Folder containing fastq files (default: ".")
  --delete              set to remove original files after
                        concatenation
  -o, --output-folder OUTPUT-FOLDER
                        destination folder for concatenated files
                        (default: "./")
  --dry-run             print messages but do nothing
  --debug               Show @debug level logging messages
  -v, --verbose         Show @info level logging messages
  -q, --quiet           Only show @error logging messages
  -l, --log LOG         Write logs to this file. Default - no file
  -h, --help            show this help message and exit
```

To run, assuming you have a folder containing the `fastq.gz` files
from the sequencing facility in a folder called `rawfastq/`:

```
$ julia --project=@. bin/catfastq.jl -i rawfastq/ -o catfastq/ -vl batch001_concatenate.log
```

### Running the bioBakery workflow

The bioBakery `wmgx` workflow was used on all samples
with the following software versions:

- Anadama2 v0.5
- bioBakery Workflows v0.10
- kneaddata v0.7
- metaphlan2 v2.2
- humann2 v0.11.1

The following command was run on the `regal` compute cluster from Harvard
research computing (a Centos7 environment).

```
$ biobakery_workflows wmgx --input catfastq/ --output analysis/ \
  --bypass-strain-profiling --pair-identifier .1 \
  --grid-jobs 96 --grid-partition serial_requeue,shared,240 # grid options
```

## 16S data

<!-- TODO: Add info about 16S analysis -->
