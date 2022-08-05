# NucProcess

NucProcess is a Python program to process DNA read-pair data obtained from Hi-C experiments, with a primary emphasis on handling datasets derived from single-cells. Processing single-cell Hi-C contacts, which relate to only one underlying 3D geometry, enables the calculation of whole genome structures, e.g using [NucDynamics](https://github.com/tjs23/nuc_dynamics). NucProcess may be used for canonical, multi-cell Hi-C (see -p option), but by default any input data is assumed to represent only one cell.

This software takes paired FASTQ sequence read files, a reference genome sequence and knowledge of experimental parameters (restriction enzyme type and the range of size of DNA fragments sequenced in the library) to create processed, Hi-C contact files. By default, the output is generated in Nuc Chromatin Contact (NCC) format, which is a simple text based format described in the README and which is directly compatible with [NucDynamics](https://github.com/tjs23/nuc_dynamics). A processing report document and a contact map will also be generated for each run; in SVG format that can be viewed in most web browsers.

## Versions

There are two main versions of NucProcess, a stable release that corresponds to the primary publications and a development version where new features are being added and tested. We recommend using the stable release under normal circumstances.

The current stable version is available here: [Release 1.3](https://github.com/tjs23/nuc_processing/tree/release_1.3)

The development version is available here: [Development master](https://github.com/tjs23/nuc_processing/tree/master)  

## Documentation

See the [wiki pages](https://github.com/tjs23/nuc_processing/wiki) for documentation on how to install and run NucProcess,
including options for the following commands:

* [nuc_process](https://github.com/tjs23/nuc_processing/wiki/nuc_process)
* [nuc_contact_map](https://github.com/tjs23/nuc_processing/wiki/nuc_contact_map) 
* [nuc_contact_probability](https://github.com/tjs23/nuc_processing/wiki/nuc_contact_map)
* [split_fastq_barcodes](https://github.com/tjs23/nuc_processing/wiki/split_fastq_barcodes)

## Citations

If you use NucProcess in published work, please cite the following reference:

  Stevens et al. Nature. 2017 Apr 6;544(7648):59-64 [PMID:28289288](https://www.ncbi.nlm.nih.gov/pubmed/28289288)

## Licensing

NucProcess is licensed under [GNU Lesser General Public License v3.0](https://github.com/tjs23/nuc_processing/blob/release_1.0/COPYING.LESSER)
