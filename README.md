# Mobile network base station scanner on Raspberry Pi

This repository the cod for my bachelor thesis "Mobile network base station scanner on Raspberry Pi" done at Faculty of Electrical Engineering and Communication at the Brno University of Technology.

## Abstract

This thesis aims to create a simple LTE base station scanner that uses Raspberry Pi as a main platform connected to a SDR. The goal is to create an application which will output the physical cell identity of a nearby cell using the cell's synchronization signals.

## Dependencies

To use this application you are going to need these python libraries

``` python
pyadi-iio
numpy
scipy
matplotlib
```

## Usage

To see the full list of commands, please either refer to the thesis or run the shell script without any parameters. Make sure to change the file permissions to make the shell script executable

``` shell
# Example command
./lte-cellscan -f 806e6,796e6 -n 3 -d
```
