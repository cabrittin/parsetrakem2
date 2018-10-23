# parsetrakem2
Python code for parsing [TrakEM2](https://imagej.net/TrakEM2) files. 
The code is intended to link TrakEM2 data with [wormwiring.org](http://wormwiring.org/) data.

## Installation
Clone or download respository. Make sure to maintain the relative paths.

### Prerequisites
Python3 will need to be installed along with the following python packages 
```
Numpy v1.15.2
Scipy v1.1.0
Matplotlib v3.0.0
lxml v4.2.5
multiprocessing-on-dill v3.5.0a4
```
**Note that multiprocessing appears to work on Python 3.5 but not on 3.6 and 3.7. So keep this in mind when choosing a python version. 
## Usage

The following outline the included scripts and provide instructions of basic usage. 
For more detailed documentation, run scripts with -h flag or consult the source code documentation.
Main scripts are located in trakem2/.  Test scripts are located in test/.

### Measure adjacency
To measure adjacency of segmented TrakEM2 file use measure_adjacency.py:
```
python measure_adjacency.py /path/to/trakem2.xml /path/to/output.xml
```
Adjacency data is written to an xml file which can be easily updated. Script can be run in parallel over muliple CPU(s).

### Convert xml output to csv
For convenience, the xml2csv.py will convert the output xml from measure_adjacency.py to csv format:
```
python xml2csv.py /path/to/xml /path/to/csv
```
Format of the csv is as follows:
```
cell_1,cell_2,index_1,index_2,layer_name,adjacency_length
```

## Author

* **Christopher Brittin** - *Initial work* - [cabrittin](https://github.com/cabrittin)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


