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

### Set area list colors
To change colors of area lists in a TrakEM2 file.
```
python set_class_colors.py -t ../data/trakem2_n2u/n2u_vol.xml -n ../mat/neuron_category_pub.txt -c ../mat/color_code.txt -o ../data/trakem2_n2u/n2u_brittin.xml
```
where the neuron category files has the format
```
cell_name,cell_class
```
and the color code file has the format
```
cell_class,html_color
```
**Note that the output xml file will not have the proper header for TrakEM2 to read the file directly. In linux you can cat the header as follows
```
cat /mat/header.txt output.xml > trakem2_readable.xml
```
where output.xml is the output file produced by this script and
trakem2_readable is the file read by trakem2.

## Author

* **Christopher Brittin** - *Initial work* - [cabrittin](https://github.com/cabrittin)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


