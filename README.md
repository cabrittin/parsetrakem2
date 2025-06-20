# parsetrakem2
Python code for parsing [TrakEM2](https://imagej.net/TrakEM2) files. 

This is the accompanying code for

Brittin, C. A. , Cook, S. J., Hall, D.H., Emmons, S. W., Cohen. N. A multiscale brain map derived from whole-brain volumetric reconstructions. Nature (2021). [[paper](https://dx.doi.org/10.1038/s41586-021-03284-x)] [[preprint](https://doi.org/10.1101/2020.05.24.112870)]

Source TrakEM2 data: [zenodo](https://zenodo.org/record/4383277#.X-wK5tZOk-I) or [wormwiring.org](http://wormwiring.org/) data.

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

### 20250531 UPDATE:
Use of the aux module is deprecated. Now using pycsvparser. In local environment pip install with

```
pip3 install pycsvparser @ git+https://github.com/cabrittin/pycsvparser.git@7a809bb74bd0313d7cd09b7a52e3eca3f9c1d926
```

## Usage

The following outline the included scripts and provide instructions of basic usage. 
For more detailed documentation, run scripts with -h flag or consult the source code documentation.
Main scripts are located in trakem2/.  Test scripts are located in test/.

### Measure adjacency
To measure adjacency of segmented TrakEM2 file use measure_adjacency.py:
```
python scripts/measure_adjacency.py /path/to/trakem2.xml /path/to/output.xml
```
Adjacency data is written to an xml file which can be easily updated. Script can be run in parallel over muliple CPU(s).

### Convert xml output to csv
For convenience, the xml2csv.py will convert the output xml from measure_adjacency.py to csv format:
```
python scripts/xml2csv.py /path/to/xml /path/to/csv
```
Format of the csv is as follows:
```
cell_1,cell_2,index_1,index_2,layer_name,adjacency_length
```

### Set area list colors
To change colors of area lists in a TrakEM2 file.
```
python scripts/set_class_colors.py -t ../data/trakem2_n2u/n2u_vol.xml -n ../mat/neuron_category_pub.txt -c ../mat/color_code.txt -o ../data/trakem2_n2u/n2u_vol_modified.xml
```
where the neuron category files has the format
```
cell_name,cell_class
```
and the color code file has the format
```
cell_class,html_color,color_intensity
```
**Note that the output xml file will not have the proper header for TrakEM2 to read the file directly. In linux you can cat the header as follows
```
cat /mat/header.txt output.xml > trakem2_readable.xml
```
where output.xml is the output file produced by this script and
trakem2_readable is the file read by trakem2.

### Extract centroid and area of each neurite segment 
```
python scripts/extract_segmentation_stats.py trakem2_xml_file output file
```

### Extract subvolumes from the volumetric data 
```
python scripts/modify_rendering.py mat/config_modify_rendering_example.ini
```
The example config file mat/config_modify_rendering_example.ini will generate Fig 1a from Brittin et al.


## Author

* **Christopher Brittin** - *Initial work* - [cabrittin](https://github.com/cabrittin)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


