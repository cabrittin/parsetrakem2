"""
measure_adjacency.py

Measures adjacency of TrakEM2 files. For each layers, extracts all of the 
adjacency lists. Fills in any gaps in the boundary path to make the boundary 
continuous (TrakEM2 only includes a subset of points on the boundary). Determines 
the bounding box for each area list. Identifes which bounding boxes overlap. 
For each pair of overlapping bounding boxies determines if the associated boundaries 
are adjacent. Adjacency is determined by calculating pairwise distances between all 
the points in boundary 1 with all the points in boundary 2. The distances are stored 
in a matrix A of dimension [i,j] where i and j are the number of points in boundary 1 
and 2, respectively. The matrix A is then binarized such that matrix element B[i,j] = 1 
if  A[i,j] <= pixel_radius and B[i,j] = 0 otherwise, where pixel_radius is the max 
distance two pixels can be separated and still be considered adjacent (default = 10). 
The length of adjacency is determined by counting the number of rows (r) with at at 
least one column value equal to 1 and the number of columns (c) with at least one row 
value equal to 1. The adjacency length between boundary 1 and boundary 2 is defined 
as min(r,c). 

Adjacency data is written to an xml file, which can be easily updated. It's possible 
to look at only specific layers rather than all of the layers. In this case, only the 
specified layers will be updated. 

To speed up processing, layers can be processed over multiple CPUs.    

Code has only been tested on Linux OS. If running on Windows there may be formatting
issues with reading/writing to files.  

Brittin, Cook, Hall, Cohen, Emmons. 'Volumetric reconstruction of 
Caenorhabditis elegans nerve ring supports combinatorial CAM expression 
model for synaptic specificity'. (2018) Under review. 

created: Christopher Brittin
date: 17 October 2018

Required 3rd party packages:
   lxml
   argparse
   multiprocessing
   numpy
   scipy

Synposis:
   python measure_adjaceny.py trakem2 fout [OPTIONS]

Parameters:
    trakem2 (str):  The file location of the trakem2 file
    fout (str): The file location of the xml file to which data will be written
    -p, --pixel_radius (int): Pixel radius to classify adjacent boundary points
                (default is 10)
    -t, --area_thresh (int): Arear lists smaller than the area thresh are excluded
                from processing (default is 200 px^2)
    -s, --scale_bounding_box (float): Scales the bounding box. Set to greater than 1
                to ensure that all adjacent boundaries are identified in the preprocessing
                step of looking for overlapping boundary boxes. (default is 1.1)
    -n, --nprox (int): Number of CPU(s) used to process each layer. (default is 1)
    -l, --layers (str): Specify which layers to process. Separate multiple layers
                 by a ','. Make sure to use the layer names in the trakem2 file. 
                 If not specified, then all layers will be processed. 


Examples:
  General use:
     python measure_adjacency.py /path/to/trakem2 /path/to/xml

  Increase number of CPUs
     python measure_adjacency.py /path/to/trakem2 /path/to/xml -n 2

  Specify layers to be processed
     python measure_adjacency.py /path/to/trakem2 /path/to/xml -l LAYER1,LAYER2,LAYER3
   

"""
import sys
import os
import argparse
import multiprocessing as mp
import time
from lxml import etree

from parsetrakem2.parsetrakem2 import ParseTrakEM2

def time_string(_seconds):
    day = _seconds // (24 * 3600)
    _seconds = _seconds % (24 * 3600)
    hour = _seconds // 3600
    _seconds %= 3600
    minutes = _seconds // 60
    _seconds %= 60
    seconds = _seconds
    return "%d:%d:%d:%d" % (day, hour, minutes, seconds)

def submit_batch(P,o):
    adj = P.batch_compute_adjacency(o)
    return adj
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('trakem2',
                        action="store",
                        help="TrakEM2 file")

    parser.add_argument('fout',
                        action = 'store',
                        help = "Output file")

    parser.add_argument('-p','--pixel_radius',
                        dest = 'pixel_radius',
                        action="store",
                        required = False,
                        default = 10,
                        type = int,
                        help = ("Boundaries separated by distances less than or "
                                "equal to the pixel radius are classified as "
                                "adjacent. DEFAULT = 10."))
    
    parser.add_argument('-t','--area_thresh',
                        dest = 'area_thresh',
                        action = 'store',
                        required = False,
                        default = 200,
                        type = int,
                        help = ("Area lists less than area_thresh are not "
                                "considered in the adajancency analysis. "
                                "DEFAULT = 200. "))

    parser.add_argument('-s','--scale_bounding_box',
                        dest = 'scale_bounding_box',
                        action = 'store',
                        required = False,
                        default = 1.1,
                        type = float,
                        help = ("Adjusts the search radius by scaling the "
                                "area list bounding boxes. DEFAULT = 1.1. "))
    
    parser.add_argument('-n','--nproc',
                        dest = 'nproc',
                        action = 'store',
                        required = False,
                        default = 1,
                        type = int,
                        help = ("Number of jobs if running "
                                "in multiprocessor mode. DEFAULT = 1.")
                        )
    
    parser.add_argument('-l','--layers',
                        dest = 'layers',
                        action = 'store',
                        required = False,
                        default = None,
                        help = ("Specifiy which layers to analyze. "
                                "Separate layers names by ',' e.g. LAYER1,LAYER2,.. "
                                "Must use layer name specified in "
                                "//t2_patch/@title in TrakEM2 file."))

    
    params = parser.parse_args()

    print('TrakEM2 file: %s' %params.trakem2)
    print('Writing to file: %s' %params.fout)
    print('Running %d jobs' %params.nproc) 
    print('Loading TrakEM2 file...')
    P = ParseTrakEM2(params.trakem2)
    P.get_layers()
    print('Extracted %d layers.' %(len(P.layers)))
    P.get_area_lists()
    print('Extracted %d area lists.' %(len(P.area_lists)))    
    if params.layers:
        print('Analyzing layers: %s' %params.layers)
        layers = params.layers.split(',')
    else:
        print('Analyzing all layers.')
        layers = sorted(P.layers.keys())
    

    #Set up xml if file if it does not exist
    if not os.path.isfile(params.fout):
        data = etree.Element('data')
        xml_out = etree.tostring(data,pretty_print=False)
        with open(params.fout,'wb') as fout:
            fout.write(xml_out)
            
    #Open xml file
    tree = etree.parse(params.fout)
    root = tree.getroot()

    #Add layers not previously analyzed
    curr_layers = [l.get('name') for l in root.findall('layer')]
    for l in layers:
        if l in curr_layers:
            xlayer = root.find("layer[@name='%s']" %l)
            root.remove(xlayer)
        _l = etree.SubElement(root,'layer')
        _l.set('name',l)
        root.append(_l) 
    xml_out = etree.tostring(tree,pretty_print=False)
    with open(params.fout,'wb') as fout:
            fout.write(xml_out)

    

    print('Processing layers...')
    N = len(layers)
    idx = 0
    _end = '\r'
    time0 = time.time()
    for l in layers:
        time1 = time.time()
        B = P.get_boundaries_in_layer(l,area_thresh = params.area_thresh,
                                      scale_bounding_box = params.scale_bounding_box)
        overlap = P.get_overlapping_boundaries(B)        

        if params.nproc == 1:
            adj = P.batch_compute_adjacency(overlap)
        else:
            overlap_split = [overlap[i::params.nproc] for i in range(params.nproc)]
            pool = mp.Pool(processes = params.nproc)
            results = [pool.apply_async(submit_batch,args=(P,o,))
                       for o in overlap_split]
            adj = [o for p in results for o in p.get()]

          
        xlayer = root.find("layer[@name='%s']" %l)
        for (b1,b2,_adj) in adj:
            xarea = etree.SubElement(xlayer,'area')
            cell1 = etree.SubElement(xarea,'cell1')
            cell1.text = b1.name
            cell2 = etree.SubElement(xarea,'cell2')
            cell2.text = b2.name
            idx1 = etree.SubElement(xarea,'index1')
            idx1.text = str(b1.index)
            idx2 = etree.SubElement(xarea,'index2')
            idx2.text = str(b2.index)     
            xadj = etree.SubElement(xarea,'adjacency')
            xadj.text = str(_adj)             
            
        idx += 1
        if idx == N: _end = '\n'
        proc_time = time_string(time.time() - time0)
        print("Processed %d/%d layers. Last layer processed: %s. "
              "Found %d adjacencies. " 
              "Time to process last layer: %2.3f sec. "
              "Total processing time: %s. "
              %(idx,N,l,len(adj),time.time() - time1,proc_time),end=_end)
    
    xml_out = etree.tostring(tree,pretty_print=False)
    with open(params.fout,'wb') as fout:
            fout.write(xml_out)   
    print('Finished!')

    
