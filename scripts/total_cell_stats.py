"""
total_cell_stats.py

Computes the total volume and surface area of cells.
Takes as input the output of extract_segmentation_stats.py

created: Christopher Brittin
date: 21 October 2018

"""

import argparse
from lxml import etree
import csv

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
 
    parser.add_argument('xml',
                        action="store",
                        help="XML file")

    parser.add_argument('fout',
                        action = 'store',
                        help = "Output csv file")
    
    parser.add_argument('--scale_adult',
                        dest='scale_adult',
                        action = 'store_true',
                        required=False,
                        default=True,
                        help='Apply scaling for the adult')

    params = parser.parse_args()

    tree = etree.parse(params.xml)
    root = tree.getroot()
    
    data = {}
    for l in root.findall('layer'):
        segs = l.findall('segment')
        _scale = 1
        if params.scale_adult:
            lname = l.get('name')
            if 'vc' in lname or 'VC' in lname:
                _scale = 2
        
        for s in segs:
            name = s.find('name').text
            area = int(s.find('area').text)
            length = int(s.find('length').text)
            if name not in data: data[name] = {'volume':0,'sa':0}
            data[name]['volume'] += area
            data[name]['sa'] += length*_scale

    with open(params.fout,'w') as fout:
        for i in sorted(data.keys()):
            tmp = ','.join([i,str(data[i]['sa']),str(data[i]['volume'])])
            fout.write(tmp + '\n')
