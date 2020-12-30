"""
measure_adjacency.py

Conversts xml file to csv file.

created: Christopher Brittin
date: 21 October 2018

"""

import argparse
from lxml import etree
import csv

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
                                     ("Measure adjaceny in segmented" 
                                     "TrakEM2 file"))
    parser.add_argument('xml',
                        action="store",
                        help="XML file")

    parser.add_argument('fout',
                        action = 'store',
                        help = "Output csv file")

    params = parser.parse_args()

    tree = etree.parse(params.xml)
    root = tree.getroot()
    layers = sorted([l.get('name') for l in root.findall('layer')])
    data = []
    for _l in layers:
        l = root.find("layer[@name='%s']" %_l)
        areas = l.findall('area')
        for a in areas:
            c1 = a.find('cell1').text
            c2 = a.find('cell2').text
            i1 = a.find('index1').text
            i2 = a.find('index2').text
            adj = a.find('adjacency').text
            data.append([c1,c2,i1,i2,_l,adj])
    

    with open(params.fout, "w") as f:
        writer = csv.writer(f)
        writer.writerows(data)
    
