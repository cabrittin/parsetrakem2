"""
xml2csv_area.py

Conversts xml output from extract_segmentation_stats.py to csv 

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
    data = [['layer_name','object_name','object_index','center_x','center_y','area','boundary_length']]
    for _l in layers:
        l = root.find("layer[@name='%s']" %_l)
        segments = l.findall('segment')
        for s in segments :
            name = s.find('name').text
            index = s.find('index').text
            centx = s.find('centx').text
            centy = s.find('centy').text
            area = s.find('area').text
            length = s.find('length').text 
            data.append([_l,name,index,centx,centy,area,length])
    

    with open(params.fout, "w") as f:
        writer = csv.writer(f)
        writer.writerows(data)
