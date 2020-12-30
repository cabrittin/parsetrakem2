"""
measure_adjacency.py

Conversts xml adjacency file to a neighborhood file.

In this formtat, each segment is listed with each of its neighbors

created: Christopher Brittin
date: 21 October 2018

"""

import os
import argparse
from lxml import etree
import csv
from tqdm import tqdm

class Segment:
    def __init__(self,cell):
        self.cell = cell
        self.neighbors = {}

    def add_neighbor(self,layer,idx,neighbor,adj):
        if layer not in self.neighbors: self.neighbors[layer] = {}
        if idx not in self.neighbors[layer]: self.neighbors[layer][idx] = {}
        if neighbor not in self.neighbors[layer][idx]: self.neighbors[layer][idx][neighbor] = 0
        self.neighbors[layer][idx][neighbor] += adj

    def write_xml(self,root):
        for l in self.neighbors:
            layer = etree.SubElement(root,'layer')
            layer.set('name',l)
            for i in self.neighbors[l]:
               idx = etree.SubElement(layer,'idx')
               idx.set('value',str(i))
               for (k,v) in self.neighbors[l][i].items():
                   neigh = etree.SubElement(idx,'neighbor')
                   neigh.set('name',k)
                   neigh.set('adjacency',str(v))

    
    

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
    N = {}
    
    for _l in tqdm(layers,desc="Loading xml data"):
        l = root.find("layer[@name='%s']" %_l)
        areas = l.findall('area')
        for a in areas:
            c1 = a.find('cell1').text
            c2 = a.find('cell2').text
            i1 = int(a.find('index1').text)
            i2 = int(a.find('index2').text)
            adj = int(a.find('adjacency').text)
            if c1 not in N: N[c1] = Segment(c1)
            if c2 not in N: N[c2] = Segment(c2)
            N[c1].add_neighbor(_l,i1,c2,adj)
            N[c2].add_neighbor(_l,i2,c1,adj)

    #Set up xml if file if it does not exist
    data = etree.Element('data')
    xml_out = etree.tostring(data,pretty_print=False)
    with open(params.fout,'wb') as fout:
        fout.write(xml_out)
            
    #Open xml file
    tree = etree.parse(params.fout)
    root = tree.getroot()

    for (cell,v) in tqdm(N.items(),desc="Writing neighbors to xml"):
        c = etree.SubElement(root,'cell')
        c.set('name',cell)
        v.write_xml(c)

    xml_out = etree.tostring(tree,pretty_print=False)
    with open(params.fout,'wb') as fout:
        fout.write(xml_out)
