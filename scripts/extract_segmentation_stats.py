"""
extract_segmentation_stats.py

Extracts the centroid and area of each segmentation.

created: Christopher Brittin
date: 01 November 2018

Required 3rd party packages:
   lxml
   argparse
   multiprocessing
   numpy
   scipy

Synposis:
   python extract_segmentation_stats.py trakem2 fout [OPTIONS]

Parameters:
    trakem2 (str):  The file location of the trakem2 file
    fout (str): The file location of the xml file to which data will be written

"""
import os
import argparse
from lxml import etree

from parsetrakem2.parse import ParseTrakEM2

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('trakem2',
                        action="store",
                        help="TrakEM2 file")

    parser.add_argument('fout',
                        action = 'store',
                        help = "Output file")

    parser.add_argument('-t','--area_thresh',
                        dest = 'area_thresh',
                        action = 'store',
                        required = False,
                        default = 200,
                        type = int,
                        help = ("Area lists less than area_thresh are not "
                                "considered in the adajancency analysis. "
                                "DEFAULT = 200. "))    
    
    parser.add_argument('--huge_tree',
                        dest='huge_tree',
                        action='store_true',
                        default=False,
                        required=False,
                        help=("Set flag if loading a large xml file, i.e. if lxml "
                            "throws a 'use XML_PARSE_HUGE option' error.")
                        )

    params = parser.parse_args()

    print('TrakEM2 file: %s' %params.trakem2)
    print('Writing to file: %s' %params.fout)
    print('Loading TrakEM2 file...')
    P = ParseTrakEM2(params.trakem2, huge_tree=params.huge_tree)
    P.get_layers()
    print('Extracted %d layers.' %(len(P.layers)))
    P.get_area_lists()
    
    #Set up xml if file if it does not exist
    if not os.path.isfile(params.fout):
        data = etree.Element('data')
        xml_out = etree.tostring(data,pretty_print=False)
        with open(params.fout,'wb') as fout:
            fout.write(xml_out)
            
    #Open xml file
    tree = etree.parse(params.fout)
    root = tree.getroot()

    print('Processing layers...')
    layers = sorted(P.layers.keys())
    for l in layers:
        print('Processed layer: %s' %l)
        B = P.get_boundaries_in_layer(l,area_thresh = params.area_thresh)
        xlayer = etree.SubElement(root,'layer')
        xlayer.set('name',l)
        root.append(xlayer)
        for _name in B.keys():
            for idx in B[_name]:
                B[_name][idx].set_centroid()
                B[_name][idx].set_boundary_length()
                xseg = etree.SubElement(xlayer,'segment')
                cell = etree.SubElement(xseg,'name')
                cell.text = B[_name][idx].name
                index = etree.SubElement(xseg,'index')
                index.text = str(B[_name][idx].index)
                centx = etree.SubElement(xseg,'centx')
                centx.text = str(int(B[_name][idx].cent[0]))
                centy = etree.SubElement(xseg,'centy')
                centy.text = str(int(B[_name][idx].cent[1]))
                area = etree.SubElement(xseg,'area')
                area.text = str(int(B[_name][idx].area))
                length = etree.SubElement(xseg,'length')
                length.text = str(int(B[_name][idx].boundary_length))

    xml_out = etree.tostring(tree,pretty_print=False)
    with open(params.fout,'wb') as fout:
        fout.write(xml_out)
    print('Finished!')

    
