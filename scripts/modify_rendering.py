"""
@name: modify_rendering.py

@description:

Modifies the rendering by keeping only the desired segments.

!!!!IMPORTANT!!!!
The output xml file will not have the appropriate header in order
to be read by TrakEM2. You will need insert the header yourself.
An accompanying header file is provided in mat/. In linux you
can combine the files with the 'cat' command e.g.

cat /mat/header.txt output.xml > trakem2_readable.xml

where output.xml is the output file produced by this script and
trakem2_readable is the file read by trakem2. 

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2019-12-05

Synopsis:
 python modify_rendering.py [CONFIG_FILE]

Parameters:
   CONFIG_FILE: path tho configuration (.ini) file. 



"""

import os
import argparse
from configparser import ConfigParser,ExtendedInterpolation
from lxml import etree

import aux
from parsetrakem2.parse import ParseTrakEM2
from parsetrakem2.tree import *


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('config',
                        action = 'store',
                        help = 'Desc')

    params = parser.parse_args()

    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(params.config)
        
    cells = aux.read.into_list(cfg['input']['cells'])
    #P = ParseTrakEM2(cfg['trakem2'])

    tree = etree.parse(cfg['input']['trakem2'])
    fix_calibration(tree)

    #Setup directory
    if cfg.getboolean('output','make_dir'):
        if not os.path.exists(cfg['output']['dname']): os.makedirs(cfg['output']['dname'])
        change_mipmaps_dir(tree,cfg['output']['mipmaps'])

    extract_area_lists(tree,cells)
    
    P = ParseTrakEM2(cfg['input']['trakem2'])
    P.get_area_lists()
    P.get_layers()
    layers = sorted(P.layers.keys())

    
    if cfg.getboolean('params','modify_Z'):
        zmin = cfg.getint('params','zmin')
        zmax = cfg.getint('params','zmax')
        layers = layers[zmin:zmax + 1]
        oids = [P.layers[l].oid for l in layers]
        restrict_segments_by_oid(tree,oids)      

    if cfg.getboolean('params','modify_ROI'):
        oids = [P.layers[l].oid for l in layers]
        roi = [int(i) for i in cfg['params']['roi'].split(',')]
        restrict_segments_by_roi(tree,roi,if_check=True)      
    
    #if cfg.getboolean('params','modify_ROI'):
    #    oids = [P.layers[l].oid for l in layers]
    #    restrict_segments_by_roi(tree,[2166,1000,4446,5472],if_check=True)      
    

    xml_out = etree.tostring(tree,pretty_print=False)
    with open(cfg['output']['fout'],'wb') as fout:
        fout.write(xml_out)

    if cfg.getboolean('params','modify_color'):
        nclass = aux.read.into_dict(cfg['input']['nclass'])
        color = aux.read.into_dict2(cfg['input']['color'])
        P = ParseTrakEM2(cfg['output']['fout'])
        P.get_area_lists()
        cols = []
        for n in tqdm(P.area_lists,desc='Setting colors'):
            if n in nclass:
                try:
                    cols.append([n,color[nclass[n]][0],float(color[nclass[n]][1])])
                except:
                    cols.append([n,color[nclass[n]][0]])
        P.set_fill(cols)
        P.xml.write(cfg['output']['fout'],pretty_print=True)

    print('Finished!')

