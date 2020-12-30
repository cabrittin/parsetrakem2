"""
@name: tree.py
@description:
Module for trakem2

Used to manipulate trakem2 file

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2019-12-05
"""
from lxml import etree
from tqdm import tqdm
import numpy as np

def fix_calibration(tree):
    root = tree.getroot()
    t2_layer_set = root.find('t2_layer_set')
    calibration = t2_layer_set.find('t2_calibration')
    calibration.set('pixelWidth','1.0')
    calibration.set('pixelHeight','1.0')
    calibration.set('pixelDepth','20.0')
    calibration.set('unit','pixel')

def extract_area_lists(tree,area_list_title):
    """
    Removes area list from TrakEM2 tree not in 
    area_lists_title

    Parameters:
    -----------
    tree: lxml.etree.parse  file
    area_list_title: (list, str) Area lists to keep

    """
    if isinstance(area_list_title,str):
        area_list_title = [area_list_title]
    
    root = tree.getroot()
    project = root.find("project")
    for neuron in project.findall('neuron'):
        title = neuron.get('title')
        if title not in area_list_title:
            project.remove(neuron)

    t2_layer_set = root.find('t2_layer_set')
    for area in t2_layer_set.findall('t2_area_list'):
        title = area.get('title')
        if title not in area_list_title:
            t2_layer_set.remove(area)

def restrict_segments_by_oid(tree,oids):
    """
    Remove segments from the tree not in the roi

    Parameters:
    -----------
    tree: lxml.etree.parse  file
    oids: layer ids
    """
    root = tree.getroot()
    t2_layer_set = root.find('t2_layer_set')
    for area in tqdm(t2_layer_set.findall('t2_area_list'),desc='Restrict oids'):
        for t2_area in area.findall('t2_area'):
            layer_id = t2_area.get('layer_id')
            if layer_id not in oids:
                area.remove(t2_area)

def restrict_segments_by_roi(tree,roi,if_check=False):
    """
    Remove segments from the tree not in the roi

    Parameters:
    -----------
    tree: lxml.etree.parse  file
    oids: layer ids
    roi: list, [xmin (int),ymin (int),xmax (int),ymax (int)]
    """
    xmin,ymin,xmax,ymax = roi
    root = tree.getroot()
    t2_layer_set = root.find('t2_layer_set')

    for area in tqdm(t2_layer_set.findall('t2_area_list'),desc='Restric rois'):
        trans = get_area_list_transform(area)
        for t2_area in area.findall('t2_area'):
            layer_id = t2_area.get('layer_id')
            #patch = root.xpath("//t2_patch/@oid=")
            t2_path = t2_area.findall('t2_path')
            num_paths = len(t2_path)
            paths_removed = 0
            for path in t2_path:
                p = path.get('d')
                p = path_to_array(p)
                p[:,0] += trans[0]
                p[:,1] += trans[1]
                pnew = [_p for _p in p.tolist() if in_roi(_p,roi)]
                if pnew:
                    pnew = np.array(pnew)
                    pnew[:,0] -= trans[0]
                    pnew[:,1] -= trans[1]
                    pnew = array_to_path(pnew)
                    path.set('d', pnew)
                else:
                    t2_area.remove(path)
                    paths_removed += 1
                #cent = np.mean(p,axis=0)
                #c1 = (cent[0] < xmin) or (cent[0] > xmax)
                #c2 = (cent[1] < ymin) or (cent[1] > ymax)
                #if c1 or c2:
                    #t2_area.remove(path)
                    #paths_removed += 1
            if paths_removed == num_paths: area.remove(t2_area)

               
def get_area_list_transform(arealist):
    transform = arealist.get('transform').replace(')','')
    transform = transform.split(',')
    return [int(float(transform[-2])),int(float(transform[-1]))]

def path_to_array(path):
    """
    Converts path to numpy array

    Parameters
    ----------
    path : path string
    
    Returns
    ---------
    path : list
      Return a tranformed boundary object path  

    """
    path = path.replace('M ','')
    path = path.replace(' z','')
    path = path.split(' L ')
    return np.array([list(map(float,p.split(' '))) for p in path])
    #path = np.array([list(map(float,p.split(' '))) for p in path])
    #path[:,0] += translate[0]
    #path[:,1] += translate[1]
    #return path

def array_to_path(parray):
    path = [' '.join(map(str,p)) for p in parray.tolist()]
    return 'M ' + ' L '.join(path) + ' z'


def in_roi(pt,roi):
    x,y = pt
    xmin,ymin,xmax,ymax = roi
    c1 = (x >= xmin) and (x <= xmax)
    c2 = (y >= ymin) and (y <= ymax)
    return c1 and c2
 
def change_mipmaps_dir(tree,dname):
    """
    Change the mimpas directory of trakem2 file

    Parameters:
    -----------
    tree: lxml.etree.parse file
    dname: str, new directory of mipmaps folder

    """
    root = tree.getroot()
    project = root.find("project")
    project.set('mipmaps_folder',dname)

def add_synapse_area_lists_to_tree(tree,syn):
    """
    Add synapse area list to tree

    Parameters:
    -----------
    tree: lxml.etree.parse file
    syn: dict,
        { 'x': x coordinate
          'y': y coordinate
          'z': image number
          'layer_id': TrakEM2 layer id
          'svg': shape of synapse in svg format
        }

    """
    root = tree.getroot()
    project = root.find("project") 
    t2_layer_set = root.find("t2_layer_set") 
    aid = 10000
    nid = 2*aid
    for (contin,s) in syn.items():
        nid += 1
        aid += 2
        oid = aid - 1
        
        #Add project neuron
        neuron = etree.SubElement(project,'neuron')
        neuron.set('id',str(nid))
        neuron.set('title',s.get_title())
        neuron.set('expand','false')
        area_list = etree.SubElement(neuron,'area_list')
        area_list.set('oid',str(oid))
        area_list.set('id',str(aid))

        #Add t2_area_list
        t2_area_list = etree.SubElement(t2_layer_set,'t2_area_list')
        t2_area_list.set('oid',str(oid))
        t2_area_list.set('width',s.width) 
        t2_area_list.set('height',s.height) 
        t2_area_list.set('transform',
                    "matrix(1.0,0.0,0.0,1.0,%d.0,%d.0)"%(s.get_min_x(),s.get_min_y()))
        t2_area_list.set('title',s.get_title())
        t2_area_list.set('links',"")
        t2_area_list.set('layer_set_id',"3")
        t2_area_list.set("fill_paint","true")
        t2_area_list.set("style","stroke:none;fill-opacitiy:1.0;fill:%s;"%s.color)
        for (o,v) in s.objects.items():
            #print(contin,s.get_title(),v.x,v.y,v.img,v.svg)
            t2_area = etree.SubElement(t2_area_list,'t2_area')
            t2_area.set('layer_id',v.layer_id)
            t2_path = etree.SubElement(t2_area,'t2_path')
            t2_path.set('d',v.svg)

