"""
@name: import_synapses.py 

@description

Import synapse data into trakem2 file

@author: Christopher Brittin
@email: 'cabrittin' + <at> + 'gmail' + '.' + 'com'
"""
import os
import argparse
from configparser import ConfigParser,ExtendedInterpolation
from lxml import etree
import networkx as nx
import numpy as np

#Local modules
import aux
from parsetrakem2.parse import ParseTrakEM2
from parsetrakem2.tree import extract_area_lists, change_mipmaps_dir, add_synapse_area_lists_to_tree
import db as db

W = "100.0"
H = "100.0"
#SVG = "M 10 0 L 0 10 L 10 20 L 19 11 L 19 10 L 20 10 L 10 0 z"
SVG = "M 10 10 L 10 90 L 90 90 L 90 10 z"
#SVG = "M 10 10 V 90 V 90 H 10 Z"

CONFIG = os.environ['CONFIG']

class Synapse:
    def __init__(self,contin,stype,partners):
        self.contin = contin 
        self.stype = stype
        self.partners = partners
        self.objects = {}
        self.color = '#ffff00'

    def add_object(self,key,obj):
        self.objects[key] = objectview(obj)

    def get_title(self):
        partners = ','.join(self.partners)
        return str(self.contin) + ': ' + self.img + ' ' +  self.stype + ' ' + partners
    
    def get_min_x(self):
        return min([v.x for (k,v) in self.objects.items()])

    def get_min_y(self):
        return min([v.y for (k,v) in self.objects.items()])

class objectview(object):
    def __init__(self, d):
        self.__dict__ = d

def clean_graph(A,G):
    edges = [(a,b) for (a,b) in G.edges()]
    for (a,b) in edges:
        if not A.has_edge(a,b):
            G.remove_edge(a,b)

def reflect_graph(G,lrdict):
    H = nx.Graph()
    if G.is_directed(): H = nx.DiGraph()
    for (a,b) in G.edges():
        H.add_edge(lrdict[a],lrdict[b],weight=G[a][b]['weight'])

    try:
        for n in H.nodes():
            H.node[n]['top_class'] = G.node[lrdict[n]]['top_class']
    except:
        pass
    
    return H

def screen_presynapses(cell,G,synapses):
    for s in synapses:
        if s[0] != cell: continue
        post = [p for p in s[1].split(',') if G.has_edge(cell,p)]
        if not post: continue
        contin = s[3]
        #print(s)
        yield Synapse(contin,'pre',post)

def screen_postsynapses(cell,G,synapses):
    for s in synapses:
        pre = s[0]
        post = s[1].split(',')
        #if cell in post: print(s)
        if cell not in post: continue
        if not G.has_edge(pre,cell):continue
        contin = s[3]
        #print(s)
        yield Synapse(contin,'post',[pre])

def screen_gap_junctions(cell,G,synapses):
    for s in synapses:
        p = []
        if s[0] == cell: p = [s[1]]
        if s[1] == cell: p = [s[0]]
        if not p: continue
        if not G.has_edge(s[0],s[1]): continue
        contin = s[3]
        yield Synapse(contin,'gap',p)

db2cfg = {'N2U':'trakem2_n2u','JSH':'trakem2_jsh'}


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('cell',
                        action = 'store',
                        help = 'Cell class')
     
    parser.add_argument('dataset',
                        action = 'store',
                        help = 'Dataset [N2U,JSH]')

    parser.add_argument('-c','--config',
                        dest = 'config',
                        action = 'store',
                        default = CONFIG,
                        required = False,
                        help = 'Config file')
    
    parser.add_argument('--bundle',
                        action='store_true',
                        default=False,
                        help='Color by intra-/inter- bundle')
    
    parser.add_argument('--pre',
                        action='store_true',
                        default=False,
                        help='Import presynapses')
    
    parser.add_argument('--post',
                        action='store_true',
                        default=False,
                        help='Import postsynapses')
    
    parser.add_argument('--gap',
                        action='store_true',
                        default=False,
                        help='Import gap junctions')
    
    
    params = parser.parse_args()

    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(params.config)

    tkey = db2cfg[params.dataset]
    
    layers = dict([(l[0],[l[1],int(l[2]),int(l[3])]) 
                        for l in aux.read.into_list2(cfg[tkey]['layer_keys'])])

    right = aux.read.into_list(cfg['mat']['right_nodes'])
    lrdict = aux.read.into_lr_dict(cfg['mat']['lrmap'])

    A = nx.read_graphml(cfg['refgraphs']['adj']%4)
    C4 = nx.read_graphml(cfg['refgraphs']['chem']%4)
    C3 = nx.read_graphml(cfg['refgraphs']['chem']%3)
    G4 = nx.read_graphml(cfg['refgraphs']['gap']%4)
    G3 = nx.read_graphml(cfg['refgraphs']['gap']%3)
    
    C = nx.compose(C3,C4)
    G = nx.compose(G3,G4)
    
    if params.cell in right:
        A = reflect_graph(A,lrdict)
        C = reflect_graph(C,lrdict)
        G = reflect_graph(G,lrdict)
    
    C.add_edge('SMBDR','RMED',weight=1)
    C.add_edge('SMBVR','RMEV',weight=1)
    con = db.connect.default(cfg[tkey]['db'])
    cur = con.cursor()

    #Load synapses
    syn = {}
   
    if params.pre or params.post:
        itype = 'prepost'
        synapses = db.mine.get_synapse_data(cur,'chemical')
        if params.pre: 
            for s in screen_presynapses(params.cell,C,synapses): syn[s.contin] = s
        if params.post:
            for s in screen_postsynapses(params.cell,C,synapses): syn[s.contin] = s
   
    if params.gap:
        itype = 'gap'
        synapses = db.mine.get_synapse_data(cur,'electrical')
        for s in screen_gap_junctions(params.cell,G,synapses): syn[s.contin] = s
    
    #Color synapses
    if params.bundle:
        itype = 'bundle'
        bundles = aux.read.into_dict(cfg['clusters']['final_lr'])
        for (contin,s) in syn.items():
            s.color = cfg['trakem2_synapses']['in_bundle']
            for p in s.partners:
                if not A.has_node(p): continue
                if bundles[params.cell] != bundles[p]:
                    s.color = cfg['trakem2_synapses']['out_bundle']
                if (bundles[params.cell] == 'Unclassified' or bundles[p] == 'Unclassified'):
                    s.color = cfg['trakem2_synapses']['unclass_bundle']
    else:
        for (contin,s) in syn.items(): s.color = cfg['trakem2_synapses'][s.stype]
    
    #Format synapses
    syn_remove = []
    for (contin,s) in syn.items():
        syn_added = False
        for (key,x,y,z,_img) in db.mine.get_contin_xyz(cur,contin):
            img = _img.replace('.tif','')
            if img not in layers: continue
            s.img = img
            obj = {'x':x + layers[img][1],
                    'y':y + layers[img][2],
                    'z':z,
                    'img':img,
                    'layer_id':layers[img][0],
                    'svg':SVG} 
            s.add_object(key,obj)
            syn_added = True
        if not syn_added: syn_remove.append(contin)
        s.width = W
        s.height = H
    
    for c in syn_remove: del syn[c]
    
    tree = etree.parse(cfg[tkey]['trakem2'])
    
    #Setup directory
    if cfg.getboolean(tkey,'make_dir'):
        if not os.path.exists(cfg[tkey]['dname']): os.makedirs(cfg[tkey]['dname'])
        change_mipmaps_dir(tree,cfg[tkey]['mipmaps'])

    extract_area_lists(tree,params.cell)
    add_synapse_area_lists_to_tree(tree,syn)

    xml_out = etree.tostring(tree,pretty_print=False)
    with open(cfg[tkey]['fout']%(params.cell,itype),'wb') as fout:
        fout.write(xml_out)
    print('%s written to %s'%(params.cell,fout))
