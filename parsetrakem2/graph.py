"""
graph.py

Classes for handling TrakEM2 graph files.

@author Christopher Brittin
@date 2019-03-20
"""

from networkx import Graph as _Graph
import lxml.etree as etree

class Graph(_Graph):
    
    def __init__(self,sif,labels):
        _Graph.__init__(self)

        self.labels = {}
        with open(labels,'r') as fin:
            for f in fin.readlines():
                [oid,cell] = f.split('\t')
                cell = cell.split(' ')[0]
                self.labels[oid] = cell
       
        with open(sif,'r') as fin:
            for f in fin.readlines():
                [pre,post] = f.split(' pd ')
                post = post.split('\n')[0]
                pre = self.labels[pre]
                post = self.labels[post]
                if pre == post: continue
                self.add_edge(pre,post)
        
