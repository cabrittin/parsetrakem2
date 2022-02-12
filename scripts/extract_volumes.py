"""
@name: extract_volumes.py
@description:
    
Extracts the volumes into 3D arrays. Each cell is stored as a separate 
array in binary (npy) format.

@author: Christopher Brittin
@email: "cabrittin"+ <at>+ "gmail"+ "."+ "com"
@date: 2019-12-05
"""
import os
import argparse
from lxml import etree
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
from multiprocessing_on_dill.managers import BaseManager
from multiprocessing_on_dill import Process,Manager,Pool
import random

from parsetrakem2.parse import ParseTrakEM2
from pycsvparser import read

def chunk_list(lst,num_chunks): 
    lst_size = len(lst)
    for (idx,i) in enumerate(range(0,lst_size,num_chunks)):
        yield idx,lst[i:i+num_chunks]

def load_data(params):
    print('TrakEM2 file: %s' %params.trakem2)
    
    BaseManager.register('ParseTrakEM2', ParseTrakEM2)
    manager = BaseManager()
    manager.start()
    inst = manager.ParseTrakEM2(params.trakem2)
    inst.get_layers()
    inst.get_area_lists()
    
    layers = sorted(inst.get_layer_list())
    print('Extracted %d layers.' %(len(layers)))
    layers = [(i,l) for (i,l) in enumerate(layers)]
    
    cells = inst.get_cells()
    cells = sorted(cells)
    cells = [(c,int(i)) for (i,c) in enumerate(cells) if c not in ['Pharynx','Phi_Marker']] 
    #cells = cells[72:]
    #cells = dict([(c,int(i)) for (i,c) in enumerate(cells)]) 
    #random.shuffle(layers)
    
    return inst,layers,cells


def mp_chunk_slices(params):
    done_slices = []
    for fname in os.listdir(params.dout):
        fname = fname.split('_')[-1]
        idx = int(fname.split('.')[0])
        done_slices.append(idx)
    done_slices = sorted(done_slices)
    inst,layers,cells = load_data(params)
    layers = [l for l in layers if l[0] not in done_slices]
    #layers=layers[60:]
    #print(layers)
    #random.shuffle(layers)
    
    #cdx = read.into_dict('data/jsh_vol_cell_index.csv')
    cdx = read.into_dict('data/n2u_vol_cell_index.csv')
    for (k,v) in cdx.items(): cdx[k] = int(v)
    
    num_chunks = len(layers) // params.nproc
    
    procs = []
    for (job_id,_layers) in chunk_list(layers,num_chunks):
        proc = Process(target=worker, args=(job_id,_layers,cdx,cells,inst,params,))
        procs.append(proc)
        proc.start()

    for proc in procs: proc.join()
    
    print("\n" * (len(procs) + 1)) 

def worker(pid,layers,cdx,cells,P,params):
    tqdm_text = f'PID {pid}: layers'
    [height,width] = P.get_layer_dims()
    V = np.zeros((height,width),dtype=np.uint8)
    cell_names = sorted([c for c in cdx.keys()])
    #cdx = dict([(c[0],int(c[1])) for c in cells])

    for (ldx,lname) in layers:
        B = P.get_boundaries_in_layer(lname,area_thresh = params.area_thresh,
                                    scale_bounding_box = params.scale_bounding_box,
                                    area_lists=cell_names)
 
        tqdm_text = f'PID {pid}: layer:#{ldx}:{lname}'
        with tqdm(total=len(cell_names), desc=tqdm_text, position=pid+1) as pbar:
            for cell in cell_names:
                try: 
                    for (k,v) in B[cell].items():
                        M = v.get_global_display_matrix(height,width)
                        V[M>0] = cdx[cell]
                except:
                    pass

                pbar.update(1)
        
        fout = f'{params.dout}JSH_slice_{ldx}.npz'
        np.savez_compressed(fout,V=V)

def mp_fix_cell(params):
    inst,layers,cells = load_data(params)
    
    random.shuffle(layers)
    cdx = read.into_dict('data/jsh_vol_cell_index.csv')
    for (k,v) in cdx.items(): cdx[k] = int(v)
    
    num_chunks = len(layers) // params.nproc

    procs = []
    for (job_id,_layers) in chunk_list(layers,num_chunks):
        proc = Process(target=worker_fix, args=(job_id,_layers,cdx,params.fix_cell,inst,params,))
        procs.append(proc)
        proc.start()

    for proc in procs: proc.join()
    
    print("\n" * (len(procs) + 1)) 

def worker_fix(pid,layers,cdx,cell_fix,P,params):
    [height,width] = P.get_layer_dims()
    cell_names = sorted([c for c in cdx.keys()])
    #cdx = dict([(c[0],int(c[1])) for c in cells])

    tqdm_text = f'PID {pid}: {cell_fix} layers:'
    with tqdm(total=len(layers), desc=tqdm_text, position=pid+1) as pbar: 
        for (ldx,lname) in layers:
            vol = f'{params.dout}JSH_slice_{ldx}.npz'
            V = np.load(vol)['V'] 
            B = P.get_boundaries_in_layer(lname,area_thresh = params.area_thresh,
                                    scale_bounding_box = params.scale_bounding_box,
                                    area_lists=[cell_fix])
            try: 
                for (k,v) in B[cell_fix].items():
                    M = v.get_global_display_matrix(height,width)
                    #print(cdx[cell_fix])
                    V[M>0] = cdx[cell_fix]
            except:
                    pass

            pbar.update(1)
            np.savez_compressed(vol,V=V)


def mp_test(params):
    P,layers,cells = load_data(params)
    [height,width] = P.get_layer_dims()
    layer = layers[0]
    (ldx,lname) = layer 
    cell_names = [c[0] for c in cells]
    cdx = dict([(c[0],int(c[1])) for c in cells])
 
    B = P.get_boundaries_in_layer(lname,area_thresh = params.area_thresh,
                                    scale_bounding_box = params.scale_bounding_box)
    
    V = np.zeros((height,width),dtype=np.uint8) 
    for cell in cell_names:
        try: 
            print(cell)
            for (k,v) in B[cell].items():
                M = v.get_global_display_matrix(height,width)
                M = M*cdx[cell] 
                V[:,:] += M.astype(np.uint8)
        except:
            pass
            #print('skipping',cell)
    fig,ax = plt.subplots(1,1,figsize=(10,10))
    ax.imshow(V)
    ax.set_title(cell)
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('trakem2',
                        action="store",
                        help="TrakEM2 file")

    parser.add_argument('dout',
                        action = 'store',
                        help = "Output directory")

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
    
    parser.add_argument('--fix_cell',
                        dest = 'fix_cell',
                        action = 'store',
                        required = False,
                        default = None,
                        help = ("Cell name to fix")
                        )


    
    params = parser.parse_args()

    
    #mp_by_subvolume(params) 
    #mp_by_cells(params) 
    #mp_by_slice(params) 
    mp_chunk_slices(params) 
    #mp_fix_cell(params) 
    
    #mp_test(params)
    """
    from pycsvparser import write
    
    inst,layers,cells = load_data(params)
    cell_names = [c[0] for c in cells]
    cdx = dict([(c[0],int(c[1])+1) for c in cells])
    
    cout = [[c,cdx[c]] for c in cell_names]
    fout = 'data/n2u_vol_cell_index.csv'
    write.from_list(fout,cout)
    """
