"""
boundary_extraction.py

Test parsetrakem2 module

created: Christopher Brittin
date: 17 October 2018

"""
import sys
sys.path.append(r'../trakem2/')
import argparse
import time
import random
import matplotlib.pyplot as plt
from itertools import combinations
import multiprocessing as mp

from parsetrakem2.parsetrakem2 import ParseTrakEM2

def batch_compute_adjacency(boundaries,P,output):
    P.batch_compute_adjacency(boundaries,output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description= "Test parsetrakem2 module")

    parser.add_argument('trakem2',
                        action='store',
                        help='TrakEM2 file')
    parser.add_argument('-v','--verbose',
                        action = 'store_true',
                        dest = 'verbose',
                        required = False,
                        default = False,
                        help = "Verbose")
    
    params = parser.parse_args()
    SCALE = 1.1
    
    print("Loading trakem2 file: %s" %params.trakem2)
    time0 = time.time()
    P = ParseTrakEM2(params.trakem2)
    print("Time to load file: %2.3f sec" %(time.time() - time0))
    P.get_layers()
    print('Extracted %d layers.' %(len(P.layers)))
    if params.verbose:
        print('Extracted the following layers:')
        for l in sorted(P.layers.keys()):
            print('\t%s' %l)
    P.get_area_lists()
    print('Extracted %d area lists.' %(len(P.area_lists)))
    if params.verbose:
        print('Extracted the following area lists:')
        for n in sorted(P.neurons):
            print('\t%s' %n)    
    L = sorted(P.layers.keys())[0]
    time0 = time.time()
    B = P.get_boundaries_in_layer(L)
    dt = time.time() - time0
    num_boundaries = 0
    for n in B: num_boundaries += len(B[n])
    print('Took %2.3f sec to extract %d boundaries in layer %s.'
          %(dt,num_boundaries,L))
    if params.verbose:
        print('Extracted the following boundaries from layer %s' %L)
        for n in sorted(B):
            print('\tCell: %s, boundaries %d' %(n,len(B[n])))

    n = random.choice(list(B.keys()))
    print('Randomly displaying bounding box(es) for %s' %n)
    for b in B[n]:
        A = B[n][b].get_display_matrix()
        plt.figure()
        plt.spy(A)
    plt.show()
    
    print('Scaling bounding boxes by %2.1f' %SCALE)
    time0 = time.time()
    for i in B:
        for b in B[i]:
            B[i][b].scale_bounding_box(SCALE)
    print('Time to scale all bounding boxes: %2.3f sec' %(time.time() - time0))
            
    for b in B[n]:
        A = B[n][b].get_display_matrix()
        plt.figure()
        plt.spy(A)
    plt.show()

    print('Pairwise overlap check')
    time0 = time.time()
    overlap = P.get_overlapping_boundaries(B)
    print('Time to complete overlap check is %2.3f sec.' %(time.time() - time0))
    print('Found %d pairs of overlaping bounding boxes.' %(len(overlap)))
    if params.verbose:
        print('Extracted overlapping bounding boxes')
        for n in sorted(overlap):
            print('\tCells: %s,%s' %n)    

    print('Computing adjacencies')
    """
    time0 = time.time()
    output = []
    for (b1,b2) in overlap:
        adj = P.compute_adjacency(b1,b2)
        if adj > 0:
            output.append((b1,b2,adj))
    print('Time to compute adjancencies %2.3f' %(time.time() - time0))
    print(len(output))
    """
    n = 2
    overlap_split = [overlap[i::n] for i in range(n)]
    time0 = time.time()
    pool = mp.Pool(processes = n)
    results = [pool.apply_async(P.batch_compute_adjacency,args=(o,))
               for o in overlap_split]
    output = [o for p in results for o in p.get()]
    print('Time to compute mp adjancencies %2.3f' %(time.time() - time0))
    print(len(output))
    
