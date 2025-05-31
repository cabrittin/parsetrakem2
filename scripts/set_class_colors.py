"""
set_class_colors.py

Sets the fill colors for area_lists in the trakem2 xml file.

!!!!IMPORTANT!!!!
The output xml file will not have the appropriate header in order
to be read by TrakEM2. You will need insert the header yourself.
An accompanying header file is provided in mat/. In linux you
can combine the files with the 'cat' command e.g.

cat /mat/header.txt output.xml > trakem2_readable.xml

where output.xml is the output file produced by this script and
trakem2_readable is the file read by trakem2. 


Code has only been tested on Linux OS. If running on Windows there may 
be formatting issues with reading/writing to files. 

created: Christopher Brittin
date: 22 February 2018

Required 3rd party packages:
  argparse


Synopsis:
 python set_class_colors.py -t /path/to/trakem2 -n /path/to/class/file -c /path/to/color/code/file -o /path/to/output/xml

Parameters:
  -t, --trakem2 (str): Trakem2 file
  -n, --nclass  (str): Class file which maps area_list names to classes
  -c, --color   (str): Color code file maps colors to classes
  -o, --output  (str): Output file

"""
import argparse

#Local modules
from pycsvparser import read
from parsetrakem2.parse import ParseTrakEM2



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-t','--trakem2',
                        dest = 'trakem2',
                        action="store",
                        required= True,
                        default = None,
                        help="TrakEM2 file"
                        )

    parser.add_argument('-n','--nclass',
                        dest = 'nclass',
                        action="store",
                        required= True,
                        default = None,
                        help=("Area list class file. Should have "
                              "format:\narealist_name,class_name")
                        )

    parser.add_argument('-c','--color',
                        dest = 'color',
                        action="store",
                        required= True,
                        default = None,
                        help=("Color code file. Should have format "
                              "\nclass_name,html_color_code")
                        )

    parser.add_argument('-o','--out',
                        dest = 'fout',
                        action="store",
                        required= True,
                        default = None,
                        help="Output xml file"
                        )
    
    params = parser.parse_args()
    nclass = read.into_dict(params.nclass)
    color = read.into_dict(params.color,multi_dim=True)
    P = ParseTrakEM2(params.trakem2)
    P.get_area_lists()
    cols = []
    for n in P.area_lists:
        if n in nclass:
            try:
                cols.append([n,color[nclass[n]][0],float(color[nclass[n]][1])])
            except:
                cols.append([n,color[nclass[n]][0]])
    P.set_fill(cols)
    P.xml.write(params.fout,pretty_print=True)
