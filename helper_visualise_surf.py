import os 
import argparse

import common_4ch.meshtools_utils as mu
import common_4ch.file_utils as fu

def main(args):
    input_surf = args.input
    input_pts = args.pts

    folder = os.path.dirname(input_surf)
    basename = os.path.basename(input_surf)
    
    basename_noext = os.path.splitext(basename)[0]

    aux_name = f'aux_{basename_noext}'

    msh = os.path.join(folder, f'{aux_name}')
    fu.mycp(input_surf, f'{msh}.elem')
    fu.mycp(input_pts, f'{msh}.pts')

    cmd = f'meshtool convert -imsh={msh} -ifmt=carp_txt -omsh={msh} -ofmt=vtk'

    print(f'Folder: {folder}')
    print(f'Base name: {basename}')
    print(f'Base name no ext: {basename_noext}')
    print(f'Mesh: {msh}')
    
    print(cmd)
    os.system(cmd)

if __name__ == '__main__' : 
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default=None, help='Input surf file')
    parser.add_argument('--pts', type=str, default=None, help='Input pts file')
    args = parser.parse_args()

    main(args)