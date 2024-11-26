import argparse
import numpy as np
import os 
from scipy.spatial import ConvexHull


from SIMULATION_library.file_utils import read_pts, write_pts, write_vtx, read_vtx
from common_4ch.mesh_utils import extract_tags
from common_4ch.file_utils import load_json

def compute_rotation_matrix(heartFolder, chamber, input_tags,biv_pts,biv_uvc_z,biv_uvc_v):
    ## Reorient the mesh such that the Z coordinate of the apex is 0 and the base is 1


    ### Find LV apex
    if chamber == "LA":
        condition = (biv_uvc_v == -1)
    elif chamber == "RA":
        condition = (biv_uvc_v == 1)

    relevant_uvc_z = biv_uvc_z[condition]

    apex_index_biv = np.argmin(relevant_uvc_z)

    apex_pts = biv_pts[condition][apex_index_biv]

    ### Find valve centre

    if chamber == "LA":
        valve = "MV"
    elif chamber == "RA":
        valve = "TV"

    tags_list_valve = extract_tags(input_tags,[valve])

    cmd = f"meshtool extract mesh -msh={heartFolder}/meshing/myocardium_OUT/myocardium -submsh={heartFolder}/surfaces_uvc/tmp/{valve} -ofmt=carp_txt -ifmt=carp_txt -tags={tags_list_valve[0]}"

    os.system(cmd)

    valve_pts = read_pts(f"{heartFolder}/surfaces_uvc/tmp/{valve}.pts")

    # valve_centre_pts = np.mean(valve_pts, axis=0)

    hull = ConvexHull(valve_pts)
    hull_pts = valve_pts[hull.vertices]  # Points forming the convex hull
    valve_centre_pts = np.mean(hull_pts, axis=0)

    write_pts(np.array([apex_pts]), f"{heartFolder}/surfaces_uvc_{chamber}/{chamber.lower()}/biv_apex.pts")
    
    # # Calculate the vector from apex to base
    vector_apba = valve_centre_pts - apex_pts


    rotation_angle = np.arctan2(vector_apba[2], vector_apba[1])

    rotation_matrix = np.array([[0, 0, 1],
                                [np.cos(rotation_angle), -np.sin(rotation_angle), 0],
                                
                                [np.sin(rotation_angle), np.cos(rotation_angle), 0]
                                
                                ])

    return rotation_matrix, valve_centre_pts

def rotate_mesh(mesh_pts, rotation_matrix, fixed_point):

    mesh_pts_rotated = np.dot(rotation_matrix, (mesh_pts - fixed_point).T).T + fixed_point

    return mesh_pts_rotated

def find_apex_SVC(input_tags, heartFolder, LV_apex_pts, RA_pts):

    tags_list_RA_SVC = extract_tags(input_tags,["RA","SVC_ring"])

    tags_RA_SVC_str = ':'.join([str(t) for t in tags_list_RA_SVC])

    cmd = f"meshtool extract surface -ofmt=carp_txt -ifmt=carp_txt -msh={heartFolder}/meshing/myocardium_OUT/myocardium -surf={heartFolder}/meshing/myocardium_OUT/tmp/RA_SVC -op={tags_RA_SVC_str}"

    os.system(cmd)
    

    RA_SVC_pts = read_pts(f"{heartFolder}/meshing/myocardium_OUT/tmp/RA_SVC.surfmesh.pts")

    distance_vector = []

    for i in range(len(RA_SVC_pts)):
        distance_vector.append(np.linalg.norm(RA_SVC_pts[i,:] - LV_apex_pts))
    
    farthest_pt = RA_SVC_pts[np.argmax(distance_vector), :]

    second_candidate_idx = find_closest_idx(RA_pts, farthest_pt)

    return second_candidate_idx

def find_apex(heartFolder,input_tags,chamber,biv_pts,biv_uvc_z,biv_uvc_v, MV_centre_pts = None):

    atrium_pts = read_pts(f"{heartFolder}/surfaces_uvc_{chamber}/{chamber.lower()}/{chamber.lower()}.pts")

    rotation_matrix, valve_centre_pts = compute_rotation_matrix(heartFolder=heartFolder, chamber=chamber, input_tags=input_tags,biv_pts=biv_pts,biv_uvc_z=biv_uvc_z,biv_uvc_v=biv_uvc_v)

    atrium_pts_rotated = rotate_mesh(mesh_pts=atrium_pts, rotation_matrix=rotation_matrix, fixed_point=valve_centre_pts)

    write_pts(atrium_pts_rotated, f"{heartFolder}/surfaces_uvc_{chamber}/{chamber.lower()}/{chamber.lower()}.rotated.pts")
    write_pts(np.array([valve_centre_pts]), f"{heartFolder}/surfaces_uvc_{chamber}/{chamber.lower()}/valve_centre.pts")

    atrium_apex_idx = np.argmax(atrium_pts_rotated[:,2])


    if chamber == "RA":
        ## Sometimes the apex of the RV is too close to the LV so the vector apex - TV is tilted towards the RA free wall. As an alternative, we compute a second candidate for apex: It will be the farthest point from the LV apex in the intersection of the SVC and the LA.

        LV_apex_pts = biv_pts[biv_uvc_v == -1][np.argmin(biv_uvc_z[biv_uvc_v == -1])]

        second_candidate_idx = find_apex_SVC(input_tags = input_tags, heartFolder = heartFolder, LV_apex_pts =LV_apex_pts, RA_pts = atrium_pts)


        ## Now we decide which one is the final one based on the min distance to the LV (using the MV centre).

        if np.linalg.norm(atrium_pts[second_candidate_idx] - MV_centre_pts) < np.linalg.norm(atrium_pts[atrium_apex_idx] - MV_centre_pts):
            atrium_apex_idx = second_candidate_idx


    write_vtx(vtx=np.asarray([atrium_apex_idx]), filename=f"{heartFolder}/surfaces_uvc_{chamber}/{chamber.lower()}/{chamber.lower()}.lvapex.vtx")

    return atrium_pts, valve_centre_pts

def find_closest_idx(array, target_point):
    distances = np.linalg.norm(array - target_point, axis=1)
    closest_index = np.argmin(distances)

    return closest_index

def select_landmarks(heartFolder, input_tags):

    # Common files used by both

    biv_pts = read_pts(f"{heartFolder}/surfaces_uvc/BiV/BiV.pts")
    biv_uvc_z = np.genfromtxt(f"{heartFolder}/surfaces_uvc/BiV/uvc/BiV.uvc_z.dat")
    biv_uvc_v = np.genfromtxt(f"{heartFolder}/surfaces_uvc/BiV/uvc/BiV.uvc_ven.dat")
    input_tags = load_json(input_tags)

    LA_pts, MV_centre_pts = find_apex(heartFolder,input_tags,"LA",biv_pts,biv_uvc_z,biv_uvc_v)
    RA_pts, _ = find_apex(heartFolder,input_tags,"RA",biv_pts,biv_uvc_z,biv_uvc_v, MV_centre_pts = MV_centre_pts)

    ##### Extract the atrial septum and find the central point
    os.makedirs(f"{heartFolder}/meshing/myocardium_OUT/tmp",exist_ok=True)

    tags_list_atria = extract_tags(input_tags,["LA","RA"])

    tags_atria_str = ':'.join([str(t) for t in tags_list_atria])

    cmd = f"meshtool extract surface -ofmt=carp_txt -ifmt=carp_txt -msh={heartFolder}/meshing/myocardium_OUT/myocardium -surf={heartFolder}/meshing/myocardium_OUT/tmp/atrial_septum -op={tags_atria_str}"

    os.system(cmd)

    atrial_septum_pts = read_pts(f"{heartFolder}/meshing/myocardium_OUT/tmp/atrial_septum.surfmesh.pts")

    atrial_septum_centre_pts = np.mean(atrial_septum_pts, axis=0)

    # To make sure we are in the surface
    atrial_septum_centre_idx =  find_closest_idx(atrial_septum_pts, atrial_septum_centre_pts)

    atrial_septum_centre_pts = atrial_septum_pts[atrial_septum_centre_idx]

    #### Find the closest idx in the la and ra to the centre point

    LA_septum_idx = find_closest_idx(LA_pts, atrial_septum_centre_pts)

    write_vtx(vtx=np.asarray([LA_septum_idx]), filename=f"{heartFolder}/surfaces_uvc_LA/la/la.rvsept_pt.vtx")

    RA_septum_idx = find_closest_idx(RA_pts, atrial_septum_centre_pts)

    write_vtx(vtx=np.asarray([RA_septum_idx]), filename=f"{heartFolder}/surfaces_uvc_RA/ra/ra.rvsept_pt.vtx")