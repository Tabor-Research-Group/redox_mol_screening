#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
import numpy as np
import os

def main():
    mono_xyz = open("monomer.xyz")
    coord, center = get_centered_coord(mono_xyz)

    #Define the distance between two furthest atoms as "diameter"
    diameter = 0
    for i in range(len(coord)):
        v = [coord[i][1],coord[i][2],coord[i][3]]
        for j in range(len(coord)):
            if j >= i:
                v_2 = [coord[j][1],coord[j][2],coord[j][3]]
                tmp_dia = sum([(v[i] - v_2[i])**2 for i in range(3)])**0.5
                if tmp_dia > diameter:
                    diameter = tmp_dia

    #Align the normal vector of molecule plane to (0,0,-1)
    coord, vector = points_to_plane(coord)
    inner = 0
    while inner < 1:
        R, inner, vector = rotation(vector, np.array([0,0,-1]))
        for i in range(len(coord)):
            v = [[coord[i][1]],[coord[i][2]],[coord[i][3]]]
            for j in range(3):
                coord[i][j+1] = round(((R)*v)[j,0],4)

    work_dir = os.path.abspath('./')

    createfolder('./random_pairs')
    os.chdir('./random_pairs')
    for i in range(args.number):
        coord_2, R, v_dis = gen_random_coord(coord,vector,diameter)
        dimer = np.concatenate((coord,coord_2))

        #Write the generated dimer xyz
        with open(f'geometry_{i}.xyz','w') as geo:
            geo.write(f'{len(coord)*2}\n\n')
            for line in dimer:
                for k in range(0,4):
                    geo.write(f'{line[k]} ')
                geo.write('\n')

        #Write the rotation matrix and translation vector that generate the dimer
        with open(f'relative_coord_{i}.txt','w') as rel:
            rel.write('Rotation Matrix :\n\n')
            for j in range(3):
                rel.write(f'{R[j,0]:.4f},{R[j,1]:.4f},{R[j,2]:.4f}\n')
            rel.write(f'\nTranslation vector : {v_dis[0]:.4f},{v_dis[1]:.4f},{v_dis[2]:.4f}\n')

    os.chdir(work_dir)

#Extract and center the xyz of the monomer
def get_centered_coord(xyz):
    start = -2
    xyz_dict = dict()

    for atoms in xyz:
        if start >= 0:
            atoms = atoms.rstrip().split()
            atoms[1:] = [float(x) for x in atoms[1:]]
            xyz_dict[start] = atoms
        start += 1

    x,y,z = 0,0,0
    for i in xyz_dict:
        x += xyz_dict[i][1]
        y += xyz_dict[i][2]
        z += xyz_dict[i][3]

    x = x/len(xyz_dict)
    y = y/len(xyz_dict)
    z = z/len(xyz_dict)
    center = [x,y,z]
    center = [round(center[i],4) for i in range(3)]
    for i in xyz_dict:
        for j in range(3):
            xyz_dict[i][j+1] = round(xyz_dict[i][j+1] - center[j],4)
    
    coord = np.array([xyz_dict[i][:] for i in xyz_dict], dtype=object)
    xyz.close()

    return coord, center
        

#ref : https://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
#Given points, return normal vector of regression plane & moving points to ax+by-z = 0 plane
def points_to_plane(coord):
    tmp_A = []
    tmp_b = []

    for atom in coord:
        if atom[0] != 'H':
            tmp_A.append([atom[1],atom[2],1])
            tmp_b.append(atom[3])

    b = np.matrix(tmp_b).T
    A = np.matrix(tmp_A)
    fit = (A.T * A).I * A.T * b

    for atom in range(len(coord)):
        coord[atom][3] = round((coord[atom][3] - fit[2,0]),4)

    vector = [fit[0,0],fit[1,0],-1]
    vector = [round(i,4) for i in vector]

    return coord, vector



#ref: https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
#Given two vectors (x,y) return rotation matrix that rotates x to y
def rotation(x, y):    
    x_length = (x[0]**2 + x[1]**2 + x[2]**2)**0.5
    n_x = [(x[i]/x_length) for i in range(3)]
    
    inner = np.inner(n_x,y)
    cross = np.cross(n_x,y)
    identity = np.identity(3)
    vx = [[0,-cross[2],cross[1]],
          [cross[2],0,-cross[0]],
          [-cross[1],cross[0],0]]
    vx = np.matrix(vx)
    
    R = identity + vx + np.dot(vx,vx)*(1/(1+inner))
    new_vec = R*(np.matrix(n_x)).T
    new_vec = [(new_vec.tolist())[i][0] for i in range(3)]

    return R, inner, new_vec

#Check if the distance between the nearest atoms smaller than tolerance distance
def check(coord1, coord2, tolerance):
    dist = 1000
    for i in coord1:
        for j in coord2:
            tmp_dist = sum([(i[k+1] - j[k+1])**2 for k in range(3)])**(0.5)
            if tmp_dist < dist:
                dist = tmp_dist
    if dist > tolerance:
        return True
    else:
        return False


def gen_random_coord(coord, vector, diameter):
    while True:
        v = np.random.uniform(-1, 1, 3)
        v_rot = v / np.linalg.norm(v)
        v = np.random.uniform(-1, 1, 3)
        length = np.random.uniform(0, diameter, 1)
        v_dis = v*length / np.linalg.norm(v)
        coord_2 = np.copy(coord)
        R, inner, vector = rotation(vector,v_rot)

        for i in range(len(coord_2)):
            v = [[coord_2[i][1]],[coord_2[i][2]],[coord_2[i][3]]]
            for j in range(3):
                coord_2[i][j+1] = round((((R)*v)[j,0] + v_dis[j]),4)

        if check(coord,coord_2,2) and (not check(coord,coord_2,3.5)):
            break

    return coord_2, R, v_dis


def createfolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Random dimer generator from monomer.xyz')
    parser.add_argument("-n", "--number", help="Number of dimer to be generated.", type=int)
    args = parser.parse_args()
    
    if len(sys.argv) < 2:
        parser.print_usage()
        sys.exit(1)

    main()


