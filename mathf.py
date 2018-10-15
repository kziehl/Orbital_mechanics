# Math module

import numpy as np
from math import *


# convert polar coordinates r and phi to cartesian coordinates x and y
# returns cartesian coordinates
def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    r = np.array( [ x , y , np.zeros( x.size ) ] )
    return r.T  # sorts x and y coordinates
# convert cartesian coordinates into polar coordinates
# returns polar coordinates
def cart2pol(x,y):
	rho = np.sqrt(x**2 + y**2)
	phi = atan2(x,y) # atan2 returns values between -pi and pi
	return (rho,phi)

# conversion from grad to radiant
def grad2rad(grad):
	return grad*2*np.pi/360.0


# conversion from radiant to grad
def rad2grad(rad):
	return rad*360.0/(2*np.pi)
	
	
# rotation matrix in 3D

def rot_by_phi( phi , r ):
	rt = []
	rot_mat = np.array( [ [ 1 , 0 , 0 ] , [ 0 , np.cos( phi ) , -np.sin( phi ) ] , [ 0 , np.sin( phi ) , np.cos( phi ) ] ] )
	for it in r:
		rt.append( rot_mat.dot( it ) ) 
	return rt
	
	
# get different coordinates of 3d position vecotor and return 3 different lists 


def get_xyz( r ):
	x , y , z = [] , [] , []
	
	for i in range( len( r ) ):
		x.append( r[i][0] )	
		y.append( r[i][1] )
		z.append( r[i][2] )
		
	return x,y,z
	
































	
