# Physics module

import numpy as np
from math import *
from mathf import *
from numpy import pi


# Constants

MU_EARTH = 3.986004418e14
SOLAR_MASS = 1.98848e30
EARTH_MASS = 5.9724e24 
#EARTH_RADIUS = ‎6.3781e6
PLUTO_MASS = 0.01303e24
GRAVI_CON = 6.67408e-11




class TOrbit:

	def __init__( self , r_init , v_init ):
		self.k_vec = np.cross( r_init , v_init )
		self.e_vec = np.cross( v_init , self.k_vec ) / MU_EARTH - r_init / np.linalg.norm( r_init )
		self.r_peri = np.dot( self.k_vec , self.k_vec ) / ( MU_EARTH * ( 1 + np.linalg.norm( self.e_vec ) ) )
		self.r_apo = np.dot( self.k_vec , self.k_vec ) / ( MU_EARTH * ( 1 - np.linalg.norm( self.e_vec ) ) )
		self.a = ( self.r_apo + self.r_peri ) / 2.0 
		self.b = self.a * np.sqrt( 1 - np.linalg.norm( self.e_vec ) )
		self.E = - MU_EARTH / ( 2 * self.a )
		self.i = np.arccos( self.k_vec[2] / np.linalg.norm( self.k_vec ) )
		self.P = 2.0 * pi * pow( self.a , 3.0 / 2.0 ) / np.sqrt( MU_EARTH ) 	

		self.e = np.linalg.norm( self.e_vec ) 
	

	def get_pol_eq ( self ):
		phi = np.arange( 0.0 , 2*pi , 0.01 )
		r = np.dot( self.k_vec , self.k_vec ) / ( MU_EARTH * ( 1 + np.linalg.norm( self.e_vec ) * np.cos( phi ) ) )
		return r , phi


	def get_true_ano ( self , r ):
		return np.arccos( np.dot( self.e_vec , r ) / ( self.e * np.linalg.norm( r ) ) )
	
	def get_ecc_ano ( self , r ):
		return np.arccos( ( self.a - np.linalg.norm( r ) ) / ( self.a * self.e ) ) 
		
	def get_mean_ano ( self , E ):
		return E - self.e * np.sin( e ) 	
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
