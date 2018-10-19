# Physics module

import numpy as np
from math import *
from mathf import *
from numpy import pi


# Constants

MU_EARTH = 3.986004418e14
MU_PLUTO = 8.81e11
MU_SUN = 1.327e20
SOLAR_MASS = 1.98848e30
EARTH_MASS = 5.9724e24 
PLUTO_MASS = 0.01303e24
GRAVI_CON = 6.67408e-11

'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
																				
     TODO: TOrbit involve parabolic, hyperbolic cases , orbital parameters still missing	
     
     add keywords to __init__ like mass,.... 
     
     MU_EARTH generalize to all planets
     
     conversion function between kilometers , meters , miles ,...
    																
																	    	
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''														


class TOrbit:

	def __init__( self , r_init , v_init ):
	
		# magnitude of r_init 
		self.r_init = np.linalg.norm( r_init ) 
		
		# magnitude of v_init 
		self.v_init = np.linalg.norm( v_init )
		
		# radial component of v 
		self.v_r_init = np.dot( r_init , v_init ) / self.r_init
		
		# reciprical of semimajor axis , alpha > 0 ellipse , alpha = 0 parabola , alpha < 0 hyperbola
		self.alpha = 2.0 / self.r_init - self.v_init**2 / MU_EARTH
		
		# specific angular momentum
		self.k_vec = np.cross( r_init , v_init )
		
		# eccentricity vector
		self.e_vec = np.cross( v_init , self.k_vec ) / MU_EARTH - r_init / np.linalg.norm( r_init )
		
		# perigee distance
		self.r_peri = np.dot( self.k_vec , self.k_vec ) / ( MU_EARTH * ( 1 + np.linalg.norm( self.e_vec ) ) )
		
		# apogee distance
		self.r_apo = np.dot( self.k_vec , self.k_vec ) / ( MU_EARTH * ( 1 - np.linalg.norm( self.e_vec ) ) )
		
		# eccentricity
		self.e = np.linalg.norm( self.e_vec ) 
		
		# semimajor axis
		self.a = ( self.r_apo + self.r_peri ) / 2.0 
		
		# semiminor axis
		self.b = self.a * np.sqrt( 1 - self.e**2 )
		
		# Energy
		self.Energy = - MU_EARTH / ( 2 * self.a )
		
		# Inclination
		self.i = np.arccos( self.k_vec[2] / np.linalg.norm( self.k_vec ) )
		
		# Period
		self.P = 2.0 * pi * pow( self.a , 3.0 / 2.0 ) / np.sqrt( MU_EARTH ) 	
		
		# magnitude of specific angular momentum 
		self.k = np.linalg.norm( self.k_vec )
		
		# node line 
		self.N_vec = np.cross( [ 0 , 0 , 1 ] , self.k_vec ) 
		self.N = np.linalg.norm( self.N_vec )
		
		# right ascension
		if self.N_vec[1] >=  0:
			self.Omega = np.arccos( self.N_vec[0] / self.N ) 
		else:
			self.Omega = 2 * pi - np.arccos( self.N_vec[0] / self.N ) 
			
		# argument of perigee
		
		if self.e_vec[2] >= 0: 
			self.omega = np.arccos( np.dot( self.N_vec , self.e_vec ) / ( self.N * self.e ) )
		else: 
			self.omega = 2 * pi - np.arccos( np.dot( self.N_vec , self.e_vec ) / ( self.N * self.e ) )
	
	
	# calculate orbital equation in plane perpendicular to k_z 
	# returns array of r and phi values

	def get_pol_eq ( self ):
		q = []
		w , r = [] , []
		phi = np.arange( 0.0 , 2*pi , 0.01 )
		
		for i in phi:
			w.append( [ np.cos( i ) , np.sin( i ) , 0 ] )
		rt = np.dot( self.k_vec , self.k_vec ) / ( MU_EARTH * ( 1 + np.linalg.norm( self.e_vec ) * np.cos( phi ) ) ) 
		for i in range( phi.size ): 
			r.append( rt[i] * np.asarray(w)[i] )
		
		v = ( MU_EARTH / self.k ) * np.array( [ - np.sin( phi ) , self.e + np.cos( phi ) , 0 ]  )
		Qt = [ [ - np.sin( self.Omega ) * np.cos( self.i ) * np.sin( self.omega ) + np.cos( self.Omega ) * np.cos( self.omega ) , - np.sin( self.Omega ) * np.cos( self.i ) * np.cos( self.omega ) - np.cos( self. Omega ) * np.sin( self.omega ) , np.sin( self.Omega ) * np.sin( self.i ) ] , [ np.cos( self.Omega ) * np.cos( self.i ) * np.sin( self.omega ) + np.sin( self.Omega ) * np.cos( self. omega ) , np.cos( self.Omega ) * np.cos( self.i ) * np.cos( self.omega) - np.sin( self.Omega ) * np.sin( self.omega ) , - np.cos( self.Omega ) * np.sin( self.i ) ] , [ np.sin( self.i ) * np.sin( self.omega ) , np.sin( self.i ) * np.cos( self.omega ) , np.cos( self.i ) ] ] 
		Q = np.matrix( Qt )
		for it in r:
			q.append( Q.dot( it ) )
		return np.squeeze(np.asarray(q)) , phi

	# determine true anomaly ( angle between perigee axis and current position )
	# TODO: What happens if r is an array? r in polar coordinates?
	# returns angle in radiant between r and perigee axis

	def get_true_ano ( self , r ):
		if self.v_r_init >= 0:
			self.true_ano = np.arccos( np.dot( self.e_vec , r ) / ( self.e * np.linalg.norm( r ) ) )
			return self.true_ano 
		else:
			self.true_ano = 2 * pi - np.arccos( np.dot( self.e_vec , r ) / ( self.e * np.linalg.norm( r ) ) )
			return self.true_ano
	
	# determine eccentric anomaly ( angle between line when r is projected on circle with radius equal to the semimajor axis and perigee axis )
	# TODO: What happens if r is in polar coordinates? 
	# returns eccentric anomaly in radiants
	def get_ecc_ano ( self , r ):
		return np.arccos( ( self.a - np.linalg.norm( r ) ) / ( self.a * self.e ) ) 
	
	# determine mean anomaly (  angle between line when r is projected on circle with radius equal to the semimajor axis going with constant speed and perigee axis	)
	def get_mean_ano ( self , E ):
		return E - self.e * np.sin( E ) 	
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
