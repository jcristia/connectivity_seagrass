# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2019, Ole Baadshaug, MET Norway & Ron Saper, Carleton University Canada.
#
# Caveat: This copyright will not interfere with the open nature of OpenDrift and OpenBerg

import logging

import numpy as np
import sys
from scipy.interpolate import interp1d

try:
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    from matplotlib import animation
except:
    logging.info('Basemap is not available, can not make plots')


from opendrift.models.basemodel import OpenDriftSimulation
from opendrift.elements.elements import LagrangianArray
from opendrift.readers.basereader import pyproj, BaseReader



# Defining the iceberg element properties
class IcebergObj(LagrangianArray):
    """Extending LagrangianArray with variables relevant for iceberg objects.

    """
    
    # We add the properties to the element class
    variables = LagrangianArray.add_variables([
        ('wind_drift_factor', {'dtype': np.float32,				# The fraction of the wind speed at which an iceberg is moved
                       	     		'unit': '%',
                       	        	'default': 0.018}),
        ('water_line_length', {'dtype': np.float32,				# Iceberg size
        							'unit': 'm',
                       	        	'default': 90.5}),
        ('keel_depth', {'dtype': np.float32,					# Iceberg keel depth
        							'unit': 'm',
                       	        	'default': 60})])

class OpenBerg(OpenDriftSimulation):
    """The Deterministic iceberg model in the OpenDrift framework.

        Advects an iceberg with the current at all available depths and
        as a function of the wind vector according to the above and below
        water cross-sectional profile of the object.
    """

    ElementType = IcebergObj
    
    required_variables = ['x_wind', 'y_wind',
                          'x_sea_water_velocity', 'y_sea_water_velocity',
                          'land_binary_mask'] 

    fallback_values = {'x_wind': 0.0,
                       'y_wind': 0.0,
                       'x_sea_water_velocity': 0.0,
                       'y_sea_water_velocity': 0.0}

    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'red',
                     'missing_data': 'gray', 'stranded': 'blue'}

    # Configuration
    def __init__(self, d=None, label=None, *args, **kwargs):

        self.name = 'OpenBerg'
        self.label=label                                                                
                
        self.required_profiles = ['x_sea_water_velocity',
                                  'y_sea_water_velocity']  # Get vertical current profiles
                                  
        self.required_profiles_z_range = [-190, 0] # [min_depth, max_depth]

        # Calling general constructor of parent class
        super(OpenBerg, self).__init__(*args, **kwargs)


    def update(self):
        """Update positions and properties of icebergs."""     
		# Move icebergs at given % of wind speed
        self.advect_wind()

        # Move icebergs as per Anne Barker et al, 2004, with weighted average of current  
        # down to draft of berg.
        # For current 10n meters down, where n is the depth index from 0 to 19

        # RETRIEVE CURRENT SPEED:
        y_sea_water_vel = self.environment_profiles['y_sea_water_velocity']
        x_sea_water_vel = self.environment_profiles['x_sea_water_velocity']        

		# FIND THE WEIGHTED AVERAGE CURRENT SPEED ACROSS THE ICEBERG KEEL        
        net_x_swv = np.zeros(x_sea_water_vel.shape[1])
        net_y_swv = np.zeros(y_sea_water_vel.shape[1])
        for indx in range(len(self.uw_weighting)):
            net_y_swv = net_y_swv +y_sea_water_vel[indx,:]*self.uw_weighting[indx]
            net_x_swv = net_x_swv +x_sea_water_vel[indx,:]*self.uw_weighting[indx]
            
        self.update_positions(net_x_swv,net_y_swv)
            

    def prepare_run(self):
        """	Model spesific preparations.
         	Set the weighting for modelled current depths as per Table 5 of Barker 2004,
        	'Determination of Iceberg Draft, Mass and Cross-Sectional Areas',
        	Proceedings of The Fourteenth (2004) International Offshore and 
        	Polar Engineering Conference.
        	
        	NB! This version of OpenBerg does not allow for seeding of iceberg elements
        	of different sizes.
        	
        	Also controles that the model handles readers without block data correctly.
        """
        # Retrieve profile provided in z dimension by reader  
        variable_groups, reader_groups, missing_variables = \
         	self.get_reader_groups(['x_sea_water_velocity','y_sea_water_velocity']) 
        
        if len(reader_groups) == 0:
        	# No current data -> fallback values used
        	self.uw_weighting = np.array([1])
        	return 
        	
		# Obtain depth levels from reader:
        reader_name = reader_groups[0][0]
        profile = np.abs(np.ma.filled(self.readers[reader_name].z))
        
       # Make sure that interpolation is left to reader if no blocks are used . 
       # NB! This is a workaround, two additional points should be removed in basereader!      
        if self.readers[reader_name].return_block == False:
        	self.use_block = False
        
        # If current data is missing in at least one dimension, no weighting is performed:
        if len(missing_variables) > 0:
        	logging.warning('Missing current data, weigthing array set to [1]')
        	self.uw_weighting = np.array([1])
        	return       
        	  		
        # No need to create weighting array if only one z-level is provided:
        if len(profile) == 1:
        	self.uw_weighting = np.array([1])
        	return
        	        	
        # Make copies to prevent outside value to be modified
        water_line_length = self.elements_scheduled.water_line_length.copy()
        depth = self.elements_scheduled.keel_depth.copy()

        # Calculate the weighting array corresponding to the iceberg profile
        uw_weighting = self.composite_iceberg(water_line_length=water_line_length, depth=depth)
             

            
        #### Interpolate weighting array to match z-dimension of current reader ###
        
        # Array of "Barker-depths" (10n)m, , where n is the depth index from 0 to 19
        x = np.linspace(0,(len(uw_weighting)-1)*10,len(uw_weighting)) 
        
        # Create a linear interpolator
        interpol = interp1d(x,uw_weighting) 
        
        # Obtain corresponding depth levels from reader:        
        self.reader_z_profile = profile[(profile >= 0) & (profile <= x.max())]
        
        if len(profile[profile < 0]) > 0:
        	logging.warning('Current reader containing currents above water!'
        						'Weighting of current profile may not work!')
        
        # Interpolate and normalize weights:
        
        ###### NB!  Interpolator only includes values from reader within the range of the 
        ######		Barker-depth array. Meaning values just outside is not used for interpolation.
        ######		Ex.: If x.max=195 and the current reader includes data for 
        ######		the z-profile: [0,3,10,15,25,50,75,100,150,200] only data from the 
        ######		z-levels [0,3,10,15,25,50,75,100,150] are used.
        
        interpol_weight = interpol(self.reader_z_profile) 
        normalized_weight = interpol_weight/interpol_weight.sum()
        
        if normalized_weight.all():
            self.uw_weighting = normalized_weight
                              	
        super(OpenBerg, self).prepare_run()

    def composite_iceberg(self, water_line_length=90.5, depth=60):

    	"""This function creates a weigthing array for the current across the keel of an
	    	iceberg based on waterline length and keel depth. The function uses the parameters
    		in table 5 from Barker et. al.(2004).
    	"""
		# Parameters from Baker et. al.:
		
    	a_param = [9.5173,11.1717,12.4798,13.6010,14.3249,13.7432,13.4527,15.7579,
    				14.7259,11.8195,11.3610,10.9202,10.4966,10.0893,9.6979,9.3216,8.9600,
    				8.6124,8.2783,7.9571]
				
    	b_param = [-25.94,-107.50,-232.01,-344.60,-456.57,-433.33,-519.56,-1111.57,-1125.00,
    					-852.90,-931.48,-1007.02,-1079.62,-1149.41,-1216.49,-1280.97,
    					-1342.95,-1402.52,-1459.78,-1514.82]
				
    	d = int(round(depth/10))

    	if d < 0: #Incase depth is given as a negative value
    		d = -d 


    	if d > len(a_param):
    		d = len(a_param)
    		print('##### OpenBerg does not support icebergs with keel depths greater than 200m!\n' +
				'Using a composite iceberg with given waterline length and keel depth 200m')  

    	area_list=[]

    	for i in range(0,d):
	
    		A = self.lin_func(a_param[i],b_param[i],water_line_length)
    		area_list.append(A)

    	A_list = np.array(area_list)

    	# Normalize array such that it sums to 1:

    	weigthing_array = np.array(area_list)/sum(np.array(area_list))

    	return weigthing_array
	
	
    def lin_func(self,a,b,L):
    	"""Returns value of linear function A=aL+b."""
    	A = a*L + b
	
    	return A
