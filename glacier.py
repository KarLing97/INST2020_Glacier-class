import numpy as np
import matplotlib.pyplot as plt
import random as rd

class Glacier_block:
    def __init__(self, albedo = 0.5, density = 2008, heat_capacity = 0.0053 * 4.184 * 100/10**(2), latent_fusion = 333.55*1000, position = (), temp = 265, debris = (0, 0), snow = (0, 0), mass = 1, mass_change = (), neighbors = (), clock = 0):
        self.albedo = albedo #0.5, https://nsidc.org/cryosphere/seaice/processes/albedo.html
        self.density = density #2008kg/m^3
        self.heat_capacity = heat_capacity #0.0053 * 4.184 * 100/10**(2)
        self.latent_fusion = latent_fusion #333.55*1000

        self.position = list(position) #x- (elevation), y-coordinate
        self.temp = temp 
        self.debris = list(debris) #first value: presence/absence of debris, second value: depth
        self.snow = list(snow) #first value: presence/absence of snow, second value: depth
        self.mass = mass #unit mass
        self.mass_change = list(mass_change)
        self.neighbors = list(neighbors) #[left, right, front, behind]
        self.clock = clock #time tracker for snow to turn into ice

    def get_pos(self):
        return self.position

    def get_ele(self):
        return self.position[0]

    def get_temp(self):
        return self.temp

    def get_mass(self):
        return self.mass

    def get_mass_change(self):
        return self.mass_change[-1]

    def get_cover_state(self):
        return 'Debris depth = ' + str(self.debris[1]) + ' while snow depth = ' + str(self.snow[-1])
    
    def set_temp(self, temp): #input integer
        self.temp = temp
        
    def set_location(self, pos): #input list of 2
        self.position = pos

    def update_neighbors(self, l, r, f, b):
        self.neighbors = l, r, f, b
        
    def cover(self, cover_type, snow_depth, debris_depth):
        if cover_type == 'd':
            self.debris = [1, debris_depth]
            self.snow = [0, 0]
        else:
            self.debris = [0, 0]
            self.snow = [1, snow_depth]

    def mass_loss_dueto_flow(self): #loss 10% of its mass in each cycle
        self.mass -= 0.1 * self.mass #loss to the block one elevation unit below 
        self.mass += 0.1 * neighbors[-1].mass #get mass from the block one elevation unit above                    

    #ablation
    def change_temp(self, solar_constant, thermal_conductivity_of_rock, debris_temp_btm, thermal_conductivity_of_snow, snow_temp_btm):
        if self.debris[0] == 0 and self.snow[0] == 0:
            change_ice_temp = ((1 - self.albedo) * solar_constant) / (self.density * self.heat_capacity)
        elif self.debris[0] == 1:
            change_ice_temp = thermal_conductivity_of_rock * (debris_temp_btm[-1] - self.temp)/(self.density * self.heat_capacity)
        elif self.snow[0] == 1:
            change_ice_temp = thermal_conductivity_of_snow * (snow_temp_btm[-1] - self.temp)/(self.density * self.heat_capacity)
        self.temp = self.temp + change_ice_temp
        
    def melt(self, solar_constant, thermal_conductivity_of_rock, debris_temp_btm, snow_temp_btm): #if ice_temp > 273
        if self.debris[0] == 0 and self.snow[0] == 0:
            ice_mass_loss = -((1 - self.albedo) * solar_constant)/self.latent_fusion
        elif self.debris[0] == 1:
            ice_mass_loss = - thermal_conductivity_of_rock * (debris_temp_btm[-1] - self.temp)/self.latent_fusion
        elif self.snow[0] == 1:
            ice_mass_loss = - thermal_conductivity_of_snow * (snow_temp_btm[-1] - self.temp)/self.latent_fusion
        self.mass_change.append(ice_mass_loss)
        self.mass += ice_mass_loss
        
    #accumulation
    def accumulate(self, snow_input):
        self.snow[0] = 1
        self.snow[1] += snow_input

    def update_track(self):
        self.clock += 1

    def check_snow_turn_ice(self):
        if self.clock >= 3:
            ice_from_snow = self.snow[1] / 2.33 #density of ice = 2.33 * density of snow
            self.mass += ice_from_snow 
            self.mass_change.append(ice_from_snow)
            self.snow = [0, 0]
            print('Snow turned into ice!')
        else:
            print('Snow is still snow')
        
class Glacier:

    def __init__(self, ela = 0, temp_canvas = (), canvas = (), max_ele = 0):
        self.ela = ela
        self.temp_canvas = list(temp_canvas)
        self.canvas = list(canvas)
        self.max_elevation = max_ele #assume start from 0

    def show_glacier(self):
        return self.canvas

    def show_block_temp(self):
        temp = []
        for r in self.canvas:
            for c in r:
                temp.append(c.get_temp())
        t = np.array(temp)
        t.reshape(len(self.canvas), len(self.canvas[0]))
        return t
    
    def show_block_mass(self):
        masses = []
        for r in self.canvas:
            for c in r:
                masses.append(c.get_mass())
        m = np.array(masses)
        m.reshape(len(self.canvas), len(self.canvas[0]))
        return m

    def build_glacier(self, height, width, ice_temp, initial_snow_depth, initial_debris_depth):
        for i in range(height * width):
            self.temp_canvas.append(Glacier_block())
        self.temp_canvas = np.array(self.temp_canvas)
        self.canvas = self.temp_canvas.reshape(height, width)

        #update temp of all blocks
        for blocks in self.canvas:
            for b in blocks:
                b.set_temp(ice_temp)
                
        #update_block_location
        for i in range(0, len(self.canvas)):
            for j in range (0, len(self.canvas[0])):
                self.canvas[i][j].set_location([i, j])

        #update_cover_type
        self.max_elevation = len(self.canvas) #update max elevation
        self.ela = round(len(self.canvas) / 2) #assume ela at half of the elevation
        for elevation in self.canvas:
            for blocks in elevation:
                if blocks.get_ele() >= self.ela:
                    blocks.cover('s', initial_snow_depth, initial_debris_depth)
                else:
                    blocks.cover('d', initial_snow_depth, initial_debris_depth)
    
    def update_ice_temp(self, solar_constant, thermal_conductivity_of_rock, debris_temp_btm, thermal_conductivity_of_snow, snow_temp_btm):
        for r in self.canvas:
            for c in r:
                if c.get_temp() <= 273:
                    c.change_temp(solar_constant, thermal_conductivity_of_rock, debris_temp_btm, thermal_conductivity_of_snow, snow_temp_btm)
                else:
                    c.melt(solar_constant, thermal_conductivity_of_rock, debris_temp_btm, snow_temp_btm)

    #for plotting
    def avg_temp(self):
        temp = []
        for r in self.canvas:
            for c in r:
                temp.append(c.get_temp())
        return np.average(temp)

    def total_mass(self):
        mass = 0
        for r in self.canvas:
            for c in r:
                mass += c.get_mass()
        return mass        

#Initial Conditions
debris_temp_surface = [275]
debris_temp_btm = [268]
snow_temp_surface = [273]
snow_temp_btm = [268]
ice_temp = 265
night_temp_min = 265
night_temp_max = 289
glacier_dim = [5, 5] #[height, width]

#Constants
#Heatcap are reduced by a factor of 100 for now to make calculations faster
density_rock = 2243 #Based on Marl
heatcap_rock = 2000/10**(2)#Average value
density_snow = 200 #Average between new snow and wind packed snow (https://www.sciencelearn.org.nz/resources/1391-snow-and-ice-density)
heatcap_snow = 2000/10**(2) #dummy value

#Value obtained from https://pubs.usgs.gov/of/1988/0441/report.pdf
thermal_conductivity_of_rock = 2.390 * 10**(-3) * 4.184 * 100  #Conversion of cal/cm to joules/m
heat_transfer_coefficient_air = 500 #Between 0.5 to 1000 Wm^(-2)K^(-1)
heat_transfer_coefficient_water = 1450 #Between 50 to 3,000 Wm^(-2)K^(-1)
thermal_conductivity_of_snow = 0.1 #approximation based on https://www.climate-policy-watcher.org/snow/thermal-properties.html
heat_transfer_coefficient_water = 500 #to 10,000 Wm^(-2)K

soil_albedo = 0.17
snow_albedo = 0.90 #assume pure snow, can easily drop to 0.5 as it ages
solar_constant = 1400

#Lists
debris_depth = list(range(0, 10))
snow_depth = list(range(0, 10))
R_rock = debris_depth[-1]/thermal_conductivity_of_rock
R_snow = snow_depth[-1]/thermal_conductivity_of_rock
surface_temp_gradient_debris = surface_temp_gradient_snow = []
bottom_temp_gradient_debris = bottom_temp_gradient_snow = []
avg_temp_gradient_debris = avg_temp_gradient_snow = []
avg_temp_debris = avg_temp_snow = [] 
air_temp = []
glacier_avg_temp = []

#Initialize
g = Glacier()
g.build_glacier(glacier_dim[0], glacier_dim[1], ice_temp, max(snow_depth), max(debris_depth))
glacier_mass = [0, g.total_mass()]
print('Created glacier!')

#Use glacier total mass change to determine continue or stop
while glacier_mass[-1] > 0 and ((glacier_mass[-1] - glacier_mass[-2]) >= 1 or (glacier_mass[-1] - glacier_mass[0]) < 300):
#Step 1: Update debris and snow surface temp
    if len(air_temp)%2 == 0: #day
            air_temp.append(273)
            change_debris_surface_temp = (1 - soil_albedo) * solar_constant/(density_rock * heatcap_rock)
            change_snow_surface_temp = (1 - snow_albedo) * solar_constant/(density_snow * heatcap_snow)
            print('day ' + str(change_debris_surface_temp))
    else: #night
        air_temp.append(rd.randint(night_temp_min, night_temp_max))
        #Convection between cold air and warmer top of the debris layer
        change_debris_surface_temp = heat_transfer_coefficient_air * (air_temp[-1] - debris_temp_surface[-1])
        change_snow_surface_temp = heat_transfer_coefficient_air * (air_temp[-1] - snow_temp_surface[-1])
        print('night ' + str(change_debris_surface_temp))
            
        debris_temp_surface.append(debris_temp_surface[-1] + change_debris_surface_temp)
        snow_temp_surface.append(snow_temp_surface[-1] + change_snow_surface_temp)
        
        #Step 2: The radiation also heats up the bottom of the debris, but at a factor reduced by the depth
        debris_temp_btm.append(debris_temp_btm[-1] + change_debris_surface_temp/max(debris_depth))
        snow_temp_btm.append(snow_temp_btm[-1] + change_snow_surface_temp/max(snow_depth))
        print('Finished step 2!')
        #Step 3: Two thermal gradients form. One from the updated surface temperature to the bottom, and another from the botton (which is at ice temperature) to the surface. Take the average as the actual gradient.
    for i in range(len(debris_depth)):
        surface_temp_gradient_debris.append(-thermal_conductivity_of_rock * debris_depth[i] + debris_temp_surface[-1])
        bottom_temp_gradient_debris.append(-thermal_conductivity_of_rock * debris_depth[i] + debris_temp_btm[-1])

    for i in range(len(snow_depth)):
        surface_temp_gradient_snow.append(-thermal_conductivity_of_snow * debris_depth[i] + snow_temp_surface[-1])
        bottom_temp_gradient_snow.append(-thermal_conductivity_of_snow * debris_depth[i] + snow_temp_btm[-1])

    for i in range(len(surface_temp_gradient_debris)):
        avg_temp_gradient_debris.append((surface_temp_gradient_debris[i] + bottom_temp_gradient_debris[i])/2)

    for i in range(len(surface_temp_gradient_snow)):
        avg_temp_gradient_snow.append((surface_temp_gradient_snow[i] + bottom_temp_gradient_snow[i])/2)

    avg_temp_debris.append(np.average(avg_temp_gradient_debris))
    avg_temp_snow.append(np.average(avg_temp_gradient_snow))
    print('Finished step 3!')
    #Step 4: Ice temperature increases
    g.update_ice_temp(solar_constant, thermal_conductivity_of_rock, debris_temp_btm, thermal_conductivity_of_snow, snow_temp_btm)

    #Step 5: Update glacier total mass and average temp
    glacier_mass.append(g.total_mass())
    glacier_avg_temp.append(g.avg_temp())

fig, axs = plt.subplots(2,2)

axs[0,0].plot(range(len(glacier_mass)),glacier_mass)
axs[0,0].set_title("Total glacier mass per iteration")

axs[0,1].plot(range(len(glacier_avg_temp)), glacier_avg_temp)
axs[0,1].set_title("Change in ice temperature until iterations break")

#axs[1,0].plot(range(len(avg_temp)), avg_temp)
#axs[1,0].set_title("Change in average debris layer temperature per iteration")

#axs[1,1].plot(range(len(surface_temp_gradient)), surface_temp_gradient)
plt.show()
