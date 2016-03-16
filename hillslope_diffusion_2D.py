
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 20:16:47 2016

@author: jackybaughman

This was an attempt at 2D Hillslope diffusion...but it definitely doesn't show
anything as is...
The newest to python...would love any suggestions!!
"""

from landlab import RasterModelGrid
import pylab, time
import numpy as np

def main():
      
        # INITIALIZE

        # User-defined parameter values
        rows = 20          # number of rows in the grid
        cols = 30        # number of columns in the grid
        dx = 10.0             # grid cell spacing
        k = 0.01             # diffusivity coefficient, in m2/yr
        num_time_steps = 10000 # number of time steps in run
        Rho_r = 2700 # rock density kg/m3
        Rho_s = 1600 # soil density kg/m3
        kappa=0.003 # m^2/yr 
        k=kappa*Rho_s # efficientcy 
        max_elev = 1000 # m
        wdotnaught = 0.005 # m/yr
        Hstar = 0.4 # m
        hill_slope = .5
        
        # Derived parameters
        dt = 0.1*dx**2 / k    # time-step size set by CFL condition

        # Create and initialize a raster model grid
        mg = RasterModelGrid(rows, cols, dx)

        # Set the boundary conditions
        mg.set_closed_boundaries_at_grid_edges(False, False, True, True)

        # Set up scalar values
        zb = mg.add_zeros('node', 'Elevation')            
        H = mg.add_zeros('node', 'regolith_thickness')
        z_H = mg.add_empty('node', 'regolith_elevation')
        dzda = mg.add_zeros('link', 'surface_slope')


        zb[:] = max_elev - mg.node_x * hill_slope # bedrock height
        z_H[:] = zb + H # bedrock + regolith
        
        # Get a list of the core cells
        core_nodes = mg.core_nodes

        # Display a message, and record the current clock time
        print( 'Running diffusion_with_model_grid.py' )
        print( 'Time-step size has been set to ' + str( dt ) + ' years.' )
        start_time = time.time()

        # RUN

        # Main loop
        for i in range(0, num_time_steps):

                # Calculate the weathering rate
                wdot = wdotnaught * np.exp(-H/Hstar)
                
                #slope
                dzda = mg.calculate_gradients_at_active_links(z_H)
                
                # Flux
                Q = -k*dzda

                # flux gradient
                dqds = mg.calculate_flux_divergence_at_nodes(Q)

                # rate of regolith thickness change
                dHdt = ((Rho_r/Rho_s)*wdot) - ((1/Rho_s)* dqds)
                
                # regolith thickness
                H = H + (dHdt*dt)
                
                # bedrock
                zb = zb-(wdot*dt)
                
                #Don't know syntax for placing channel boundary conditions...
#zb(2:end-1) = zb(2:end-1)-(wdot(2:end-1)*dt); % change of bedrock due to weathering
#zb(1)=zb(1)-(edot*dt); % fixed boundary condition channel erosion
#zb(end)=zb(end)-(edot*dt); % fixed boundary condition channel erosion
                
                
                # Update the elevations
                z_H[core_nodes] = zb[core_nodes] + H[core_nodes]
                

        # FINALIZE

        # Get a 2D array version of the elevations
        zr = mg.node_vector_to_raster(H)

        # Create a shaded image
        pylab.close()  # clear any pre-existing plot
        im = pylab.imshow(zr, cmap=pylab.cm.RdBu, extent=[0,cols*dx,0,rows*dx],
                                          origin='lower')
        # add contour lines with labels
        cset = pylab.contour(zr, extent=[0,cols*dx,rows*dx,0], hold='on',
                                                 origin='image')
        pylab.clabel(cset, inline=True, fmt='%1.1f', fontsize=10)

        # add a color bar on the side
        cb = pylab.colorbar(im)
        cb.set_label('Thickness in meters')

        # add a title and axis labels
        pylab.title('Regolith thickness')
        pylab.xlabel('Distance (m)')
        pylab.ylabel('Distance (m)')

        # Display the plot
        pylab.show()
        print('Run time = '+str(time.time()-start_time)+' seconds')
        
                # Get a 2D array version of the elevations
        zr = mg.node_vector_to_raster(z_H)

        # Create a shaded image
        pylab.close()  # clear any pre-existing plot
        im = pylab.imshow(zr, cmap=pylab.cm.RdBu, extent=[0,cols*dx,0,rows*dx],
                                          origin='lower')
        # add contour lines with labels
        cset = pylab.contour(zr, extent=[0,cols*dx,rows*dx,0], hold='on',
                                                 origin='image')
        pylab.clabel(cset, inline=True, fmt='%1.1f', fontsize=10)

        # add a color bar on the side
        cb = pylab.colorbar(im)
        cb.set_label('Elevation in meters')

        # add a title and axis labels
        pylab.title('Simulated topography')
        pylab.xlabel('Distance (m)')
        pylab.ylabel('Distance (m)')

        # Display the plot
        pylab.show()
        print('Run time = '+str(time.time()-start_time)+' seconds')

if __name__ == "__main__":
        main()