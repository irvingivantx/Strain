import matplotlib.pyplot as plt
import numpy as np
import csv

from scipy.interpolate import CubicSpline, splrep, BSpline,splev

#READING DATA
Data=open('/Users/Irving/Programming/Python Projects/Strain/Dent/Axial/x_coordinates_dent_80_inches.csv')
x_data_axial=np.genfromtxt(Data,delimiter=",",dtype=float)

Data=open('/Users/Irving/Programming/Python Projects/Strain/Dent/Axial/y_coordinates_dent_80.csv')
y_data_axial=np.genfromtxt(Data,delimiter=",",dtype=float)

# Wall thickness
t = 0.281

# Pipe Radius
Ro=15

#CREATE ARRAY FOR MAX STRAIN VALUES.
Strain_iterations=10
step=-1/Strain_iterations
smoothing_factor=np.arange(0.1,0.01,step,dtype=float)
max_strain=np.array(range(len(smoothing_factor)),dtype=float)
for master_counter in range(len(smoothing_factor)):
    #CALCULATING SPLINE

    dent_profile_splinecoefficients=splrep(x_data_axial,y_data_axial,k=3,s=smoothing_factor[master_counter])
    dent_profile_spline=BSpline(*dent_profile_splinecoefficients,extrapolate=True)(x_data_axial)

    #CALCULATING DERVIATIVE OF SPLINE

    dent_profile_spline_deriv=splev(x_data_axial,dent_profile_splinecoefficients,der=1)
    dent_profile_spline_deriv2=splev(x_data_axial,dent_profile_splinecoefficients,der=2)

    #CREATING ARRAY OF RADIUS, STRAIN VALUES
    radius=np.array(range(len(x_data_axial)),dtype=float) # MAKE ARRAY TYPE FLOAT
    strain_long=np.array(range(len(x_data_axial)),dtype=float) # MAKE ARRAY


    # radius_circum=np.array(range(len(x_circum)),dtype=float) # MAKE ARRAY TYPE FLOAT FOR RADIUS - CIRCUMFERENTIAL
    # strain_circum=np.array(range(len(x_circum)),dtype=float) # MAKE ARRAY TYPE FLOAT FOR RADIUS - CIRCUMFERENTIAL
    for derivator in range(len(dent_profile_spline_deriv)):
        radius[derivator]=(1+(dent_profile_spline_deriv[derivator])**2)**(3/2)/abs(dent_profile_spline_deriv2[derivator])
        strain_long[derivator] = t / (2 * (radius[derivator]))
        # if dent_profile_spline_deriv[derivator-1]>0:
        #     if dent_profile_spline_deriv[derivator]<0:
        #         print("x -distance", x_data_axial[derivator])
        #         print("int", derivator)
        #         print("radius", radius[derivator])
        #         print("strain", strain_long[derivator])
        #         print("---------------")
        if dent_profile_spline_deriv2[derivator-1]>0:
            if dent_profile_spline_deriv2[derivator]<0:
                print("Inflection Point x -distance", x_data_axial[derivator])

    max_strain[master_counter]=np.amax(strain_long)


fig,ax=plt.subplots()
# fig, (ax, ax2)=plt.subplots(2,1)
#ax2.set_yscale('log')
# ax2.set_xlabel("Smoothing Factor")
# ax2.set_ylabel("Longitudinal Strain")
# #ax.scatter(smoothing_factor,max_strain,marker=".",color="blue")
# ax.plot(x_data_axial,dent_profile_spline,linestyle="-",color="red")
# ax2.plot(x_data_axial,dent_profile_spline_deriv,linestyle='-')
# ax2.plot(x_data_axial,dent_profile_spline_deriv2,linestyle='-')
# ax3.plot(x_data_axial,strain_long,linestyle='-')`dr4
ax.plot(x_data_axial,dent_profile_spline,linestyle="-",color="red")
ax.plot(x_data_axial,dent_profile_spline_deriv,linestyle='-')
ax.plot(x_data_axial,dent_profile_spline_deriv2,linestyle='-')
print(x_data_axial)
plt.show()