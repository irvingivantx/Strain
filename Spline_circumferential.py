import matplotlib.pyplot as plt
import numpy as np
import csv

from scipy.interpolate import CubicSpline, splrep, BSpline,splev


#READING DATA
Data=open('/Users/Irving/Programming/Python Projects/Strain/Dent/Circumferential/x_circumferential_coordinates_dent_80.csv')
x_data_circumferential=np.genfromtxt(Data,delimiter=",",dtype=float)

Data=open('/Users/Irving/Programming/Python Projects/Strain/Dent/Circumferential/y_circumferential_coordinates_dent_80.csv')
y_data_circumferential=np.genfromtxt(Data,delimiter=",",dtype=float)
# Wall thickness
t = 0.281

# Pipe Radius
Ro=15

#CREATE ARRAY FOR MAX STRAIN VALUES.
Strain_iterations=10000
step=-1/Strain_iterations
smoothing_factor=np.arange(1,.09,step,dtype=float)
max_strain=np.array(range(len(smoothing_factor)),dtype=float)
x_data_circumferential_transformed=np.arange(x_data_circumferential[0],x_data_circumferential[len(x_data_circumferential)-1],0.1)

for master_counter in range(len(smoothing_factor)):
    #CALCULATING SPLINE
    dent_profile_splinecoefficients=splrep(x_data_circumferential,y_data_circumferential,k=3,s=smoothing_factor[master_counter])
    dent_profile_spline=BSpline(*dent_profile_splinecoefficients,extrapolate=True)(x_data_circumferential_transformed)

    #CALCULATING DERVIATIVE OF SPLINE

    dent_profile_spline_deriv=splev(x_data_circumferential_transformed,dent_profile_splinecoefficients,der=1)
    dent_profile_spline_deriv2=splev(x_data_circumferential_transformed,dent_profile_splinecoefficients,der=2)

    #CREATING ARRAY OF RADIUS, STRAIN VALUES
    radius=np.array(range(len(x_data_circumferential_transformed)),dtype=float) # MAKE ARRAY TYPE FLOAT
    strain_long=np.array(range(len(x_data_circumferential_transformed)),dtype=float) # MAKE ARRAY



    # radius_circum=np.array(range(len(x_circum)),dtype=float) # MAKE ARRAY TYPE FLOAT FOR RADIUS - CIRCUMFERENTIAL
    # strain_circum=np.array(range(len(x_circum)),dtype=float) # MAKE ARRAY TYPE FLOAT FOR RADIUS - CIRCUMFERENTIAL
    for derivator in range(len(dent_profile_spline_deriv)):
        radius[derivator]=(1+(dent_profile_spline_deriv[derivator])**2)**(3/2)/abs(dent_profile_spline_deriv2[derivator])
        if dent_profile_spline_deriv2[derivator]<0:
            radiusofcurvature=-1*radius[derivator]
        else:
            radiusofcurvature=radius[derivator]
        strain_long[derivator] = (t/2)*((1/Ro)-(1/radiusofcurvature))


        # if dent_profile_spline_deriv[derivator-1]>0:
        #     if dent_profile_spline_deriv[derivator]<0:
        #         print("x -distance", x_data_circumferential_transformed[derivator])
        #         print("int", derivator)
        #         print("radius", radius[derivator])
        #         print("strain", strain_long[derivator])
        #         print("---------------")
    max_strain[master_counter]=np.max(strain_long)

fig, (ax,ax2,ax3)=plt.subplots(3,1)


#ax2.scatter(x_data_circumferential_transformed,dent_profile_spline_deriv,marker=".",color="blue")
ax.scatter(x_data_circumferential,y_data_circumferential,marker=".",color="blue")

ax.plot(x_data_circumferential_transformed,dent_profile_spline)
ax.plot(x_data_circumferential_transformed,dent_profile_spline_deriv)
ax2.scatter(smoothing_factor,max_strain,marker=".",color="blue")
ax3.plot(x_data_circumferential_transformed,strain_long,linestyle="-")
#LABELS
ax.set_xlabel("Channel Displacement [in]")
ax.set_ylabel("Deflections [in]")
ax2.set_xlabel("Smoothing Factor")
ax2.set_ylabel("Circumferential Strain")

plt.show()