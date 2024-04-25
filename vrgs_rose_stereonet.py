#-------------------------------------------------------------------------
#         Stereonets and Rose Diagrams from VRGS Orientation Data 
#               ------------------------------------------
#
#  Functions to work with VRGS orientation measurements:
#  1) Analyse _Orientation_ interpreations taken in VRGS
#  2) Create Equal Area stereonet plots. Planes and Poles.
#  3) Estimate KMeans of calculated poles.
#  4) Create rose diagram plots from strike orientation data
#       with the option to create single hemisphere plots
#
#   These methods leverage the awesome `mplstereonet` package by 
#    Joe Kington - https://github.com/joferkington/mplstereonet
#-------------------------------------------------------------------------

# import appropriate packages
import vrgs
import numpy as np
import pandas as pd
import mplstereonet as mpl
import matplotlib.pyplot as plt

def orientation_analysis(df, bedding=False, joint=False):
    """
    Function to calculate strike from azimuth/dip data,
    and calculate the KMeans centre (mode) of orientation data. 

    Parameters
    -----------
    df: Pandas Dataframe 
        Pandas Dataframe of orientation measurements.
    bedding: boolean
        Specify if plot is created for bedding orienation data.
    joint: boolean
        Specify if plot is created for joint orientation data.
        
    Returns
    -----------
    Arrays that contain strike, dip, and KMeans estimation 
    
    array[0] strike value in degrees
    array[1] dip value in degrees
    array[2] KMeans estimation of strike/dip
    """
    
    if joint == True and bedding == False:
        df2 = df[df['group'].str.contains("joint")==True].copy()
        df2['strike'] = df2['azimuth'] - 90
        df2.loc[df2.strike <0, 'strike'] += 360
        strike = df2['strike'].values
        dip = df2['dip'].values

        #KMeans estimation
        centres = mpl.kmeans(strike, dip, num=2, measurement='poles') # we pass 'poles' because we hare using sequence of strike and dip

        return strike, dip, centres
    
    elif bedding == True and joint == False:
        df1 = df[df['group'].str.contains("bedding")==True].copy()
        df1['strike'] = df1['azimuth'] - 90
        df1.loc[df1.strike <0, 'strike'] += 360
        strike = df1['strike'].values
        dip = df1['dip'].values

        # KMeans estimation
        centres = mpl.kmeans(strike, dip, num=1, measurement='poles') # we pass 'poles' because we hare using sequence of strike and dip

        return strike, dip, centres
    
    else:
        raise Exception("Please choose one orientation analysis method.")

def rose_plot(strike, hemisphere=False, bedding=False, joint=False, fault=False, bin_size=10):
    """
    Function to calculate 'rose' histogram of strike values.

    Parameters
    -----------
    strike: sequence of numbers
        The calculated strike in degrees. 
    bin_size: int or sequence of numbers (optional)
        The number of bins or a sequence of bin edges to use. Default is 10.
    hemisphere: boolean
        Specify if only northern hemisphere (half) polar plot is displayed.
    bedding: boolean
        Specify if plot is created for bedding orienation data.
    joint: boolean
        Specify if plot is created for joint orientation data.
    fault: boolean
        Specify if plot is created for fault orientation data.
        
    Returns
    -----------
    A matplotlib PatchCollection. 
    """
    
    #define constants
    circle_degrees = 360
    
    # calculate histogram of strike values
    bin_edges = np.arange(-5, circle_degrees + bin_size, bin_size)
    number_of_strikes, bin_edges = np.histogram(strike, bin_edges)
    
    # wrap the data around the plot
    number_of_strikes[0] += number_of_strikes[-1]
    
    # split and concatenate strike values
    half = np.sum(np.split(number_of_strikes[:-1], 2), 0)
    two_halves = np.concatenate([half, half])

    # create rose diagram plot
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111, polar='True')
    colors = plt.cm.winter_r(two_halves / bin_size)
    
    if hemisphere:
        ax.set_thetamin(-90)
        ax.set_thetamax(90)
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_xticks(np.arange(np.radians(-90),np.radians(91),np.radians(45))) 
        ax.set_xticklabels(["270$\\degree$",'315$\\degree$','0$\\degree$','45$\\degree$','90$\\degree$'],fontsize='medium')            # This argument controls azimuthal tick labels
        ax.grid(ls='--')
    else:
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.grid(ls='--')
    
    # create the histogram
    ax.bar(np.deg2rad(np.arange(0, circle_degrees, bin_size)), 
           two_halves, 
           width=np.deg2rad(bin_size), 
           bottom=0.0, color=colors, 
           edgecolor='k', alpha=0.8)

    if bedding:
        ax.set_title("Rose Diagram for Bedding", y=1.1)
    elif joint:
        ax.set_title("Rose Diagram for Joint(s)", y=1.1)
    elif fault:
        ax.set_title("Rose Diagram for Fault(s)", y=1.1)

    plt.show()

def stereonet_plot(strike, dip, centres, bedding=False, joint=False, fault=False, kmeans=False):
    """
    Function to create stereonets of bedding orientation measurements.
    Kamb contouring method with exponential smoothing is used(See Vollmer, 1995).
    Filled contours created by default.

    Parameters
    -----------
    strike: sequence of numbers
        The calculated strike in degrees.
    dip: sequence of numbers
        The observed dip values in degrees.
    bedding: boolean
        Specify if plot is created for bedding orienation data.
    joint: boolean
        Specify if plot is created for joint orientation data.
    fault: boolean
        Specify if plot is created for fault orientation data.
    kmeans: boolean
        Specify if KMeans plot is created for orientation data.

    Returns
    -----------
    A matplotlib PatchCollection

    References
    -----------
    [1] Vollmer, 1995. C Program for Automatic Contouring of Spherical
       Orientation Data Using a Modified Kamb Method. Computers &
       Geosciences, Vol. 21, No. 1, pp. 31--49.

    [2] Kamb, 1959. Ice Petrofabric Observations from Blue Glacier,
       Washington, in Relation to Theory and Experiment. Journal of
       Geophysical Research, Vol. 64, No. 11, pp. 1891-1909.
    """

    fig = plt.figure(figsize=(17,6))
    
    ax1 = fig.add_subplot(131, projection='equal_area_stereonet')
    ax1.plane(strike, dip, c='k', linewidth=0.5, alpha=0.7)
    ax1.grid()
    
    ax2 = fig.add_subplot(132, projection='equal_area_stereonet')
    ax2.pole(strike, dip, c='k', markersize=4, alpha=0.7, label='Pole of the Planes')
    ax2.density_contourf(strike, dip, measurement='poles', cmap='viridis', alpha=0.67)
    ax2.density_contour(strike, dip, measurement='poles', colors='k', linewidths=0.5)
    ax2.grid()
    
    if kmeans:
        # KMeans plot
        ax3 = fig.add_subplot(133, projection='stereonet')        
        ax3.pole(strike, dip, c='k', markersize=4, alpha=0.7, label='Pole of the Planes')
        ax3.density_contourf(strike, dip, measurement='poles', cmap='viridis', alpha=0.67)
        ax3.density_contour(strike, dip, measurement='poles', colors='k', linewidths=0.5)
        ax3.grid()
        ax3.set_title('KMeans of Poles', y=1.1)
        
        # find the two modes
        strike_cent, dip_cent = mpl.geographic2pole(*zip(*centres))
        ax3.pole(strike_cent, dip_cent, 'ro', markeredgecolor='k',ms=7)
        
        # label the modes
        for (x0, y0) in centres:
            s, d = mpl.geographic2pole(x0, y0)
            x, y = mpl.pole(s, d) # otherwise, we may get the antipode!
        
            if x > 0:
                kwargs = dict(xytext=(50, -25), ha='left')
            else:
                kwargs = dict(xytext=(-25, 30), ha='right')
        
            ax3.annotate('{:03.0f}/{:03.0f}'.format(s[0], d[0]), xy=(x, y),
                        xycoords='data', textcoords='offset pixels',
                        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=.5'),
                        color='blue',**kwargs)
            
    if bedding:
        ax1.set_title("Planes of Bedding Measurements", y=1.1)
        ax2.set_title("Density Contour of Bedding Poles", y=1.1)
    elif joint:
        ax1.set_title("Planes of Joint Measurements", y=1.1)
        ax2.set_title("Density Contour of Joint Poles", y=1.1)
    elif fault:
        ax1.set_title("Planes of Fault Measurements", y=1.1)
        ax2.set_title("Density Contour of Fault Poles", y=1.1)
    
    plt.show()

# process data and call functions

# retrieve orientations from interpreations menu tree
orient_list = vrgs.orientations()
#print(orient_list)

# convert list to a pandas dataframe for ease of use
df = pd.DataFrame(orient_list)
df.columns =['name','group', 'id', 'x', 'y', 'z', 'dip', 'azimuth']
print(df)

# analyse the data to extract strike and dip measurements
#bed_strike, bed_dip, bed_centre = orientation_analysis(df, bedding=True) # set bedding to True to analyse bedding data 
#joint_strike, joint_dip, joint_centre = orientation_analysis(df, joint=True) # set joint to True to analyse joint data

# plot stereonets of bedding orientation data
#stereonet_plot(bed_strike, bed_dip, bed_centre, bedding=True) # add the returned array components from orientation_analysis() as function args

#plot stereonets of joint orientation data
#stereonet_plot(joint_strike, joint_dip, joint_centre, joint=True, kmeans=True) # add orientation_analysis() arrays as function args and set kmeans to True

# rose diagram of the joint measurements
#rose_plot(joint_strike, joint=True)

# Use hemisphere=True to plot only "northern" hemisphere polar plot
#rose_plot(joint_strike, joint=True, hemisphere=True)