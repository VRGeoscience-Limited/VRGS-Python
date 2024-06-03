#----------------------------------------------------------------------------------
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
#
# The MIT License (MIT)
#
# Copyright (c) 2024 David and Brian
#
# Permission is hereby granted, free of charge, to any person obtaining 
# a copy of this software and associated documentation files (the "Software"), 
# to deal in the Software without restriction, including without limitation the 
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
# sell copies of the Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in 
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
# TORT OR OTHERWISE, ARISING FROM, # OUT OF OR IN CONNECTION WITH THE SOFTWARE 
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#----------------------------------------------------------------------------------

# import appropriate packages
import vrgs
import numpy as np
import pandas as pd
import mplstereonet as mpl
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk

# set up Tkinter window
window = tk.Tk()
window.title("Stereonets in VRGS")

# get screen width and height and calulate x and y for gui window
w,h = 500,300
ws = window.winfo_screenwidth() # width of screen
hs = window.winfo_screenheight() # height of screen
x = (ws/2) - (w/2)
y = (hs/2) - (h/2)
# put window in centre of screen
window.geometry('%dx%d+%d+%d' % (w, h, x, y))
window.resizable(0,0)
window.configure(bg="#d5ddf3")

# label choices
def choice_selected(event):
    data_label = tk.Label(window, bg="#d5ddf3", text=combo_drop_menu1.get())
    data_label.option_clear
    data_label.place(x=45,y=100)
    plot_label = tk.Label(window, bg="#d5ddf3", text=combo_drop_menu2.get())
    plot_label.place(x=250,y=100)

# dropdown menu options
plot_options = ("Stereonet (with KMeans)", "Rose Diagram", "Northern Hemisphere Rose Diagram")
data_options = ("Joint", "Bedding")

# setting the string variables
plot_chosen = tk.StringVar(value = plot_options)
plot_chosen.set([0]) # set the initla plot style string
data_chosen = tk.StringVar(value = data_options)
data_chosen.set([0]) # set the initial data string

# create the dropdown menus
choice = tk.LabelFrame(text="Choose the data to analyse", bg = '#d5ddf3')
choice.grid(row=0, column=0, padx=10)
combo_drop_menu1 = ttk.Combobox(choice, textvariable=data_chosen, 
                                values=data_options, width=25, state='readonly')
combo_drop_menu1.set("Select the data to analyse")
#combo_drop_menu1.bind("<<ComboboxSelected>>", choice_selected)
combo_drop_menu1.grid(row=0, column=0, padx=10, pady=10)

choice2 = tk.LabelFrame(text="Choose the analyis and plot style", bg = '#d5ddf3')
choice2.grid(row=0, column=2, padx=10, pady=10)
combo_drop_menu2 = ttk.Combobox(choice2, textvariable=plot_chosen, 
                                values=plot_options, width=35, state='readonly')
combo_drop_menu2.set("Choose the analyis and plot style")
#combo_drop_menu2.bind("<<ComboboxSelected>>", choice_selected)
combo_drop_menu2.grid(row=0, column=1, padx=10, pady=10)

# create progress bar
progress_bar = ttk.Progressbar(window, orient='horizontal',
                               length = 400, mode='determinate')
progress_bar.place(x=50,rely=0.75)

# create button to interact with dropdown menu
calc_button = tk.Button(text='Analyse and Plot', bg='#F7F9FE', 
                        height=2, width=20, relief="flat", 
                        default='active', command = lambda:data_selected()) 
calc_button.place(relx=0.5, rely=0.5, anchor='center')

# function to simulate progress
def update_progress(step, total_steps):
    progress_bar['value'] = (step / total_steps) * 100
    window.update_idletasks()

# function calls based on user selection
def data_selected():
    """
    Method to call other data processing and plotting methods based on user input
    
    """
    # retrieve orientations from interpreations menu tree=
    orient_list = vrgs.orientations()
    df = pd.DataFrame(orient_list)
    df.columns =['name','group', 'id', 'x', 'y', 'z', 'dip', 'azimuth']
    
    #progress bar steps
    total_steps = 5  # Total steps for the process
    update_progress(1, total_steps)

    if combo_drop_menu1.get() == "Joint" and combo_drop_menu2.get() == "Stereonet (with KMeans)":
        print("You've chosen Joint and Stereonet (with KMeans)")
        strike, dip, centres = orientation_analysis(df, joint=True) # set joint to True to analyse joint data
        update_progress(2, total_steps)  # update progress bar
        stereonet_plot(strike, dip, centres, joint=True, kmeans=True) # add the returned array components from orientation_analysis() as function args
        

    elif combo_drop_menu1.get() == "Joint" and combo_drop_menu2.get() == "Rose Diagram":
        print("You've chosen Joint and Rose Diagram")
        strike, dip, centres = orientation_analysis(df, joint=True) # set joint to True to analyse joint data
        update_progress(2, total_steps)  
        rose_plot(strike, joint=True)
    
    elif combo_drop_menu1.get() == "Bedding" and combo_drop_menu2.get() == "Stereonet (with KMeans)":
        print("You've chosen Bedding and Stereonet (with KMeans)")
        strike, dip, centres = orientation_analysis(df, bedding=True) # set joint to True to analyse joint data
        update_progress(2, total_steps)  # update progress bar
        stereonet_plot(strike, dip, centres, bedding=True, kmeans=True) # add the returned array components from orientation_analysis() as function args
        
    elif combo_drop_menu1.get() == "Bedding" and combo_drop_menu2.get() == "Rose Diagram":
        print("You've chosen Joint and Rose Diagram")
        strike, dip, centres = orientation_analysis(df, bedding=True) # set joint to True to analyse joint data
        update_progress(2, total_steps)  # update progress bar
        rose_plot(strike, bedding=True)
        
    elif combo_drop_menu1.get() == "Joint" and combo_drop_menu2.get() == "Northern Hemisphere Rose Diagram":
        print("You've chosen Joint and Northern Hemisphere Rose Diagram")
        strike, dip, centres = orientation_analysis(df, joint=True) # set joint to True to analyse joint data
        update_progress(2, total_steps)  # update progress bar
        rose_plot(strike, joint=True, hemisphere=True) 
        

    elif combo_drop_menu1.get() == "Bedding" and combo_drop_menu2.get() == "Northern Hemisphere Rose Diagram":
        print("You've chosen Joint and Northern Hemisphere Rose Diagram")
        strike, dip, centres = orientation_analysis(df, bedding=True) # set joint to True to analyse joint data
        update_progress(2, total_steps)  # update progress bar
        rose_plot(strike, bedding=True, hemisphere=True)
        
    update_progress(5, total_steps)  # Process Completed

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


    update_progress(4, 5)  # Plot Generated
    update_progress(5, 5)  # Process Completed
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
    

    update_progress(4, 5)  # Plot Generated
    update_progress(5, 5)  # Process Completed
    plt.show()

window.mainloop()