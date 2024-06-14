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
#  These methods leverage the awesome `mplstereonet` package by 
#   Joe Kington - https://github.com/joferkington/mplstereonet
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

# Import appropriate libraries
import vrgs
import numpy as np
import pandas as pd
import mplstereonet as mpl
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk

# Retrieve orientations from interpreations menu tree=
orient_list = vrgs.orientations()
df = pd.DataFrame(orient_list)
df.columns =['name','group', 'id', 'x', 'y', 'z', 'dip', 'azimuth']

# Extract the last part of the 'group' column
df['group'] = df['group'].str.extract(r'Orientation/(.*?)/')[0]

def orientation_analysis(df, group):
    """
    Function to calculate strike from azimuth/dip data.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing orientation measurements.
    group : str
        Group to analyze (e.g., 'bedding', 'joint', 'fault').
    
    Returns
    -------
    tuple
        Arrays containing strike and dip values.
    """
    
    df2 = df[df['group'].str.contains(group, case=False)].copy()
    df2['strike'] = (df2['dip_azimuth'] - 90) % 360
    strike = df2['strike'].values
    dip = df2['dip'].values
    
    return strike, dip

def stereonet_plot_planes(strike, dip, update_progress=None, total_steps=None, current_step=None, progress_bar=None):
    """
    Plot stereonet of planes.

    Parameters
    ----------
    strike : array
        Array of strike values.
    dip : array
        Array of dip values.
    update_progress : function, optional
        Function to update the progress bar.
    total_steps : int, optional
        Total number of steps for progress bar.
    current_step : int, optional
        Current step for progress bar.
    progress_bar : ttk.Progressbar, optional
        Progress bar widget.
    """

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='equal_area_stereonet')
    ax.plane(strike, dip, c='k', linewidth=0.5, alpha=0.7)
    ax.grid()
    ax.set_title("Stereonet of Planes", y=1.1)

    if update_progress and progress_bar:
        update_progress(current_step + 1, total_steps, progress_bar)
    
    plt.show()

def stereonet_plot_contoured_poles(strike, dip, update_progress=None, total_steps=None, current_step=None, progress_bar=None):
    """
    Plot stereonet of contoured poles.

    Parameters
    ----------
    strike : array
        Array of strike values.
    dip : array
        Array of dip values.
    update_progress : function, optional
        Function to update the progress bar.
    total_steps : int, optional
        Total number of steps for progress bar.
    current_step : int, optional
        Current step for progress bar.
    progress_bar : ttk.Progressbar, optional
        Progress bar widget.

    References
    -----------
    [1] Vollmer, 1995. C Program for Automatic Contouring of Spherical
    Orientation Data Using a Modified Kamb Method. Computers &
    Geosciences, Vol. 21, No. 1, pp. 31--49.

    [2] Kamb, 1959. Ice Petrofabric Observations from Blue Glacier,
    Washington, in Relation to Theory and Experiment. Journal of
    Geophysical Research, Vol. 64, No. 11, pp. 1891-1909.
    """

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='equal_area_stereonet')
    ax.pole(strike, dip, c='k', markersize=4, alpha=0.7, label='Pole of the Planes')
    ax.density_contourf(strike, dip, measurement='poles', cmap='viridis', alpha=0.67)
    ax.density_contour(strike, dip, measurement='poles', colors='k', linewidths=0.5)
    ax.grid()
    ax.set_title("Stereonet of Contoured Poles", y=1.1)
    
    if update_progress and progress_bar:
        update_progress(current_step + 1, total_steps, progress_bar)
    
    plt.show()

def stereonet_plot_kmeans(strike, dip, num_clusters, update_progress=None, total_steps=None, current_step=None, progress_bar=None):
    """
    Plot stereonet of contoured poles.

    Parameters
    ----------
    strike : array
        Array of strike values.
    dip : array
        Array of dip values.
    num_clusters : int
        Number of clusters for KMeans.
    update_progress : function, optional
        Function to update the progress bar.
    total_steps : int, optional
        Total number of steps for progress bar.
    current_step : int, optional
        Current step for progress bar.
    progress_bar : ttk.Progressbar, optional
        Progress bar widget.

    References
    -----------
    [1] Vollmer, 1995. C Program for Automatic Contouring of Spherical
    Orientation Data Using a Modified Kamb Method. Computers &
    Geosciences, Vol. 21, No. 1, pp. 31--49.

    [2] Kamb, 1959. Ice Petrofabric Observations from Blue Glacier,
    Washington, in Relation to Theory and Experiment. Journal of
    Geophysical Research, Vol. 64, No. 11, pp. 1891-1909.
    """

    # Calcualte Kmeans centroids 
    centres = mpl.kmeans(strike, dip, num=num_clusters, measurement='poles')
    
    # Set up figures
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='equal_area_stereonet')
    ax.pole(strike, dip, c='k', markersize=4, alpha=0.7, label='Pole of the Planes')
    ax.density_contourf(strike, dip, measurement='poles', cmap='viridis', alpha=0.67)
    ax.density_contour(strike, dip, measurement='poles', colors='k', linewidths=0.5)
    ax.grid()

    # Plot centroids
    strike_cent, dip_cent = mpl.geographic2pole(*zip(*centres))
    ax.pole(strike_cent, dip_cent, 'ro', markeredgecolor='k', ms=7)

    # Label the centroids
    for (x0, y0) in centres:
        s, d = mpl.geographic2pole(x0, y0)
        x, y = mpl.pole(s, d)  # otherwise, we may get the antipode!

        kwargs = dict(xytext=(50, -25), ha='left') if x > 0 else dict(xytext=(-25, 30), ha='right')
        ax.annotate('{:03.0f}/{:03.0f}'.format(s[0], d[0]), xy=(x, y),
                    xycoords='data', textcoords='offset pixels',
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=.5'),
                    color='blue', **kwargs)

    ax.set_title(f"Stereonet of Contoured Poles with {num_clusters} KMeans Clusters", y=1.1)
    
    # Update the progress bar
    if update_progress and progress_bar:
        update_progress(current_step + 1, total_steps, progress_bar)
    
    plt.show()


def rose_plot(strike, hemisphere=False, bedding=False, joint=False, fault=False, bin_size=10, update_progress=None, total_steps=None, current_step=None, progress_bar=None):
    """
    Plot equal area rose diagram of strike values. 

    Parameters
    ----------
    strike : array
        Array of strike values.
    hemisphere : bool, optional
        Whether to plot only the northern hemisphere. Default is False.
    bedding : bool, optional
        Whether the plot is for bedding data. Default is False.
    joint : bool, optional
        Whether the plot is for joint data. Default is False.
    fault : bool, optional
        Whether the plot is for fault data. Default is False.
    bin_size : int, optional
        Size of bins for the histogram. Default is 10.
    update_progress : function, optional
        Function to update the progress bar.
    total_steps : int, optional
        Total number of steps for progress bar.
    current_step : int, optional
        Current step for progress bar.
    progress_bar : ttk.Progressbar, optional
        Progress bar widget.
    """
    # Calculate histogram for the rose plot
    circle_degrees = 360
    bin_edges = np.arange(-5, circle_degrees + bin_size, bin_size)
    number_of_strikes, bin_edges = np.histogram(strike, bin_edges)
    number_of_strikes[0] += number_of_strikes[-1]
    half = np.sum(np.split(number_of_strikes[:-1], 2), 0)
    two_halves = np.concatenate([half, half])

    # Normalize the heights to represent equal areas
    area = np.sum(two_halves)
    heights = np.sqrt(two_halves / area) * np.sqrt(circle_degrees / bin_size)

    # Set up the figure
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, polar=True)
    colors = plt.cm.winter_r(two_halves / bin_size)

    if hemisphere:
        ax.set_thetamin(-90)
        ax.set_thetamax(90)
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_xticks(np.arange(np.radians(-90), np.radians(91), np.radians(45)))
        ax.set_xticklabels(["270$\\degree$", '315$\\degree$', '0$\\degree$', '45$\\degree$', '90$\\degree$'], fontsize='medium')
        ax.grid(ls='--')
    else:
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.grid(ls='--')

    ax.bar(np.deg2rad(np.arange(0, circle_degrees, bin_size)), 
           heights, width=np.deg2rad(bin_size), 
           bottom=0.0, color=colors, edgecolor='k', alpha=0.8)

    # Plot titles based on plot chosen
    if bedding:
        ax.set_title("Rose Diagram for Bedding", y=1.1)
    elif joint:
        ax.set_title("Rose Diagram for Joint(s)", y=1.1)
    elif fault:
        ax.set_title("Rose Diagram for Fault(s)", y=1.1)

    if update_progress and progress_bar:
        update_progress(current_step + 1, total_steps, progress_bar)

    plt.show()

def setup_gui(df):
    """
    Setup and run the Tkinter GUI.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing orientation measurements.
    """

    def update_dropdown_menu(df, combo1, combo2):
        """
        Update the dropdown menu options based on the DataFrame.

        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame containing orientation measurements.
        combo1 : ttk.Combobox
            Combobox for data to analyze.
        combo2 : ttk.Combobox
            Combobox for analysis and plot style.
        """

        # Fill drop down menus 
        unique_groups = df['group'].unique().tolist()
        combo1['values'] = unique_groups
        combo2['values'] = plot_options
        set_combobox_width(combo2, plot_options)

    def set_combobox_width(combobox, values):
        """
        Set the width of the combobox to fit the longest item.

        Parameters
        ----------
        combobox : ttk.Combobox
            Combobox widget.
        values : list
            List of values in the combobox.
        """
        # Calculate max width for drop down menu based on string length
        max_width = max(len(str(value)) for value in values) + 1  # Add a little extra width for padding
        combobox.config(width=max_width)

    def data_selected():
        """
        Handle the event when the "Analyse and Plot" button is clicked.
        """
        # Add an error message
        error_label.config(text="", fg='red')  # Clear previous error message
        
        # Update to progress bar
        total_steps = 3
        analysis_type = combo_drop_menu1.get()
        plot_type = combo_drop_menu2.get()

        # Set analysis type based on selection
        if analysis_type == "Select the data to analyse" or plot_type == "Choose the analysis and plot style":
            error_label.config(text="Please choose options from both drop down menus.", fg='red')
            return
        
        if plot_type == "Stereonet of Contoured Poles with KMeans":
            if cluster_entry.get() == "":
                error_label.config(text="Please enter the number of clusters for KMeans.", fg='red')
                return
            try:
                num_clusters = int(cluster_entry.get())
            except ValueError:
                error_label.config(text="Please enter a valid integer for the number of clusters.", fg='red')
                return

        # Update progress bar
        update_progress(1, total_steps, progress_bar)

        # Conduct actual anlysis based on chosen type
        strike, dip = orientation_analysis(df, group=analysis_type)
        
        # Update progress bar again
        update_progress(2, total_steps, progress_bar)

        #  Create plots based on the chosen anlysis type
        if plot_type == "Stereonet of Planes":
            stereonet_plot_planes(strike, dip, update_progress=update_progress, total_steps=total_steps, current_step=2, progress_bar=progress_bar)
        elif plot_type == "Stereonet of Contoured Poles":
            stereonet_plot_contoured_poles(strike, dip, update_progress=update_progress, total_steps=total_steps, current_step=2, progress_bar=progress_bar)
        elif plot_type == "Stereonet of Contoured Poles with KMeans":
            stereonet_plot_kmeans(strike, dip, num_clusters, update_progress=update_progress, total_steps=total_steps, current_step=2, progress_bar=progress_bar)
        elif plot_type == "Rose Diagram":
            rose_plot(strike, update_progress=update_progress, total_steps=total_steps, current_step=2, progress_bar=progress_bar)
        elif plot_type == "Northern Hemisphere Rose Diagram":
            rose_plot(strike, hemisphere=True, update_progress=update_progress, total_steps=total_steps, current_step=2, progress_bar=progress_bar)
        else:
            error_label.config(text="Invalid combination of analysis type and plot type selected.", fg='red')
            return
        
        # Final update to progress bar
        update_progress(3, total_steps, progress_bar)

    def update_progress(step, total_steps, progress_bar):
        """
        Update the progress bar.

        Parameters
        ----------
        step : int
            Current step in the process.
        total_steps : int
            Total number of steps in the process.
        progress_bar : ttk.Progressbar
            Progress bar widget.
        """
        progress_bar['value'] = (step / total_steps) * 100
        window.update_idletasks()

    def on_plot_type_change(event):
        """
        Function to enable 'Number of Clusters' text box when 
        'Stereonet of Contoured Poles with KMeans' is selected.
        Disabled by default.

        Parameters
        ----------
        event : tkinter.Event
            The event object.
        """
        if combo_drop_menu2.get() == "Stereonet of Contoured Poles with KMeans":
            cluster_entry.config(state='normal')
        else:
            cluster_entry.config(state='disabled')
            cluster_entry.delete(0, tk.END)


    # Create GUI window and add components
    window = tk.Tk()
    window.title("Stereonets in VRGS")
    w, h = 700, 350
    ws = window.winfo_screenwidth()
    hs = window.winfo_screenheight()
    x = (ws / 2) - (w / 2)
    y = (hs / 2) - (h / 2)
    window.geometry('%dx%d+%d+%d' % (w, h, x, y))
    window.resizable(0, 0)
    window.configure(bg="#d5ddf3")

    # Set plot options for drop down menu
    plot_options = ["Stereonet of Planes", "Stereonet of Contoured Poles", "Stereonet of Contoured Poles with KMeans", "Rose Diagram", "Northern Hemisphere Rose Diagram"]

    # Data selection drop down menu
    choice = tk.LabelFrame(window, text="Select data to analyse", bg='#d5ddf3')
    choice.grid(row=0, column=0, padx=10, pady=10)
    combo_drop_menu1 = ttk.Combobox(choice, state='readonly')
    combo_drop_menu1.set("Select data to analyse")
    combo_drop_menu1.grid(row=0, column=0, padx=10, pady=10)

    # Anlasis style drop down menu
    choice2 = tk.LabelFrame(window, text="Choose the analysis and plot style", bg='#d5ddf3')
    choice2.grid(row=0, column=1, padx=10, pady=10)
    combo_drop_menu2 = ttk.Combobox(choice2, state='readonly')
    combo_drop_menu2.set("Choose the analysis and plot style")
    combo_drop_menu2.bind("<<ComboboxSelected>>", on_plot_type_change)
    combo_drop_menu2.grid(row=0, column=0, padx=10, pady=10)

    # Number of clusters text box
    choice3 = tk.LabelFrame(window, text="Number of Clusters for KMeans", bg='#d5ddf3')
    choice3.grid(row=0, column=2, padx=10, pady=10)
    cluster_entry = ttk.Entry(choice3, state='disabled')
    cluster_entry.grid(row=0, column=0, padx=10, pady=10)

    # Add progress bar
    progress_bar = ttk.Progressbar(window, orient='horizontal', length=int(w * 0.75), mode='determinate')
    progress_bar.place(x=int(w * 0.125), rely=0.65)

    # Error message
    error_label = tk.Label(window, text="", bg='#d5ddf3', fg='red')
    error_label.place(x=int(w * 0.125), y=int(h * 0.8), width=int(w * 0.75))

    # Add button to run analysis and plot the data
    calc_button = tk.Button(window, text='Analyse and Plot', bg='#F7F9FE', height=2, width=20, default='active', relief="flat", command=data_selected)
    calc_button.place(relx=0.5, rely=0.5, anchor='center')

    update_dropdown_menu(df, combo_drop_menu1, combo_drop_menu2)

    window.mainloop()

# Call the function to setup and display the GUI
setup_gui(df)
