#----------------------------------------------------------------------------------
#              Fracture / Fault Slip Susceptibility Analysis 
#               -------------------------------------------
#
# Equations & code from:
# Morris et al., 1996. Geology. Slip-tendency analysis and fault reactivation
# Ferrill et al., 1999. GSA Today. Stressed Rock Strains Groundwater at Yucca Mountain, Nevada
# Allmendinger et al., 2012 Structural Geology Algorithms, Cambridge, University Press
# Stephens et al., 2018. EPSL. Mechanical models to estimate the paleostress state from igneous intrusions
#       
# Method implementation inspired 
# and modified from FracTend - https://github.com/DaveHealy-github/FracTend
#
# The MIT License (MIT)
#
# Copyright (c) 2024 Brian and David
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
import tkinter as tk
from tkinter import ttk
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize 
import mplstereonet

class StressAnalysis:
    """Class to conduct stress analyses on geological faults and fractures."""

    def __init__(self, dips, dip_azimuths, 
                 sigma1, sigma2, sigma3, 
                 trend1, plunge1, 
                 trend2, plunge2, 
                 trend3, plunge3):
        """
        Initialise the StressAnalysis class with the provided parameters.

        Parameters:
        dips (array-like): Dip angles of the fractures.
        dip_azimuths (array-like): Dip azimuth angles of the fractures.
        sigma1 (float): Magnitude of the maximum principal stress.
        sigma2 (float): Magnitude of the intermediate principal stress.
        sigma3 (float): Magnitude of the minimum principal stress.
        trend1 (float): Trend of sigma1.
        plunge1 (float): Plunge of sigma1.
        trend2 (float): Trend of sigma2.
        plunge2 (float): Plunge of sigma2.
        trend3 (float): Trend of sigma3.
        plunge3 (float): Plunge of sigma3.
        """
        self.dips = dips
        self.dip_azimuths = dip_azimuths
        self.sigma1 = sigma1
        self.sigma2 = sigma2
        self.sigma3 = sigma3
        self.trend1 = trend1
        self.plunge1 = plunge1
        self.trend2 = trend2
        self.plunge2 = plunge2
        self.trend3 = trend3
        self.plunge3 = plunge3
        self.sigma_rotated = self.rotate_stress_tensor()
        self.df_tendencies = self.calculate_tendencies()

    @staticmethod
    def deg_to_rad(deg):
        """
        Convert degrees to radians.

        Parameters:
        deg (float): Angle in degrees.

        Returns:
        float: Angle in radians.
        """
        return np.radians(deg)

    @staticmethod
    def orientation_vector(trend, plunge):
        """
        Compute the orientation vector from trend and plunge.

        The orientation vector are the cosines of the angles 
        between the line and the three coordinate axes (x, y, z). 

        Parameters:
        trend (float): The trend of the line in degrees.
        plunge (float): The plunge of the line in degrees.

        Returns:
        numpy.ndarray: A 3-element array containing the orientation vector.
        """

        # Convert trend and plunge from degrees to radians
        trend = StressAnalysis.deg_to_rad(trend) # azimuthal angle in the horizontal plane, measured clockwise from the north 
        plunge = StressAnalysis.deg_to_rad(plunge) # The angle between the line and the horizontal plane, measured downward
        
        # Compute the orientation vectors

        # along the x-axis (east direction)
        # project the line horizontally (onto a plane) then into the east direction
        east_line = np.cos(trend) * np.cos(plunge) 
        
        # along the y-axis (north direction)
        #project the line horizontally (onto a plane) then into the north direction
        north_line = np.sin(trend) * np.cos(plunge)  
        
        # along the z-axis (downward direction)
        # vertical projection of the line
        down_line = np.sin(plunge) 

        orientation_vectors = np.array([east_line, north_line, down_line])

        return orientation_vectors

    @staticmethod
    def dip_azimuth_to_strike_dip(dip, dip_azimuth):
        """
        Convert dip and dip azimuth to strike and dip.

        Parameters:
        dip (float): Dip angle in degrees.
        dip_azimuth (float): Dip azimuth angle in degrees.

        Returns:
        tuple: Strike and dip angles in degrees.
        """
        strike = (dip_azimuth - 90) % 360
        return strike, dip

    @staticmethod
    def trend_plunge_to_plane(trend, plunge):
        """
        Convert trend and plunge to strike and dip planes.

        Parameters:
        trnd (float): trend angle in degrees.
        plunge (float): plunge angle in degrees.

        Returns:
        tuple: Strike and dip angles in degrees.
        """
        strike = (trend + 90) % 360
        dip = 90 - plunge
        return strike, dip

    @staticmethod
    def fault_normal_vector(dip, dip_azimuth):
        """
        Compute the normal vector to the fault plane given dip and dip azimuth.

        Uses Cartesian coordinate system (z is positive upwards). 
        The normal vector in the z-direction is negative when the plane dips downward.

        Parameters:
        dip (float): Dip angle in degrees.
        dip_azimuth (float): Dip azimuth angle in degrees.

        Returns:
        numpy.ndarray: A 3-element array representing the normal vector to the fault plane.
        """
        # Convert dip and dip azimuth from degrees to radians
        dip_rad = StressAnalysis.deg_to_rad(dip)
        dip_azimuth_rad = StressAnalysis.deg_to_rad(dip_azimuth)
        
        n = [
            np.sin(dip_rad) * np.cos(dip_azimuth_rad), # along the x-axis (east direction)
            np.sin(dip_rad) * np.sin(dip_azimuth_rad), # along the y-axis (north direction)
            -np.cos(dip_rad)                           # along the z-axis (downward direction)
        ]

        n_arr = np.array(n)

        return n_arr

    def rotate_stress_tensor(self):
        """
        Rotate the stress tensor based on the trends and plunges of the principal stresses.

        This method transforms the stress tensor from the principal stress
        coordinates to the geographic coordinates using the trends and plunges of the
        principal stresses.

        This method rotates the stress tensor to align it with the geographic coordinate system 
        based on the given trends and plunges of the principal stresses (sigma1, sigma2, sigma3).

        1. Create a diagonal matrix from the principal stress magnitudes. 
        
        2. Calculate the orientation vector for the principal stresses using 
        their trends and plunges. The orientation vectors are used to form a rotation matrix. 
        
        3. The rotation matrix is applied to the stress tensor, transforming it 
        into geographic coordinates.

        Returns:
        numpy.ndarray: Rotated stress tensor in geographic coordinates.
        """
        # Creates a diagonal matrix (tensor) from the magnitudes of the principal stresses
        principal_stress_tensor = np.diag([self.sigma1, self.sigma2, self.sigma3])
        
        # Calculate the orientation vector for the maximum, intermediate, and minimum principal stresses.
        # Vectors represent the directions of the principal stresses in 3D space.
        orientation_vector_sigma1  = self.orientation_vector(self.trend1, self.plunge1)
        orientation_vector_sigma2  = self.orientation_vector(self.trend2, self.plunge2)
        orientation_vector_sigma3  = self.orientation_vector(self.trend3, self.plunge3)

        # Construct the rotation matrix by stacking the orientation vectors to form a 3x3 matrix.
        stacked_orientation_vectors  = np.vstack((orientation_vector_sigma1, orientation_vector_sigma2, orientation_vector_sigma3))         

        # Tranpose the stacked matrix so vectors are columns
        rotation_matrix = stacked_orientation_vectors.T
                
        # Rotate the stress tensor from the principal stress coordinate system to the geographic coordinate system.
        # sigma_rotated = rotation_matrix * intermediate tensor = (principal_stress_tensor * rotation_matrix^T)
        intermediate_tensor = np.dot(principal_stress_tensor, stacked_orientation_vectors)

        sigma_rotated = np.dot(rotation_matrix, intermediate_tensor)


        return sigma_rotated
    
    def convert_stress_axes_to_poles(self):
        """
        Convert principal stress axes to planes and then to poles for mplstereonet.

        This method transforms the principal stress axes (sigma1, sigma2, sigma3) 
        defined by their trends and plunges into plane representations. It then 
        converts these planes into their corresponding poles for stereonet plotting.

        Returns:
        tuple: Two lists containing the longitudes and latitudes of the poles 
            representing the principal stress axes on a stereonet.
        """
        
        # Convert trends and plunges of the principal stress axes into planes
        stress_poles = [
            self.trend_plunge_to_plane(self.trend1, self.plunge1), 
            self.trend_plunge_to_plane(self.trend2, self.plunge2),
            self.trend_plunge_to_plane(self.trend3, self.plunge3)
        ]
        
        # Convert the planes into poles (longitudes and latitudes) for stereonet plotting
        stress_lons, stress_lats = zip(*[mplstereonet.pole(strike, dip) for strike, dip in stress_poles])
        
        return stress_lons, stress_lats

    @staticmethod
    def resolve_stresses_on_fault(sigma_rotated, normal_vector):
        """
        Resolve the normal and shear stresses on the fault plane using the rotated stress tensor.

        This method calculates the normal and shear stresses acting on a fault plane using 
        the rotated stress tensor and the normal vector to the fault plane. 
        
        1. Compute the normal stress (sigma_n) by projecting the stress tensor 
        onto the normal vector. 
        
        2. Calculate the full stress vector acting on the fault plane. 
        The shear stress (tau) is determined by finding the component of this stress vector 
        that is parallel to the fault plane.

        Parameters:
        sigma_rotated (numpy.ndarray): Rotated stress tensor.
        normal_vector (numpy.ndarray): Normal vector to the fault plane.

        Returns:
        tuple: Normal stress and shear stress on the fault plane.
        """

        # Calculate the stress vector acting on the fault/fracture plane
        # This is done by multiplying the rotated stress tensor by the normal vector to the fault plane
        stress_on_plane = np.dot(sigma_rotated, normal_vector)
        
        # Calculate the normal stress on the fault/fracture plane
        # This is the projection of the stress tensor along the normal vector
        sigma_n = np.dot(normal_vector, stress_on_plane)
        
        # Calculate the magnitude of the shear stress on the fault/fracture plane
        # Subtract normal component from the stress vector to derive shear component
        # Compute magnitude of shear component using Euclidean norm 
        tau = np.linalg.norm(stress_on_plane - sigma_n * normal_vector)
        
        return sigma_n, tau

    @staticmethod
    def slip_tendency(normal_stress, shear_stress):
        """
        Calculate the slip tendency from Morris et al., 1996.

        Parameters:
        normal_stress (float): Normal stress on the fault plane.
        shear_stress (float): Shear stress on the fault plane.

        Returns:
        float: Slip tendency.
        """
        return shear_stress / normal_stress

    def calculate_tendencies(self):
        """
        Calculate various stress-related tendencies for each fracture orientation and combine the results into a single DataFrame.

        Returns:
        pandas.DataFrame: DataFrame containing the following columns:
                          - 'Dip': Dip angle of the fracture.
                          - 'Dip Azimuth': Dip azimuth angle of the fracture.
                          - 'Strike': Strike angle of the fracture.
                          - 'Normal Stress': Normal stress on the fault plane.
                          - 'Shear Stress': Shear stress on the fault plane.
                          - 'Ts': Slip tendency.
                          - 'Ts/Tsmax': Normalised slip tendency.
        """
        results = []
        for dip, dip_azimuth in zip(self.dips, self.dip_azimuths):
            
            # Convert dip and dip azimuth to strike and dip
            strike, _ = self.dip_azimuth_to_strike_dip(dip, dip_azimuth)
            
            # Calculate normal vector to fault/fracture
            n = self.fault_normal_vector(dip, dip_azimuth)
            
            # Calculate normal stress and shear stress
            sigma_n, tau = self.resolve_stresses_on_fault(self.sigma_rotated, n)
            
            # Calculate slip tendency 
            slip_tendency_value = self.slip_tendency(sigma_n, tau)
            
            results.append([dip, dip_azimuth, strike, sigma_n, tau, slip_tendency_value])

        df = pd.DataFrame(results, columns=['Dip', 'Dip Azimuth', 'Strike', 'Normal Stress', 'Shear Stress', 'Ts'])

        df['Ts/Tsmax'] = df['Ts'] / df['Ts'].max()
        return df

class StereonetPlotter:
    """
    Class for plotting Stereonets from tendency analyses.

    This class provides methods to generate and plot stereonets for various tendency analyses, 
    including slip tendency.
    """

    def __init__(self, stress_analysis):
        self.stress_analysis = stress_analysis

    def generate_dip_azimuth_grid(self, dip_range, dip_azimuth_range):
        """
        Generate grid of dip and dip azimuth values, and compute normalised slip tendency

        Parameters:
        dip_range (list): List of dip values.
        dip_azimuth_range (list): List of dip azimuth values.

        Returns:
        tuple: Lists of dips, dip azimuths, normalised slip tendencies
        """
        all_dips = []
        all_dip_azimuths = []
        all_slip_tendencies = []

        for dip in dip_range:
            for dip_azimuth in dip_azimuth_range:
                normal_vector = self.stress_analysis.fault_normal_vector(dip, dip_azimuth)
                normal_stress, shear_stress = self.stress_analysis.resolve_stresses_on_fault(self.stress_analysis.sigma_rotated, normal_vector)
                slip_tendency_value = self.stress_analysis.slip_tendency(normal_stress, shear_stress)

                all_dips.append(dip)
                all_dip_azimuths.append(dip_azimuth)
                all_slip_tendencies.append(slip_tendency_value)

        max_slip_tendency = max(all_slip_tendencies)
        normalised_slip_tendencies = [st / max_slip_tendency for st in all_slip_tendencies]
        return all_dips, all_dip_azimuths, normalised_slip_tendencies

    @staticmethod
    def plot_stereonet(ax, lons, lats, values, user_lons, user_lats, stress_lons, stress_lats, levels, cmap, title):
        """
        Plot stereonet with contours and data points.

        Parameters:
        ax (matplotlib.axes._subplots.AxesSubplot): The axis to plot on.
        lons (array): Longitudes of grid points.
        lats (array): Latitudes of grid points.
        values (array): Values at each grid point.
        user_lons (array): Longitudes of provided data points.
        user_lats (array): Latitudes of provided data points.
        stress_lons (array): Longitudes of stress axes.
        stress_lats (array): Latitudes of stress axes.
        levels (array): Contour levels.
        cmap (str): Colormap for the plot.
        title (str): Title of the plot.

        Returns:
        contour: The contour plot object.
        """
        contour = ax.tricontourf(lons, lats, values, levels=levels, cmap=cmap)
        ax.scatter(user_lons, user_lats, c='w', edgecolors='k', s=50, label='Orientation Data',zorder=3)
        ax.scatter(stress_lons[0], stress_lats[0], c='k', edgecolors='w', s=200, marker='s', label=r'$\sigma_1$', zorder=3, clip_on=False)
        ax.scatter(stress_lons[1], stress_lats[1], c='k', edgecolors='w', s=200, marker='d', label=r'$\sigma_2$', zorder=3, clip_on=False)
        ax.scatter(stress_lons[2], stress_lats[2], c='k', edgecolors='w', s=200, marker='^', label=r'$\sigma_3$', zorder=3, clip_on=False)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.04), ncol=4, handletextpad=0.3)
        ax.set_azimuth_ticks([])
        ax.set_azimuth_ticklabels([])
        ax.set_title(title, fontsize=14, y=1.05)

        return contour

    def add_colourbar(self, fig, contour, levels, label, cmap, fontsize=12):
        """
        Add a colorbar to the plot with specified font size.

        Parameters:
        fig (matplotlib.figure.Figure): The figure to add the colorbar to.
        contour (QuadContourSet): The contour set for which the colorbar is created.
        levels (array): Contour levels.
        label (str): Label for the colorbar.
        cmap (str): Colormap used in the plot.
        fontsize (int): Font size for the colorbar label.
        """
        norm = Normalize(vmin=min(levels), vmax=max(levels))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(contour, ax=fig.axes, orientation='vertical', pad=0.1, aspect=30, shrink=1)
        cbar.set_label(label, fontsize=fontsize)
        cbar.set_ticks(levels[::2])
        cbar.set_ticklabels([f"{val:.2f}" if max(levels) <= 1 else f"{int(val)}" for val in levels[::2]])

    def create_tendency_plot(self, dip_range, dip_azimuth_range, plot_type='slip'):
        """
        Create and display a tendency plot on a stereonet.

        Parameters:
        dip_range (list): List of dip values to generate grid.
        dip_azimuth_range (list): List of dip azimuth values to generate grid.
        plot_type (str): Type of plot to create ('slip').
        """
        
        df_tendencies = self.stress_analysis.df_tendencies
        user_dips = df_tendencies['Dip'].values
        user_dip_azimuths = df_tendencies['Dip Azimuth'].values

        all_dips, all_dip_azimuths, normalised_slip_tendencies= self.generate_dip_azimuth_grid(dip_range, dip_azimuth_range)

        all_strikes = [self.stress_analysis.dip_azimuth_to_strike_dip(dip, da)[0] for dip, da in zip(all_dips, all_dip_azimuths)]
        lons, lats = mplstereonet.pole(all_strikes, all_dips)

        user_strikes = [self.stress_analysis.dip_azimuth_to_strike_dip(dip, da)[0] for dip, da in zip(user_dips, user_dip_azimuths)]
        user_lons, user_lats = mplstereonet.pole(user_strikes, user_dips)

        stress_lons, stress_lats = self.stress_analysis.convert_stress_axes_to_poles()

        if plot_type == 'slip':
            values = normalised_slip_tendencies
            levels = np.linspace(0, 1, 21)
            cmap = 'viridis'
            title = 'Estimated Fracture Slip Tendency - $T_s$'
            colorbar_label = 'Normalised Slip Tendency ($T_s$)'

        fig, ax = plt.subplots(1, 1, figsize=(8, 6), subplot_kw=dict(projection='stereonet'))

        contour = self.plot_stereonet(ax, lons, lats, values, user_lons, user_lats, stress_lons, stress_lats, levels, cmap, title)
        self.add_colourbar(fig, contour, levels, colorbar_label, cmap)
    
    def show_all_plots(self):
        plt.show()

class StressGUI:
    """
    A GUI application for plotting Stereonets from stress tendency analyses."""

    def __init__(self, root, df):
        """
        Initialise the StressGUI class with the provided parameters.

        Parameters:
        root (tk.Tk): Root window of the GUI.
        df (pd.DataFrame): DataFrame containing the input data.
        """
        self.root = root
        self.df = df
        
        self.centre_window(725, 425)
        
        self.root.title("Stress and Tendency Analysis")
        self.root.configure(background='#d5ddf3') #border color

        # Create a style and configure the background color for ttk.Frame and checkbox
        style = ttk.Style()
        style.configure('My.TFrame', background='#d5ddf3')
        style.configure('My.TCheckbutton', background='#d5ddf3')
        
        self.create_sigma_input_frame()
        self.create_group_selection_frame()
        self.create_plot_options_frame()
        self.create_generate_button()
        self.create_error_label()
    
    def centre_window(self, width, height):
        """
        Center the window on the screen.

        Parameters:
        width (int): Width of the window.
        height (int): Height of the window.
        """
        ws = self.root.winfo_screenwidth()   # width of the screen
        hs = self.root.winfo_screenheight()  # height of the screen
        x = (ws / 2) - (width / 2)
        y = (hs / 2) - (height / 2)
        self.root.geometry('%dx%d+%d+%d' % (width, height, x, y))
        self.root.resizable(0, 0)
    
    def create_sigma_input_frame(self):
        """
        Create input fields for sigma values and trends/plunges.

        This method creates input fields for the user to enter the values of sigma1, sigma2, and sigma3.
        The input fields are arranged in a grid layout within a ttk.Frame.
        """
        input_frame = ttk.Frame(self.root, padding="10 10 10 10", style='My.TFrame')
        input_frame.grid(row=0, column=0, padx=10, pady=5, sticky='ew')

        # Create the grid layout for input boxes
        for i in range(5):
            input_frame.columnconfigure(i, weight=1)
    

        self.sigma1_entry = self.create_label_and_entry(input_frame, "Sigma\u2081 (σ₁) :", 0, 0)
        self.sigma2_entry = self.create_label_and_entry(input_frame, "Sigma\u2082 (σ₂) :", 1, 0)
        self.sigma3_entry = self.create_label_and_entry(input_frame, "Sigma\u2083 (σ₃) :", 2, 0)
        self.trend_sigma1_entry = self.create_label_and_entry(input_frame, "Trend σ₁:", 0, 2)
        self.plunge_sigma1_entry = self.create_label_and_entry(input_frame, "Plunge σ₁:", 0, 4)
        self.trend_sigma2_entry = self.create_label_and_entry(input_frame, "Trend σ₂:", 1, 2)
        self.plunge_sigma2_entry = self.create_label_and_entry(input_frame, "Plunge σ₂:", 1, 4)
        self.trend_sigma3_entry = self.create_label_and_entry(input_frame, "Trend σ₃:", 2, 2)
        self.plunge_sigma3_entry = self.create_label_and_entry(input_frame, "Plunge σ₃:", 2, 4)

    def create_label_and_entry(self, frame, text, row, column):
        """
        Helper function to create a label and entry widget.

        Parameters:
        frame (ttk.Frame): Frame to contain the label and entry.
        text (str): Text for the label.
        row (int): Row position in the grid.
        column (int): Column position in the grid.

        Returns:
        ttk.Entry: The created entry widget.
        """

        # Create and place the label
        label = ttk.Label(frame, text=text, background='#d5ddf3')
        label.grid(row=row, column=column, padx=10, pady=5, sticky='w')

        # Create and place the entry
        entry = ttk.Entry(frame)
        entry.grid(row=row, column=column + 1, padx=10, pady=5, sticky='w')
        return entry
    
    def create_plot_options_frame(self):
        """
        Create checkboxes for plot options.
        """

        plot_options_frame = tk.LabelFrame(self.root, text="Select Analysis and Data Display Option", 
                                           font=('Helvetica', 10, 'bold'), labelanchor='n',bg='#d5ddf3')
        plot_options_frame.grid(row=2, column=0, padx=10, pady=40, ipady=10, sticky='ew')

        # Ensure the frame itself stretches across the entire width of the window
        self.root.columnconfigure(0, weight=1)
        
        # Add weight to columns to distribute space evenly
        for i in range(3):
            plot_options_frame.columnconfigure(i, weight=1)

        self.plot_slip_var = tk.BooleanVar()
        self.plot_table_var = tk.BooleanVar()  

        plot_slip_check = ttk.Checkbutton(plot_options_frame, text='Slip Tendency (Ts)', variable=self.plot_slip_var, style='My.TCheckbutton')
        plot_table_check = ttk.Checkbutton(plot_options_frame, text='Show Data Table', variable=self.plot_table_var, style='My.TCheckbutton')  

        # Place the checkboxes on the grid and spread them evenly
        plot_slip_check.grid(row=1, column=0, padx=(0,10), pady=(10, 0))
        plot_table_check.grid(row=1, column=2, padx=(0,10), pady=(10, 0))

        # Add weight to columns to distribute space evenly
        for i in range(3):
            plot_options_frame.columnconfigure(i, weight=1)


    def create_group_selection_frame(self):
        """
        Create dropdown menu for selecting the group.

        This method creates a LabelFrame containing a Combobox for the user to select the group of data
        to analyze. The Combobox is populated with unique group names from the input DataFrame.
        """
        group_frame = tk.LabelFrame(self.root, text="Select Data to Analyse", bg='#d5ddf3')
        group_frame.grid(row=1, column=0, padx=10, pady=10)

        group_frame.columnconfigure(0, weight=1)
        group_frame.columnconfigure(1, weight=1)

        self.group_combo = ttk.Combobox(group_frame, state='readonly', width=30)
        self.group_combo.set("Select Data")
        self.group_combo.grid(row=0, column=1, padx=10, pady=5, sticky='w')

        self.update_dropdown_menu()

    def update_dropdown_menu(self):
        """
        Update the dropdown menu options based on the DataFrame.

        This method populates the Combobox with unique group names extracted from the 'group' column
        of the input DataFrame.
        """
        unique_groups = self.df['group'].unique().tolist()
        self.group_combo['values'] = unique_groups

    def create_generate_button(self):
        """
        Create the Generate Table button.

        This method creates a button that, when clicked, triggers the generation of the analysis table.
        """
        self.generate_button = tk.Button(self.root, background='#F7F9FE', default='active',
                                         height=2, width=25, relief="flat", 
                                         text="Estimate Stress and Tendency", command=self.generate_plots_table)
        self.generate_button.grid(row=3, column=0, padx=10,  pady=0)

    def create_error_label(self):
        """
        Create a label for error messages.

        This method creates a label widget that displays error messages to the user.
        """
        self.error_label = ttk.Label(self.root, text="", foreground="red", background='#d5ddf3')
        self.error_label.grid(row=4, column=0, padx=0, pady=(5, 5))

    def display_dataframe(self, df):
        """
        Display the DataFrame in a new Toplevel window.

        Parameters:
        df (pd.DataFrame): DataFrame to display.
        """
        df_window = tk.Toplevel(self.root)
        df_window.title("Estimated Stress and Slip Tendency")
        df_window.geometry("700x400")

        df_text = tk.Text(df_window, wrap='none')
        df_text.pack(expand=True, fill='both')

        # Add vertical scroll bar if necessary
        scrollbar_y = ttk.Scrollbar(df_window, orient='vertical', command=df_text.yview)
        scrollbar_y.pack(side='right', fill='y')
        df_text.configure(yscrollcommand=scrollbar_y.set)

        # Add horizontal scroll bar if necessary
        scrollbar_x = ttk.Scrollbar(df_window, orient='horizontal', command=df_text.xview)
        scrollbar_x.pack(side='bottom', fill='x')
        df_text.configure(xscrollcommand=scrollbar_x.set)

        df_text.insert(tk.END, df.to_string(index=False))

    def generate_plots_table(self):
        """
        Generate the analysis table based on the selected group and display it.

        This method retrieves the selected group, filters the input DataFrame based on the group,
        extracts the relevant data.
        
        Method creates stereonet plots for slip tendency. 
        A table of estimations is also created.
        """
        selected_group = self.group_combo.get()
        if selected_group == "Select group to analyze":
            self.error_label.config(text="Please select data to analyze.")
            return 

        # Filter the DataFrame based on the selected group
        filtered_df = self.df[self.df['group'].str.contains(selected_group, case=False)]
        dips = filtered_df['dip'].values
        dip_azimuths = filtered_df['dip_azimuth'].values


        # Validate input values
        try:
            sigma1 = float(self.sigma1_entry.get())
            sigma2 = float(self.sigma2_entry.get())
            sigma3 = float(self.sigma3_entry.get())
            trend1 = float(self.trend_sigma1_entry.get())
            plunge1 = float(self.plunge_sigma1_entry.get())
            trend2 = float(self.trend_sigma2_entry.get())
            plunge2 = float(self.plunge_sigma2_entry.get())
            trend3 = float(self.trend_sigma3_entry.get())
            plunge3 = float(self.plunge_sigma3_entry.get())

            # Check if any trend and plunge values are in correct ranges
            if not (0 <= plunge1 <= 90 and 0 <= plunge2 <= 90 and 0 <= plunge3 <= 90):
                self.error_label.config(text="Error: Plunge values must be between 0 and 90 degrees.")
                return
            if not (0 <= trend1 <= 360 and 0 <= trend2 <= 360 and 0 <= trend3 <= 360):
                self.error_label.config(text="Error: Trend values must be between 0 and 360 degrees.")
                return
        
        except ValueError:
            self.error_label.config(text="Error: Please add data to all input fields.")
            return

        # Clear error label if validation passes
        self.error_label.config(text="")

        stress_analysis = StressAnalysis(
            dips, dip_azimuths, sigma1, sigma2, sigma3,
            trend1, plunge1,
            trend2, plunge2,
            trend3, plunge3,
        )
        
        # Initialise the StereonetPlotter
        stereonet_plotter = StereonetPlotter(stress_analysis)

        # Define the dip and azimuth ranges
        dip_range = np.linspace(0, 90, 100)
        dip_azimuth_range = np.linspace(0, 360, 100)

        # Create the plots based on selected options
        if self.plot_slip_var.get():
            stereonet_plotter.create_tendency_plot(dip_range, dip_azimuth_range, plot_type='slip')
        
        # Display the dataframe if the table checkbox is selected
        if self.plot_table_var.get():
            self.display_dataframe(stress_analysis.df_tendencies)

        stereonet_plotter.show_all_plots()    
    
    def orientation_analysis(self, group):
        """
        Filter the DataFrame based on the selected group and return the relevant data.

        Parameters:
        group (str): The group to filter the data by.

        Returns:
        tuple: Filtered dip and dip azimuth values.
        """
        df2 = self.df[self.df['group'].str.contains(group, case=False)].copy()
        dip_azimuths = df2['dip_azimuth'].values
        dip = df2['dip'].values
        return dip, dip_azimuths

# Retrieve orientations from interpreations menu tree=
orient_list = vrgs.orientations()
df = pd.DataFrame(orient_list)
df.columns =['name','group', 'id', 'x', 'y', 'z', 'dip', 'dip_azimuth']

# Extract the last part of the 'group' column
df['group'] = df['group'].str.extract(r'Orientation/(.*?)/')[0]

# Build the Tkinter GUI App
root = tk.Tk()
app = StressGUI(root, df)
root.mainloop()
