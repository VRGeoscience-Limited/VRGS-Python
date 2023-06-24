import vrgs
import numpy as np
import matplotlib.cm as cm
import mplstereonet
import matplotlib.pyplot as plt

# Define constants
BIN_SIZE = 10
DEGREES_IN_CIRCLE = 360
HALF_CIRCLE = 180

def plot_pole_of_planes(strike, dip):
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='stereonet')
    ax.pole(strike, dip, c='k', label='Pole of the Planes')
    ax.density_contourf(strike, dip, measurement='poles', cmap='Reds')
    ax.set_title('Density coutour of the Poles', y=1.10, fontsize=15)
    ax.grid()
    plt.show()

def plot_rose_diagram(two_halves):
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='polar')
    ax.bar(np.deg2rad(np.arange(0, DEGREES_IN_CIRCLE, BIN_SIZE)), two_halves, 
           width=np.deg2rad(BIN_SIZE), bottom=0.0, color='.8', edgecolor='k')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_thetagrids(np.arange(0, DEGREES_IN_CIRCLE, BIN_SIZE), labels=np.arange(0, DEGREES_IN_CIRCLE, BIN_SIZE))
    ax.set_rgrids(np.arange(1, two_halves.max() , 10), angle=0, weight= 'black')
    ax.set_title('Rose Diagram of Strikes', y=1.10, fontsize=15)
    plt.show()

orient_list = vrgs.orientations()
print(orient_list) 
name, id, x_pos, y_pos, z_pos, dip, azimuth = zip(*orient_list)

# Compute strike from azimuth
strike = np.asarray(azimuth) - HALF_CIRCLE

# Histogram of strikes
bin_edges = np.arange(-5, DEGREES_IN_CIRCLE + BIN_SIZE, BIN_SIZE)
number_of_strikes, bin_edges = np.histogram(strike, bin_edges)

# Wrap around
number_of_strikes[0] += number_of_strikes[-1]

# Split and concatenate
half = np.sum(np.split(number_of_strikes[:-1], 2), 0)
two_halves = np.concatenate([half, half])

plot_pole_of_planes(strike, dip)
plot_rose_diagram(two_halves)