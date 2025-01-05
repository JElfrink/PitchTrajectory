"""
This script is for the creation of a streamlit app that will allow users to explore pitch trajectories based on data that they can input themselves.

"""
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import math
from Trajectory_Function import simulate_pitch_trajectory

# Title
st.title('Pitch Trajectory Analysis')

### Sidebar
st.sidebar.title('Input Parameters')

## Create sliders for Input Parameters
st.sidebar.markdown('**Pitch Parameters**')
# Release Location parameters
x0 = st.sidebar.slider('Horizontal Release Location (x0 in ft)', min_value=-5.0, max_value=5.0, value=0.0, step=0.1, help = 'in statcast x0 = release_pos_x' )
y0 = st.sidebar.slider('Release Extension (y0 in ft)', min_value=50.0, max_value=66.0, value=55.0, step=0.1, help = 'in statcast y0 = 60.5-release_extension')
z0 = st.sidebar.slider('Release Height (z0 in ft)', min_value=5.0, max_value=7.0, value=6.0, step=0.1, help = 'in statcast z0 = release_pos_z')

# Pitch Parameters
release_speed = st.sidebar.slider('Release Speed (mph)', min_value=65.0, max_value=105.0, value=90.0, step=0.1, help = 'If using statcast -> release_speed = math.sqrt(vxR**2+vyR**2+vzR**2)*0.6818 # in mph')
release_angle = st.sidebar.slider('Release Angle (degrees)', min_value=-10.0, max_value=10.0, value=-3.0, step=0.1, help = 'If using statcast -> release_angle = math.atan(vzR/math.sqrt(vxR**2+vyR**2))*(180/math.pi) # in degrees')
release_direction = st.sidebar.slider('Release Direction (degrees)', min_value=-8.0, max_value=8.0, value=0.0, step=0.1, help = 'If using statcast -> release_direction = math.atan(-vxR/vyR)*180/(math.pi) # in degrees')
Handedness = st.sidebar.selectbox('Pitcher Handedness', ['R', 'L'], help = 'R for Right Handed Pitcher, L for Left Handed Pitcher')
pitch_type = st.sidebar.selectbox('Pitch Type', ['Fastball', 'Curveball', 'Slider', 'Changeup', 'Cutter', 'Sinker', 'Splitter', 'Knuckleball', 'Eephus', 'Other'])
backspin = st.sidebar.slider('Backspin (rpm)', min_value=-3000, max_value=3000, value=2000, step=10, help = 'sometimes referred to as wb')
sidespin = st.sidebar.slider('Sidespin (rpm)', min_value=-3000, max_value=3000, value=0, step=10, help = 'sometimes referred to as ws')
gyrospin = st.sidebar.slider('Gyrospin (rpm)', min_value=-2000, max_value=2000, value=0, step=10, help = 'sometimes referred to as wg')

# Run simulation
df = simulate_pitch_trajectory(x0, y0, z0, release_speed, release_angle, release_direction, Handedness, pitch_type, backspin, sidespin, gyrospin)

### Plot the trajectory in a 3D plot
import plotly.graph_objs as go
import math
import numpy as np
import streamlit as st
import time 

# Create data
x = df['X']
y = df['Y']
z = df['Z']



# Create a line for the trajectory
trajectory_line = go.Scatter3d(
    x=x,
    y=y,
    z=z,
    mode='markers',
    marker=dict(
        color='indianred',
        size = 6,
        opacity = 0.5
    )
)

# Create a line for the trajectory
trajectory_line_k_zone = go.Scatter3d(
    x=x,
    y=y,
    z=z,
    mode='markers',
    marker=dict(
        color='indianred',
        size = 12,
        opacity = 0.5
    )
)




# Create  the strike zone
strike_zone = go.Scatter3d(
    x=[8.5/12, -8.5/12, -8.5/12,    8.5/12,     8.5/12,     8.5/12, -8.5/12,    -8.5/12, 8.5/12,    8.5/12, -8.5/12,    -8.5/12,    -8.5/12, -8.5/12, 8.5/12, 8.5/12],
    y=[17/12,   17/12,   0,         0,          17/12,      17/12,  17/12,      0,       0,         17/12,  17/12,      17/12,      0,      0, 0, 0],
    z=[6*0.3,   6*0.3,   6*0.3,     6*0.3,      6*0.3,     6*0.77, 6*0.77,     6*0.77,  6*0.77,     6*0.77,  6*0.77,     6*0.3,     6*0.3,  6*0.77, 6*0.77, 6*0.3],
    mode = 'lines',
    line = dict(
        color = 'blue',
        width = 2
    )
)

# Create the home plate
home_plate = go.Scatter3d(
    x=[8.5/12, -8.5/12, -8.5/12, 0, 8.5/12, 8.5/12, 8.5/12],
    y=[17/12,17/12, 8.5/12,0, 8.5/12, 17/12, 8.5/12],
    z=[0,0,0,0,0,0],
    mode = 'lines',
    line = dict(
        color = 'white',
        width = 2
    )
)

# Create the mound
# Generate a grid of points in spherical coordinates
phi = np.linspace(0, np.pi / 2, 50)  # azimuthal angle
theta = np.linspace(0, 2 * np.pi, 100)  # polar angle
phi_grid, theta_grid = np.meshgrid(phi, theta)

# Convert spherical coordinates to Cartesian coordinates
x = 9 * np.sin(phi_grid) * np.cos(theta_grid)
y = 9 * np.sin(phi_grid) * np.sin(theta_grid) + 60.5
z = 1 * np.cos(phi_grid)

# Clay color command
color1 = 'coral'
# Create a filled surface plot
mound = go.Surface(
    x=x,
    y=y,
    z=z,
    colorscale=[[0, color1], [1,color1]],  # Adjust the colorscale as desired
    showscale = False
)

# Create a pitcher plate
pitcher_plate = go.Scatter3d(
    x=[1,1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,1],
    y=[60.5, 61,61,60.5,60.5,60.5,61,61,61,61,61,61,60.5,60.5,60.5,60.5],
    z=[1, 1,1,1,1,1.2,1.2,1,1.2,1.2,1,1.2,1.2,1,1.2,1.2],
    mode = 'lines',
    line = dict(
        color = 'white',
        width = 2
    )
)

# Create the figure and display
fig = go.Figure(data=[pitcher_plate, mound, home_plate, strike_zone, trajectory_line])



# Show the full field
fig.update_layout(
    scene = dict(
        zaxis = dict(
            range=[0, 10],
            # green ground color

            backgroundcolor="rgb(0, 128, 0)",
            
            showbackground=True
        ),
        xaxis = dict(
            range=[-2, 2]

        ),
        yaxis = dict(
            range=[-2,80],
            backgroundcolor="cadetblue",
            showbackground=True
            # scale of the y axis is 1 ft = 1 unit
            

            

        ),
        # create the figure to be a rectangle
        aspectratio=dict(x=0.5, y=2, z=1),
        
    ),
    showlegend = False
    
)

st.write('Pitch Trajectory')
st.plotly_chart(fig)



k_zone = go.Figure(data=[pitcher_plate, mound, home_plate, strike_zone, trajectory_line_k_zone])
k_zone.update_layout(
    scene = dict(
        zaxis = dict(
            range=[0, 6*0.77+0.5]
        ),
        xaxis = dict(
            range=[-9.5/12-0.5, 9.5/12+0.5]
        ),
        yaxis = dict(
            range=[0, 17/12]
        ),
        aspectmode = 'cube',
        
    ),
    showlegend = False
)
st.write('Strike Zone')
st.plotly_chart(k_zone)

# Display the release parameters
st.write('Release Parameters')
st.write('Release Location (x0): ', x0)
st.write('Release Extension (y0): ', y0)
st.write('Release Height (z0): ', z0)

# Display the data of the pitch where it crosses the strike zone so where y = 0
if st.checkbox('Display Pitch Data at the Strike Zone'):
    df_strike_zone = df[df['Y'] <= 0]
    st.write('Pitch Data at the Strike Zone')
    st.write(df_strike_zone)

# Give the option to display the data in a table
if st.checkbox('Display Full Data'):
    st.write('Full Pitch Data')
    st.write(df)


# Create a button to store the data in a csv file
if st.button('Save Pitch Data'):
    df.to_csv('pitch_data.csv')
    st.write('Data Saved as pitch_data.csv')
