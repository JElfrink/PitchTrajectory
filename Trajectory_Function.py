""" 
Script that creates the ball trajectory of a thrown baseball pitch based on the release parameters

Script is based on Alan Nathan's Trajectory Calculator (http://baseball.physics.illinois.edu/trajectory_calculator-new3D.html) and the paper "Analytical Solution of the Pitching Trajectory Equations with Laminar Air Flow" by Alan M. Nathan

The calculation uses an approach that recalculates the acceleration at each time step based on the current velocity and position of the ball. This allows for a more accurate simulation of the ball's trajectory.
There are more complex models that take into account the Magnus- and/or Seam-Shifted effects, but this model is a good approximation for most pitches.

Script created by: Jasper Elfrink

"""

import math
import pandas as pd
import numpy as np



def simulate_pitch_trajectory(x0, y0, z0, release_speed, release_angle, release_direction, Handedness, pitch_type, backspin, sidespin, gyrospin):
 # Timeseries properties
    tau = 10000 # in seconds
    dt = 0.001 # in seconds




    ## Constants necessary

    # Baseball Properties
    mass = 5.125 # Mass of baseball in oz
    circumference = 9.125 # Circumference of baseball in inches
    cd0 = 0.3008 # Drag coefficient
    cdspin = 0.0292 # Spin coefficient

    cl0 = 0.583 # Lift coefficient
    cl1 = 2.333 # Lift coefficient
    cl2 = 1.120 # Lift coefficient


    # Environmental properties
    tempF = 70  # in degrees Fahrenheit
    temperature = (5/9)*(tempF-32) # Temperature in degrees Celsius
    elevFt = 15 # in feet
    elevation = elevFt/3.2808 # Elevation in meters
    vwind = 5 # in mph
    phiwind = 90 # in degrees
    hwind = 10 # in feet
    relative_humidity = 50 # in %
    barometric_pressure_Hg = 29.92 # in Hg


    vxw = vwind*1.467*math.sin(phiwind*math.pi/180) # in ft/s
    vyw = vwind*1.467*math.cos(phiwind*math.pi/180) # in ft/s
    svp = 4.5841*math.exp((18.687-temperature/234.5)*temperature/(257.14+temperature)) # in mmHg
    barometric_pressure = barometric_pressure_Hg*1000/39.37 # in mmHg

    beta = 0.0001217 # actual pressure

    rho_air = 1.2929*(273/(temperature+273)*(barometric_pressure*math.exp(-beta*elevation)-0.3783*relative_humidity*svp/100)/760) # in kg/m^3
    rho = rho_air*0.06261 # in lb/ft^3
    #rho_air = 1.194 # in kg/m^3
    #rho = 0.0747    # in lb/ft^3

    reynolds_number = rho_air*44.7*(circumference/(math.pi*39.37))*(temperature+273.16+120)/(0.000001512*(temperature+273.16)**1.5) # reynolds number for v = 100mph

    c0 = 0.07182*rho*(5.125/mass)*(circumference/9.125)**2 # in lb*s^2/ft^4
    #c0 = 0.005368
    const = c0 # in lb*s^2/ft^4

    ## Parameters based on manual inputs

    # Handedness
    if Handedness == "R":           # 1 for right handed, -1 for left handed
        sign = 1
    else:
        sign = -1 

    v0 = release_speed*1.467 # in ft/s
    v0x = 1.467*release_speed*math.cos(release_angle*math.pi/180)*math.sin(release_direction*math.pi/180) # in ft/s
    v0y = -1.467*release_speed*math.cos(release_angle*math.pi/180)*math.cos(release_direction*math.pi/180) # in ft/s
    v0z = 1.467*release_speed*math.sin(release_angle*math.pi/180) # in ft/s

    wx = (-backspin*math.cos(release_direction*math.pi/180)- sidespin*math.sin(release_angle*math.pi/180)*math.sin(release_direction*math.pi/180) + gyrospin*v0x/v0)*math.pi/30 # in rad/s
    wy = (backspin*math.sin(release_direction*math.pi/180)- sidespin*math.sin(release_angle*math.pi/180)*math.cos(release_direction*math.pi/180) + gyrospin*v0y/v0)*math.pi/30 # in rad/s
    wz = (sidespin*math.cos(release_angle*math.pi/180)+ gyrospin*v0z/v0)*math.pi/30 # in rad/s


    omega = math.sqrt(backspin**2+sidespin**2)*math.pi/30+0.001 # initial spin in rpm
    romega = (circumference/2/math.pi)*omega/12 # initial spin in ft/s, presumes ball radius = circumference/(2*pi)


    total_spin = math.sqrt(backspin**2+sidespin**2+gyrospin**2)+0.001 # in rpm
    cd = cd0 + cdspin*total_spin/1000 # drag coefficient
    flag = 1 # Flag=0 means full spin is used for both spin-dependent drag and S; Flag=1 means only the transverse spin is used for both spin-dependent drag and S

    ## Set 0 parameters of all variables in the timeseries

    t = 0
    x = x0
    y = y0
    z = z0
    r = math.sqrt(x**2+y**2)
    phi = 0
    vx = v0x
    vy = v0y
    vz = v0z
    v = math.sqrt(vx**2+vy**2+vz**2)
    vmph = v/1.467
    w_perp = math.sqrt(total_spin**2-flag*((30/math.pi)*(wx*vx+wy*vy+wz*vz)/v)**2) # perpendicular spin in rpm
    r_omega_perp = (w_perp*math.pi/30)*(circumference/(2*math.pi))/12 # perpendicular spin in ft/s


    if z >= hwind:
        vw = math.sqrt((vx-vxw)**2+(vy-vyw)**2+vz**2) # in ft/s
        vxw = vxw
        vyw = vyw

    else:
        vw = v # in ft/s
        vxw = 0
        vyw = 0

    Cd = cd0 + (cdspin*w_perp/1000)*math.exp(-t/(tau*146.7/vw)) 
    S = (r_omega_perp/vw)*math.exp(-t/(tau*146.7/vw)) # in ft^2
    cl = cl2*S/(cl0+cl1*S) # lift coefficient


        

    adragx = -const*Cd*vw*(vx-vxw)
    adragy = -const*Cd*vw*(vy-vyw)
    adragz = -const*Cd*vw*vz

    w_perp_over_w = r_omega_perp/romega
    aMagx = const*(cl/omega)*vw*(wy*vz-wz*(vy-vyw))/w_perp_over_w
    aMagy = const*(cl/omega)*vw*(wz*(vx-vxw)-wx*vz)/w_perp_over_w
    aMagz = const*(cl/omega)*vw*(wx*(vy-vyw)-wy*(vx-vxw))/w_perp_over_w

    ax = adragx + aMagx
    ay = adragy + aMagy
    az = adragz + aMagz - 32.174



    ## Create timeseries
    time_series = [t]
    x_series = [x]
    y_series = [y]
    z_series = [z]
    vx_series = [vx]
    vy_series = [vy]
    vz_series = [vz]
    ax_series = [ax]
    ay_series = [ay]
    az_series = [az]
    adragx_series = [adragx]
    adragy_series = [adragy]
    adragz_series = [adragz]
    aMagx_series = [aMagx]
    aMagy_series = [aMagy]
    aMagz_series = [aMagz]
    w_perp_series = [w_perp]
    Cd_series = [Cd]
    vw_series = [vw]
    flag_series = [flag]
    vxw_series = [vxw]
    vyw_series = [vyw]
    S_series = [S]
    cl_series = [cl]

    previous_y = y

    while y >= -2:
        t = t+dt
        x = x+vx*dt + 0.5*ax*dt*dt
        y = y+vy*dt + 0.5*ay*dt*dt
        z = z+vz*dt + 0.5*az*dt*dt
        r = math.sqrt(x**2+y**2)
        phi = math.atan2(y,x)*180/math.pi
        vx = vx+ax*dt
        vy = vy+ay*dt
        vz = vz+az*dt
        v = math.sqrt(vx**2+vy**2+vz**2)
        vmph = v/1.467
        w_perp = math.sqrt(total_spin**2-flag*((30/math.pi)*(wx*vx+wy*vy+wz*vz)/v)**2) # perpendicular spin in rpm
        r_omega_perp = (w_perp*math.pi/30)*(circumference/(2*math.pi))/12 # perpendicular spin in ft/s
        
        if z >= hwind:
            vw = math.sqrt((vx-vxw)**2+(vy-vyw)**2+vz**2) # in ft/s
            vxw = vxw
            vyw = vyw

        else:
            vw = v # in ft/s
            vxw = 0
            vyw = 0

        Cd = cd0 + (cdspin*w_perp/1000)*math.exp(-t/(tau*146.7/vw)) 
        S = (r_omega_perp/vw)*math.exp(-t/(tau*146.7/vw)) # in ft^2
        cl = cl2*S/(cl0+cl1*S) # lift coefficient


        

        adragx = -const*Cd*vw*(vx-vxw)
        adragy = -const*Cd*vw*(vy-vyw)
        adragz = -const*Cd*vw*vz

        w_perp_over_w = r_omega_perp/romega
        aMagx = const*(cl/omega)*vw*(wy*vz-wz*(vy-vyw))/w_perp_over_w
        aMagy = const*(cl/omega)*vw*(wz*(vx-vxw)-wx*vz)/w_perp_over_w
        aMagz = const*(cl/omega)*vw*(wx*(vy-vyw)-wy*(vx-vxw))/w_perp_over_w

        ax = adragx + aMagx
        ay = adragy + aMagy
        az = adragz + aMagz - 32.174
        
        #flag = 1 if ((previous_y - 17/12)*(y - 17/12)) < 0 else 0
        

        # Update time series
        time_series.append(t)
        x_series.append(x)
        y_series.append(y)
        z_series.append(z)
        vx_series.append(vx)
        vy_series.append(vy)
        vz_series.append(vz)
        ax_series.append(ax)
        ay_series.append(ay)
        az_series.append(az)
        adragx_series.append(adragx)
        adragy_series.append(adragy)
        adragz_series.append(adragz)
        aMagx_series.append(aMagx)
        aMagy_series.append(aMagy)
        aMagz_series.append(aMagz)
        w_perp_series.append(w_perp)
        Cd_series.append(Cd)
        vw_series.append(vw)
        flag_series.append(flag)
        vxw_series.append(vxw)
        vyw_series.append(vyw)
        S_series.append(S)
        cl_series.append(cl)


        # Update previous_y for the next iteration
        previous_y = y

        



    data = {
        "Time": time_series,
        "X": x_series,
        "Y": y_series,
        "Z": z_series,
        "vx": vx_series,
        "vy": vy_series,
        "vz": vz_series,
        "ax": ax_series,
        "ay": ay_series,
        "az": az_series,
        "adragx": adragx_series,
        "adragy": adragy_series,
        "adragz": adragz_series,
        "aMagx": aMagx_series,
        "aMagy": aMagy_series,
        "aMagz": aMagz_series,
        "w_perp": w_perp_series,
        "Cd": Cd_series,
        "vw": vw_series,
        "flag": flag_series,
        "vxw": vxw_series,
        "vyw": vyw_series,
        "S": S_series,
        "cl": cl_series

    }

    return pd.DataFrame(data)