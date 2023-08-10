#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 21:20:02 2021

@author: DM

TODO: Non-integer hour time zones.
"""

#############################################
# Address: Lat, Lon, Elev (km), Timezone (HH)
# 5544 MapleC: +42.2371552, -77.8020070, 0.675, -5 -4Daylight
# 351 Willow : +42.0106135, -71.2165150, 0.040, -5 -4Daylight
# SCMARS:      +41.965747,  -71.186538,  0.026, -5 -4Daylight
mylat, mylon, myele, mytz = +41.965747,  -71.186538, 0.026, -4
#mylat, mylon, myele, mytz = 0, 0, 0, 0

# window_azimuth_deg is defined such that window_azimuth_deg = (0,90,-90) 
# for window facing (S,E,W) respectively
window_azimuth_deg = 0.0

# Refractive index of the material between the hole and the screen
mu = 1.495

gnomon_sz = 1.0

# Elevation cutoff for the Sun
elvcutoff_deg = 15.0

# Shadow is TOO_BIG if elevation criteria is not met.
TOO_BIG = 9.9e9


# Dates of winter and summer solstice
w_sol="2021-12-21"
s_sol="2022-06-21"
w_sol2="2022-12-21"
#############################################

#############################################
# Import required libraries
#import sys
import numpy as np
from astroquery.jplhorizons import Horizons
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#############################################

#############################################
# Globals

# viewer =  1 if the viewer is on the same side of dial's plane as the gnomon
# viewer = -1 if the viewer is on the other side of the dial's plane.
viewer = -1     # Looking from behind

myloc = {'lon': mylon, 'lat': mylat, 'elevation': myele}
lon_flag = "E"
lat_flag = "N"
if mylon < 0 :
    lon_flag = "W"
if mylat < 0 :
    lat_flag = "S"
deg = u"\N{DEGREE SIGN}"
#loc_print = 'Subtract one hour if daylight savings is NOT in effect.'
loc_print = "Customized for Wheaton\n" +str(mylat)+ deg + lat_flag + ", " 
loc_print = loc_print + str(np.abs(mylon)) + deg + lon_flag + ", "
loc_print = loc_print + "elevation " + str(int(myele*3280.84)) + "'"#, and"
#loc_print = loc_print + "\n when the Sun is higher than "
#loc_print = loc_print + str(int(elvcutoff_deg)) + deg
#loc_print = loc_print + " above horizon and less than 5 hours from local noon."
#sys.exit(0)
#############################################
    

def tilted_xy (lat, az, el, el_cutoff):
# Input: latitude of location, and azimuth and elevation of an object.
# Output: xy-coordinates of the gnomon's shadow in a 
# coordinate system where x-axis is towards the east, and y-axis is towards
# the north celestial pole, and z-axis is towards the zenith.
# NOTE: Input angles in radians, and output lengths in units of the
# gnomon's height.
#
    # Don't bother doing anything if elevation is too low
    if (el<el_cutoff):
        return TOO_BIG, TOO_BIG
    
    # Given the azimuth and elevation, the polar angle (theta) and 
    # azimuthal angle (phi) in a coordinate system where 
    # theta=0 => zenith/Z-axis, and  phi=0 =>east/X-axis and 
    # phi=90 => north/Y-axis.
    theta = np.pi/2 - el
    phi = np.pi/2 - az
    
    # The xyz-components of the unit vector in this coordinate system are
    x = np.sin(theta)*np.cos(phi)
    y = np.sin(theta)*np.sin(phi)
    z = np.cos(theta)
    
    # Rotating about the x-axis by 'lat', the xyz-prime-components in the 
    # rotated coordinate system are
    cl, sl = np.cos(lat), np.sin(lat)
    xp = x
    yp =  cl*y + sl*z
    zp = -sl*y + cl*z
    
    # The theta-prime and phi-prime coordinates are
    thetap = np.arctan2(np.sqrt(xp*xp+yp*yp), zp)
    phip = np.arctan2(yp, xp)
    #print('thetap= ',np.rad2deg(thetap), 'phip= ',np.rad2deg(phip))
    
    # Don't do anything if the polar angle from the z-prime axis is 
    # too large because the shadow will be too long
    thetap_cutoff = np.pi/2 - el_cutoff
    if (thetap>thetap_cutoff):
        return TOO_BIG, TOO_BIG
    
    # The xp and yp-components of the shadow of a gnomon of unit height,
    # located at the origin and the gnomon pointed towards +zp-axis are
    phi_shadow = phip + np.pi
    theta_shadow = np.arcsin(np.sin(thetap)/mu)
    rr = gnomon_sz * np.tan(theta_shadow)
    xx = rr * np.cos(phi_shadow)
    yy = rr * np.sin(phi_shadow)
    
    return xx, yy


def vertical_xy (waz, az, el, el_cutoff):
# Input: waz=window azimuth, and azimuth and elevation of an object.
# waz is defined such that window_az = waz = (0,pi/2,-pi/2) for window
# facing (S,E,W) respectively.
# Output: xy-coordinates of the gnomon's shadow in a 
# coordinate system where x-axis is towards the east, and y-axis is towards
# the zenith.
# NOTE: Input angles in radians, and output lengths in units of the
# gnomon's height.
#
    # Don't bother doing anything if elevation is too low
    if (el<el_cutoff):
        return TOO_BIG, TOO_BIG
    
    # Given the azimuth and elevation, the polar angle (theta) and 
    # azimuthal angle (phi) in a coordinate system where 
    # theta=0 => zenith/Z-axis, and  phi=0 =>east/X-axis and 
    # phi=90 => north/Y-axis.
    theta = np.pi/2 - el
    phi = np.pi/2 - az
    
    # The xyz-components of the unit vector in this coordinate system are
    x = np.sin(theta)*np.cos(phi)
    y = np.sin(theta)*np.sin(phi)
    z = np.cos(theta)
    
    # Rotating about the x-axis by 90-deg, the xyz-prime-components in this
    # rotated coordinate system where the +x-axis is towards East, +y towards
    # Zenith, and +z towards South, are
    cl, sl = 0,1    #cos(pi/2), sin(pi/2)
    xp = x
    yp =  cl*y + sl*z
    zp = -sl*y + cl*z
    
    # Further rotating about the y-axis by an angle window_az, where
    # window_az is defined such that window_az = waz = (0,90,-90) for window
    # facing (S,E,W) respectively. Basically rotating CCW from +z-axis towards
    # +x-axis. The double-primed coordinates are related to the primed 
    # coordinates via
    # zpp     cos(waz)   sin(waz)    0       zp
    # xpp =  -sin(waz    cos(waz)    0       xp
    # ypp        0          0        1       yp
    cwaz, swaz = np.cos(waz), np.sin(waz)
    zpp = cwaz*zp + swaz*xp
    xpp = cwaz*xp - swaz*zp
    ypp = yp
    
    # The theta'' and phi'' coordinates are
    thetapp = np.arctan2(np.sqrt(xpp*xpp+ypp*ypp), zpp)
    phipp = np.arctan2(ypp, xpp)
    
    # Don't do anything if the polar angle from the z'' axis is 
    # too large because the shadow will be too long
    thetapp_cutoff = np.pi/2 - el_cutoff
    if (thetapp>thetapp_cutoff):
        return TOO_BIG, TOO_BIG
    
    # The xp and yp-components of the shadow of a gnomon of unit height,
    # located at the origin and the gnomon pointed towards +zp-axis are
    phi_shadow = phipp + np.pi
    theta_shadow = thetapp
    rr = gnomon_sz * np.tan(theta_shadow)
    xx = rr * np.cos(phi_shadow) * viewer
    yy = rr * np.sin(phi_shadow)
    
    return xx, yy


def draw_EWNS_unit_line (xpos, ypos):
    plt.plot([xpos-0.5,xpos+0.5], [ypos,ypos], 'k-')
    plt.plot([xpos,xpos], [ypos-0.05,ypos+0.05], 'k-')


def get_hr_label_pos (x,y):
    ymax = np.amax(y)
    ylab = ymax + 0.1
    xlab = x[np.argmax(y)]
    lab_rot = 0
    return xlab, ylab, np.rad2deg(lab_rot)


def get_pos(start_day, stop_day, hr, mn):
    mystrt  = start_day + " " + str(hr).zfill(2) + ":" + str(mn).zfill(2)
    mystop  = stop_day  + " " + str(hr).zfill(2) + ":" + str(mn).zfill(2)
    myepoch = {'start': mystrt, 'stop':  mystop, 'step':'1d'}
    mysun   = Horizons(id='Sun', location=myloc, epochs=myepoch, id_type='id')
    eph     = mysun.ephemerides(refraction=True)
    ndays = len(eph['AZ'])
    myxx, myyy = [], []
    for ii in range(ndays):
        tmpx, tmpy = vertical_xy(np.deg2rad(window_azimuth_deg), 
                               np.deg2rad(eph['AZ'][ii]),  
                               np.deg2rad(eph['EL'][ii]), 
                               np.deg2rad(elvcutoff_deg))
    
        if ((tmpx != TOO_BIG) and (tmpy != TOO_BIG)):
            myxx.append(tmpx)
            myyy.append(tmpy)
            
    # To ensure the arrays have finite length, add negative bogus at the end
    if len(myxx) < 1:
        myxx.append(TOO_BIG)
        myyy.append(TOO_BIG)
        
    return np.array(myxx), np.array(myyy)


xpgsz = gnomon_sz*7.6763923241404814
ypgsz = gnomon_sz*8.053883609691486
print('\n PDF page size: ',xpgsz,' x ',ypgsz,'inches. Ratio= ', xpgsz/ypgsz)

opfile = 'wheaton_vertical_sundial.pdf'

font = {'family' : 'cursive',
        'weight' : 'bold',
        'size'   : int(gnomon_sz*10)}
matplotlib.rc('font', **font)

fig = plt.figure(figsize=(xpgsz, ypgsz), dpi=600)
ax = fig.add_subplot(111)
props = {'ha': 'center', 'va': 'center'}
ax.set_aspect('equal', adjustable='box')

# Do uttarayan
strt_day = w_sol
stop_day = s_sol
print('\n Making Uttarayan calculations')

skipquarters = 0
xmin, xmax = TOO_BIG, -TOO_BIG
ymin, ymax = TOO_BIG, -TOO_BIG
for hh in range(6, 19 + 1):
    print(' ',hh)
    xx,yy = get_pos(strt_day, stop_day, hh-mytz, 00)
    if (len(xx)>1):
        ax.plot(xx,yy,'k-', linewidth=2)
        xlab,ylab,lab_rot = get_hr_label_pos(xx, yy)
        ax.text(xlab, ylab, hh, props, rotation=lab_rot)
        
        if xmin>np.amin(xx):
            xmin=np.amin(xx)
        if ymin>np.amin(yy):
            ymin=np.amin(yy)
        if xmax<np.amax(xx):
            xmax=np.amax(xx)
        if ymax<np.amax(yy):
            ymax=np.amax(yy)

    if skipquarters == 0:
        for mm in range(15, 46, 15):
            xx,yy = get_pos(strt_day, stop_day, hh-mytz, mm)
            if (len(xx)>1):
                ax.plot(xx,yy,'k-', linewidth=0.5)
                if xmin>np.amin(xx):
                    xmin=np.amin(xx)
                if ymin>np.amin(yy):
                    ymin=np.amin(yy)
                if xmax<np.amax(xx):
                    xmax=np.amax(xx)
                if ymax<np.amax(yy):
                    ymax=np.amax(yy)

draw_EWNS_unit_line(0, 0)
circleu = plt.Circle((0, 0), 0.05, color='r', fill=False)
plt.gca().add_patch(circleu)

# Do dakshinayan
strt_day = s_sol
stop_day = w_sol2
print('\n Making Dakshinayan calculations')

# Offset dakshinayan plots downward by this amount
yoff_dk = 1.3*(ymax - ymin) 

skipquarters = 0
for hh in range(6, 19 + 1):
    print(' ',hh)
    xx,yy = get_pos(strt_day, stop_day, hh-mytz, 00)
    yy = yy - yoff_dk
    if (len(xx)>1):
        ax.plot(xx,yy,'k-', linewidth=2)
        xlab,ylab,lab_rot = get_hr_label_pos(xx, yy)
        ax.text(xlab, ylab, hh, props, rotation=lab_rot)
        
        if xmin>np.amin(xx):
            xmin=np.amin(xx)
        if ymin>np.amin(yy):
            ymin=np.amin(yy)
        if xmax<np.amax(xx):
            xmax=np.amax(xx)
        if ymax<np.amax(yy):
            ymax=np.amax(yy)

    if skipquarters == 0:
        for mm in range(15, 46, 15):
            xx,yy = get_pos(strt_day, stop_day, hh-mytz, mm)
            yy = yy - yoff_dk
            if (len(xx)>1):
                ax.plot(xx,yy,'k-', linewidth=0.5)
                if xmin>np.amin(xx):
                    xmin=np.amin(xx)
                if ymin>np.amin(yy):
                    ymin=np.amin(yy)
                if xmax<np.amax(xx):
                    xmax=np.amax(xx)
                if ymax<np.amax(yy):
                    ymax=np.amax(yy)

draw_EWNS_unit_line(0, 0-yoff_dk)
circled = plt.Circle((0, 0-yoff_dk), 0.05, color='r', fill=False)
plt.gca().add_patch(circled)

xrange = xmax - xmin
yrange = ymax - ymin
ymax = 0    # For vertical ones
xlo = xmin - 0.05*xrange 
ylo = ymin - 0.05*yrange 
xhi = xmax + 0.05*yrange
yhi = ymax + 0.05*yrange
                
ax.set_xlim([xlo, xhi])
ax.set_ylim([ylo, yhi])

xpgsz = xhi - xlo
ypgsz = yhi - ylo
print(xpgsz, ypgsz, xpgsz/ypgsz)

# Message which dial to use at what time
msg_txt = 'Use top dial from winter solstice to summer solstice. '
msg_txt = msg_txt + 'Bottom dial from summer solstice to winter solstice.\n'
msg_txt = msg_txt + 'Subtract one hour if dayligt savings is NOT in effect.'
message_xpos = 0.5*(xlo+xhi)
message_ypos = (ylo + yhi)/2 + 0.00*yrange
ax.text(message_xpos, message_ypos, msg_txt, props)


loc_print_xpos = message_xpos
loc_print_ypos = ylo + 0.05*yrange
ax.text(loc_print_xpos, loc_print_ypos, loc_print, props)

pdf_pages = PdfPages(opfile)

fig.subplots_adjust(left = 0., bottom = 0., right = 1., top = 1.)
pdf_pages.savefig(fig)
pdf_pages.close()

#sys.exit(0)
