#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 21:20:02 2021

@author: DM

TODO: Non-integer hour time zones.
"""

#############################################
# Address: Lat, Lon, Elev (km), Timezone (HH)
# 5544 MapleC: +42.2371552, -77.8020070, 0.675, -5
# 351 Willow : +42.0106135, -71.2165150, 0.040, -5
mylat, mylon, myele, mytz = +42.2371552, -77.8020070, 0.675, -4
#mylat, mylon, myele, mytz = 0, 0, 0, 0

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
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#############################################

#############################################
# Globals
myloc = {'lon': mylon, 'lat': mylat, 'elevation': myele}
lon_flag = "E"
lat_flag = "N"
if mylon < 0 :
    lon_flag = "W"
if mylat < 0 :
    lat_flag = "S"
deg = u"\N{DEGREE SIGN}"
#loc_print = 'Subtract one hour if daylight savings is NOT in effect.'
loc_print = "Customized for " +str(mylat)+ deg + lat_flag + ", " 
loc_print = loc_print + str(np.abs(mylon)) + deg + lon_flag + ", "
loc_print = loc_print + "elevation " + str(int(myele*3280.84)) + "', and"
loc_print = loc_print + "\n when the Sun is higher than "
loc_print = loc_print + str(int(elvcutoff_deg)) + deg
loc_print = loc_print + " above horizon and less than 5 hours from local noon."
#sys.exit(0)
#############################################

def tilted_xy (lat, az, el, el_cutoff):
# Input: latitude of location, and azimuth and elevation of an object.
# Output: xy-coordinates of the gnomon's shadow in a 
# coordinate system where x-axis is towards the east, and y-axis is towards
# the north celestial pole.
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
    theta_shadow = thetap
    rr = gnomon_sz * np.tan(theta_shadow)
    xx = rr * np.cos(phi_shadow)
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
        tmpx, tmpy = tilted_xy(np.deg2rad(mylat), 
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



xpgsz = 7.532894892202737
ypgsz = 4.831104907349666
print('\n PDF page size: ',xpgsz,' x ',ypgsz,'inches. Ratio= ', xpgsz/ypgsz)

opfile = 'alfred_tilted_sundial.pdf'




plt.rcParams["font.family"] = "cursive"
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


# Do dakshinayan
strt_day = s_sol
stop_day = w_sol2
print('\n Making Dakshinayan calculations')

# Offset dakshinayan plots downward by this amount
yoff_dk = 1.2*(ymax - ymin) 

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

xrange = xmax - xmin
yrange = ymax - ymin
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
msg_txt = 'Use the top dial from winter solstice to summer solstice.\n'
msg_txt = msg_txt + 'Use bottom dial from summer solstice to winter solstice.\n'
msg_txt = msg_txt + 'Subtract one hour if dayligt savings is NOT in effect.'
message_xpos = 0.5*(xlo+xhi)
message_ypos = (ylo + yhi)/2 + 0.07*yrange
ax.text(message_xpos, message_ypos, msg_txt, props)


loc_print_xpos = message_xpos
loc_print_ypos = ylo + 0.05*yrange
ax.text(loc_print_xpos, loc_print_ypos, loc_print, props)

pdf_pages = PdfPages(opfile)

fig.subplots_adjust(left = 0., bottom = 0., right = 1., top = 1.)
pdf_pages.savefig(fig)
pdf_pages.close()

#sys.exit(0)
