#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Planar sundial when Sun's elevation is greater than 'elvcutoff_deg'. 
# Assumes a gnomon of height unity. Scale all results appropriately
# if gnonom is not of unit height.
# Input Hour range in local time
#
# TODO: Non-integer hour time zones.

#############################################
# Address: Lat, Lon, Elev (km), Timezone (HH)
# 5544 MapleC: +42.2371552, -77.8020070, 0.675, -5
# 351 Willow : +42.0106135, -71.2165150, 0.040, -5
mylat, mylon, myele, mytz = +42.2371552, -77.8020070, 0.675, -5


# Elevation cutoff for the Sun
elvcutoff_deg = 15.0

# Dates of winter and summer solstice
w_sol="2021-12-21"
s_sol="2022-06-21"
w_sol2="2022-12-21"
#############################################

#############################################
# Import required libraries
import sys
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
loc_print = "Customized for " +str(mylat)+ deg + lat_flag + ", " 
loc_print = loc_print + str(np.abs(mylon)) + deg + lon_flag + ", "
loc_print = loc_print + "elevation " + str(int(myele*3280.84)) + "', "
loc_print = loc_print + "and when the Sun is higher than "
loc_print = loc_print + str(int(elvcutoff_deg)) +deg +" above horizon."
#print(loc_print)
#sys.exit(0)
#############################################


def AzEl2xy (az, el):
    theta = np.deg2rad(90.0 - az) + np.pi
    if (el <= elvcutoff_deg):
        rr = 0.0
    else:
        rr = 1./np.tan(np.deg2rad(el))        
    xx = rr * np.cos(theta)
    yy = rr * np.sin(theta)    
    return xx,yy


def get_pos(start_day, stop_day, hr, mn):
    mystrt  = start_day + " " + str(hr).zfill(2) + ":" + str(mn).zfill(2)
    mystop  = stop_day  + " " + str(hr).zfill(2) + ":" + str(mn).zfill(2)
    myepoch = {'start': mystrt, 'stop':  mystop, 'step':'1d'}
    mysun   = Horizons(id='Sun', location=myloc, epochs=myepoch, id_type='id')
    eph     = mysun.ephemerides(refraction=True)
    ndays = len(eph['AZ'])
    #myxx, myyy = np.zeros(ndays), np.zeros(ndays)
    myxx, myyy = [], []
    for ii in range(ndays):
        tmpx, tmpy = AzEl2xy(eph['AZ'][ii],  eph['EL'][ii])
        #print(tmpx,tmpy)
        #myxx.append(tmpx)
        #myyy.append(tmpy)
    
        if ((tmpx != 0.0) and (tmpy !=0)) :
            myxx.append(tmpx)
            myyy.append(tmpy)
            
    # To ensure the arrays have finite length, add 0 at the end
    myxx.append(0)
    myyy.append(0)
        
    return np.array(myxx), np.array(myyy)

def get_hr_label_pos (x,y):
    rsq = x*x + y*y
    rmax = np.amax(rsq)
    where_max = np.argmax(rsq)
    xmax = x[where_max]
    ymax = y[where_max]
    theta_max = np.arctan2(ymax, xmax)
    xlab = 1.05*np.sqrt(rmax)*np.cos(theta_max)
    ylab = 1.05*np.sqrt(rmax)*np.sin(theta_max)
    lab_rot = theta_max - np.pi/2.0
    return xlab, ylab, np.rad2deg(lab_rot)

def draw_EWNS_unit_line (xpos, ypos):
    plt.plot([xpos-0.5,xpos+0.5], [ypos,ypos], 'k-')
    plt.plot([xpos,xpos], [ypos-0.05,ypos+0.05], 'k-')



# First pass over the whole year to get axis ranges.
labpos = [[0.0, 0.0]]
for hh in range(6, 18 + 1):
    print(' ',hh, end='')
    xx,yy = get_pos(w_sol, w_sol2, hh - mytz, 00)
    xlab,ylab,lab_rot = get_hr_label_pos(xx, yy)
    if (xlab != 0.0 and ylab != 0.0):
        labpos = np.append(labpos, [[xlab, ylab]], axis=0)

# Get axis limits
llim = np.amin(labpos, axis=0)
ulim = np.amax(labpos, axis=0)
# Override limits if need be
#llim = [-2.9, -0.4]
#ulim = [+2.9, +2.6]

xrange = ulim[0] - llim[0]
yrange = ulim[1] - llim[1]
xhi = ulim[0] + 0.02*xrange
xlo = llim[0] - 0.02*xrange
yhi = ulim[1] + 0.15*yrange
ylo = llim[1] - 0.15*yrange
xpgsz = xhi - xlo
ypgsz = yhi - ylo
message_xpos = 0.5*(xlo+xhi)
message_ypos = yhi - 0.05*yrange
loc_print_xpos = message_xpos
loc_print_ypos = ylo + 0.05*yrange
print('\n PDF page size: ',xpgsz,' x ',ypgsz,'inches.')


# Do uttarayan
strt_day = w_sol
stop_day = s_sol

if strt_day == s_sol :
    print('\n Making Dakshinayan calculations')
    opfile = 'dakshinayan.pdf'
    message_text = 'Use this dial from summer solstice to winter solstice. Add one hour if daylight savings is in effect.'
else:
    print('\n Making Uttarayan calculations')
    opfile = 'uttarayan.pdf'
    message_text = 'Use this dial from winter solstice to summer solstice. Add one hour if daylight savings is in effect.'


plt.rcParams["font.family"] = "cursive"
fig = plt.figure(figsize=(xpgsz, ypgsz), dpi=600)

ax = fig.add_subplot(111)
props = {'ha': 'center', 'va': 'center'}
ax.set_xlim([xlo, xhi])
ax.set_ylim([ylo, yhi])
ax.set_aspect('equal', adjustable='box')
pdf_pages = PdfPages(opfile)

skipquarters = 0
labpos = [[0.0, 0.0]]
for hh in range(6, 18 + 1):
    print(' ',hh, end='')
    xx,yy = get_pos(strt_day, stop_day, hh-mytz, 00)
    if (len(xx)>1):
        ax.plot(xx[:-1],yy[:-1],'k-', linewidth=2)
    xlab,ylab,lab_rot = get_hr_label_pos(xx, yy)
    if (xlab != 0.0 and ylab != 0.0):
        ax.text(xlab, ylab, hh, props, rotation=lab_rot)
        labpos = np.append(labpos, [[xlab, ylab]], axis=0)
    if skipquarters == 0:
        for mm in range(15, 46, 15):
            xx,yy = get_pos(strt_day, stop_day, hh-mytz, mm)
            if (len(xx)>1):
                ax.plot(xx[:-1],yy[:-1],'k-', linewidth=0.5)


ax.text(message_xpos, message_ypos, message_text, props)
ax.text(loc_print_xpos, loc_print_ypos, loc_print, props)
draw_EWNS_unit_line(0, 0)

fig.subplots_adjust(left = 0., bottom = 0., right = 1., top = 1.)
pdf_pages.savefig(fig)
pdf_pages.close()

#sys.exit(0)

# Do dakshinayan
strt_day = s_sol
stop_day = w_sol2

if strt_day == s_sol :
    print('\n Making Dakshinayan calculations')
    opfile = 'dakshinayan.pdf'
    message_text = 'Use this dial from summer solstice to winter solstice. Add one hour if daylight savings is in effect.'
else:
    print('\n Making Uttarayan calculations')
    opfile = 'uttarayan.pdf'
    message_text = 'Use this dial from winter solstice to summer solstice. Add one hour if daylight savings is in effect.'


plt.rcParams["font.family"] = "cursive"
fig = plt.figure(figsize=(xpgsz, ypgsz), dpi=600)

ax = fig.add_subplot(111)
props = {'ha': 'center', 'va': 'center'}
ax.set_xlim([xlo, xhi])
ax.set_ylim([ylo, yhi])
ax.set_aspect('equal', adjustable='box')
pdf_pages = PdfPages(opfile)

skipquarters = 0
labpos = [[0.0, 0.0]]
for hh in range(6, 18 + 1):
    print(' ',hh, end='')
    xx,yy = get_pos(strt_day, stop_day, hh-mytz, 00)
    if (len(xx)>1):
        ax.plot(xx[:-1],yy[:-1],'k-', linewidth=2)
    xlab,ylab,lab_rot = get_hr_label_pos(xx, yy)
    if (xlab != 0.0 and ylab != 0.0):
        ax.text(xlab, ylab, hh, props, rotation=lab_rot)
        labpos = np.append(labpos, [[xlab, ylab]], axis=0)
    if skipquarters == 0:
        for mm in range(15, 46, 15):
            xx,yy = get_pos(strt_day, stop_day, hh-mytz, mm)
            if (len(xx)>1):
                ax.plot(xx[:-1],yy[:-1],'k-', linewidth=0.5)


ax.text(message_xpos, message_ypos, message_text, props)
ax.text(loc_print_xpos, loc_print_ypos, loc_print, props)
draw_EWNS_unit_line(0, 0)

fig.subplots_adjust(left = 0., bottom = 0., right = 1., top = 1.)
pdf_pages.savefig(fig)
pdf_pages.close()



sys.exit(0)
