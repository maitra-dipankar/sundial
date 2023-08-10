#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 16:03:47 2023

@author: DM

TODO: 
"""

#############################################
# Address: Lat, Lon, Elev (km), Timezone (HH)
# 5544 MapleC: +42.2371552, -77.8020070, 0.675, -5 -4Daylight
# 351 Willow : +42.0106135, -71.2165150, 0.040, -5 -4Daylight
# SCMARS:      +41.965747,  -71.186538,  0.026, -5 -4Daylight
#mylat, mylon, myele, mytz = +41.965747,  -71.186538, 0.026, -4  # SCMARS
mylat, mylon, myele, mytz = +42.237155,  -77.802007, 0.675, -4  # 5544 Mapl
#mylat, mylon, myele, mytz = 0, 0, 0, 0


#############################################

#############################################
# Import required libraries
import numpy as np
from astroquery.jplhorizons import Horizons
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#############################################

#############################################
# Globals
myloc = {'lon': mylon, 'lat': mylat, 'elevation': myele}
lon_flag = " E"
lat_flag = " N"
if mylon < 0 :
    lon_flag = " W"
if mylat < 0 :
    lat_flag = " S"
deg = u"\N{DEGREE SIGN}"
loc_print = 'Subtract one hour if daylight savings is NOT in effect.'
loc_print = "Customized for:\n" +str(mylat)+ deg + lat_flag + ", " 
loc_print = loc_print + str(np.abs(mylon)) + deg + lon_flag + ",\n"
loc_print = loc_print + "Elevation: " + str(int(myele*1000)) + "m, "
loc_print = loc_print + "Time zone: UTC" + f'{mytz:+g}'


# Dates of winter and summer solstice
s_sol="2023-06-21"
w_sol="2023-12-21"
s_sol2="2024-06-21"

# Testing
DPI = 300
pgWidth, pgAngDeg = 6, 30       # page width in inches
pa = np.deg2rad(pgAngDeg)
tipz = pgWidth * np.sin(pa)     # gnomon tip's z-coordinate
axTilt = 23.44

# Output graphics properties
opfmt = 'svg'

font = {'family' : 'cursive',
        'weight' : 'bold',
        'size'   : int(tipz*6)}
matplotlib.rc('font', **font)
props = {'ha': 'center', 'va': 'center'}

#############################################

def solarDirectionCosines (azimuth_deg, altitude_deg):
    '''
    Given azimuth and altitude of the Sun, returns direction cosines of a
    unit vector pointing towards the Sun. The coordinate system is chosen 
    such that the +x-axis points East, +y-axis points North, and the +z-axis 
    points towards zenith.
    
    The relation between the altitude and the polar angle (zenith distance) 
    theta in this case is theta = 90 - altitude
    
    Similarly the relation between azimuth and phi is phi = 90 - azimuth
    
    Since a unit vector in spherical coordinates (r=1, theta, phi) has 
    direction cosines given by
        alpha = sin(theta) cos(phi)
        beta  = sin(theta) sin(phi)
        gamma = cos(theta)
        
    In terms of azimuth and altitude, the unit vector is given by
        alpha = cos(alt) sin(azm)
        beta  = cos(alt) cos(azm)
        gamma = sin(alt)

    Parameters
    ----------
    azimuth_deg : float
        Solar azimuth, in degrees.
    altitude_deg : float
        Solar altitude, in degrees.

    Returns
    -------
    alpha, beta, gamma : float
        Direction cosines such that 
        \vec u = alpha i-hat + beta j-hat + gamma k-hat 
        points towards the Sun.

    '''
    
    azm, alt = np.deg2rad(azimuth_deg), np.deg2rad(altitude_deg)
    cosAlt = np.cos(alt)
    
    alpha = cosAlt * np.sin(azm)
    beta  = cosAlt * np.cos(azm)
    gamma = np.sin(alt)
    
    return alpha, beta, gamma


def shadowLoc (azimuth_deg, altitude_deg, page_width, page_angle_deg):
    '''
    Given solar position (azimuth and altitude) and the book sundial's 
    properties (page width and page angle), find out where the shadow of
    the gnomon's tip will land. This is done as follows:
    
    Our initial coordinate system is one where the +x-axis points towards
    East, +y towards North, and +z points towards the zenith. The 'book' 
    sundial is placed such that the spine of the book is aligned along the
    y-axis, with the gnomon's base at (0,0,0). The width of each left and 
    right page of the book is 'page_width', and each of these pages are 
    tilted page_angle_deg above the flat horizontal.
    
    First we find the z-coordinate of the gnomon's tip so that the tip is
    flush with the left and right edge of the book. This ensures we can 
    catch the shadow of the Sun all the way from sunrise to sunset.
    
    Then we use the properties of the page to compute the components of the
    unit vector normal to the left or right page. Assume the sundial reader is
    on the south side of this book sundial. In this case, the shadow of the 
    gnomon falls on the left page if the Sun is in the eastern sky (solar 
    azimuth < 180 degrees), and on the right page if the Sun is in the 
    western sky.
    Let nx, ny, nz be the components of the unit vector perpendicular to the 
    plane on which the shadow of the gnomon falls. This plane passes thru 
    (0,0,0) and contains the y-axis. Therefore ny=0, and the equation of this
    plane is given by
       [(x,y,z) - (0,0,0) ] dot (nx,ny,nz) = 0
    => nx * x + nz * z = 0        -------------------- Eq. (1)

    
    Based on the solar position, we now compute the direction cosines of the 
    line that creates the shadow of the gnomon's tip. Let alfa, beta, gama be
    these direction cosines.
    
    The line forming the shadow of the gnomon's tip goes thru (0,0,tz) and
    has direction cosines alfa, beta, gama. Therefore the parametric equation 
    of this line is given by
    (x,y,z) = (0,0,tz) + p * (alfa, beta, gamma) ----- Eq. (2)
    
    Here p is a scalar parameter.
    
    Eq. (2) implies that x = p*alfa, and z = tz + p*gama
    Eq. (1) implies that z = -nx*x/nz
    
    Combining: -nx*(p*alfa)/nz = tz + p*gama
            => -p*(gama + alfa*nx/nz) = tz
            => p = -tz/( (gama + alfa*nx/nz) )
    
    Once p is computed, the intersection of the tip's shadow line with the
    page's plane is given by the point
    (xs, ys, zs) = (p*alfa, p*beta, tz + p*gama)      [See eq. (2)]
    
    Now to draw the shadow positions on the dial's flat pages we define a
    2D xy coordinate system for the book. The y-axis of this 'book' 
    coordinate is the same as our original xyz coordinate system with the 
    +y-axis horizontal and pointing towards the north. The origin of the book
    coordinate system is at the same location as that of the original xyz
    system too. The +x-axis however points 'page angle' degrees above due 
    east and the -x-axis points 'page angle' degrees above due west. Thus the
    left page of the book (as seen by someone on the south side of the book) 
    has -ve x coordinates and the right page of the book has +ve x values.
    If the shadow is at (xs,ys,zs) in our original coordinate system, then 
    the book coordinates of this point are (xb,yb) where 
    xb = sign(xs)*sqrt(xs^2 + zs^2) and yb = ys.
    
    
    Parameters
    ----------
    azimuth_deg : float
        Solar azimuth, in degrees.
    altitude_deg : float
        Solar altitude, in degrees.
    page_width : float
        Width of the page (scalable to any unit of length).
    page_angle_deg : float
        Angle that each page makes above the horizontal.

    Returns
    -------
    Coordinates of the shadow on the dial, in 'book' coordinates (xb, yb), 
    in the same units as that of the page_width.
    '''
    
    pa = np.deg2rad(page_angle_deg)
    tipz = page_width * np.sin(pa)            # gnomon tip's z-coordinate
    
    # components of unit vector normal of plane of page
    nx, nz = np.sin(pa), np.cos(pa)  # if sun in E sky (az<=180)
    
    if azimuth_deg > 180 :
        nx = -nx
      
    # direction cosines of the line of shadow's tip
    alfa, beta, gama = solarDirectionCosines(azimuth_deg, altitude_deg)
    
    p = -tipz/( (gama + alfa*nx/nz) )
    
    # intersection of tip's shadow line with page's plane
    xs, ys, zs = p*alfa, p*beta, tipz + p*gama
    
    # convert to book coordinates
    xb, yb = np.sign(xs)*np.sqrt(xs*xs + zs*zs), ys
    
    return xb, yb


def get_pos(start_day, stop_day, hr, mn):
    mystrt  = start_day + " " + str(hr).zfill(2) + ":" + str(mn).zfill(2)
    mystop  = stop_day  + " " + str(hr).zfill(2) + ":" + str(mn).zfill(2)
    myepoch = {'start': mystrt, 'stop':  mystop, 'step':'1d'}
    mysun   = Horizons(id='Sun', location=myloc, epochs=myepoch)
    eph     = mysun.ephemerides(refraction=True)
    ndays = len(eph['AZ'])
    myxx, myyy = [], []
    for ii in range(ndays):
        myaz, myel = eph['AZ'][ii], eph['EL'][ii]
        if myel>0 :
            tmpx, tmpy = shadowLoc (myaz, myel, pgWidth, pgAngDeg)
            myxx.append(tmpx)
            myyy.append(tmpy)
            
    return np.array(myxx), np.array(myyy)


def get_hr_label_pos (x,y):
    ymax = np.amin(y)
    ylab = ymax - 0.2
    xlab = x[np.argmin(y)]
    return xlab, ylab



# Set up page dimensions
xmax = pgWidth
ymax = tipz*np.tan(np.deg2rad(mylat + axTilt))
f = 0.05          # White padding fraction around page
xlo, xhi = -(1+f)*xmax, (1+f)*xmax
ylo, yhi = -(1+f)*ymax, (1+f)*ymax
pgXsz, pgYsz = xhi-xlo, yhi-ylo

print('\n Page size (inches): ',pgXsz,' x ',pgYsz)


# Do dakshinayan
strt_day = s_sol
stop_day = w_sol
print('\n Making Dakshinayan calculations')
opname = 'sundial_book_dakshinayan'
opfile = opname + '.' + opfmt

# Get ready for plotting
fig = plt.figure(figsize=(pgXsz, pgYsz), dpi=DPI, clear=True)
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')
ax.set_xlim([xlo, xhi])
ax.set_ylim([ylo, yhi])

# Draw the page boundaries and gnomon location
ax.plot([-xmax,-xmax], [-ymax, ymax],'k-', lw=4)
ax.plot([-xmax, xmax], [ ymax, ymax],'k-', lw=4)
ax.plot([ xmax, xmax], [ ymax,-ymax],'k-', lw=4)
ax.plot([ xmax,-xmax], [-ymax,-ymax],'k-', lw=4)
ax.plot([0, 0], [-ymax, ymax],'k-', lw=4)
ax.plot([-0.01*xmax,0.01*xmax], [0,0],'k-', lw=2)

skipquarters = 0

for hh in range(5, 20 + 1):
    xx,yy = get_pos(strt_day, stop_day, (hh-mytz)%24, 00)
    if len(xx)>1 : 
        ax.plot(xx,yy,'k-',lw=2)
        xlab,ylab = get_hr_label_pos(xx, yy)
        ax.text(xlab, ylab, hh, props)

    if skipquarters == 0:
        for mm in range(15, 46, 15):
            xx,yy = get_pos(strt_day, stop_day, (hh-mytz)%24, mm)
            if (len(xx)>1):
                ax.plot(xx,yy,'k--', lw=0.5)

# Message which dial to use at when
msg_txt = 'Use this dial from summer solstice to winter solstice.\n'
msg_txt = msg_txt + 'Subtract one hour if dayligt savings is NOT in effect.'
message_xpos = xlo + 0.25*(xhi-xlo)
message_ypos = ylo + 0.15*(yhi-ylo)
ax.text(message_xpos, message_ypos, msg_txt, props)

loc_print_xpos = xlo + 0.75*(xhi-xlo)
loc_print_ypos = message_ypos
ax.text(loc_print_xpos, loc_print_ypos, loc_print, props)
ax.axis('off')
fig.subplots_adjust(left = 0., bottom = 0., right = 1., top = 1.)

# Finish plotting and write output graphics file
if opfmt == 'svg':
    fig.savefig(opfile, format=opfmt, dpi=DPI)
elif opfmt == 'pdf':
    pdf_pages = PdfPages(opfile)
    pdf_pages.savefig(fig)
    pdf_pages.close()



# Do uttarayan
strt_day = w_sol
stop_day = s_sol2
print('\n Making Uttarayan calculations')
opname = 'sundial_book_uttarayan'
opfile = opname + '.' + opfmt

# Get ready for plotting
fig = plt.figure(figsize=(pgXsz, pgYsz), dpi=DPI, clear=True)
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')
ax.set_xlim([xlo, xhi])
ax.set_ylim([ylo, yhi])

# Draw the page boundaries and gnomon location
ax.plot([-xmax,-xmax], [-ymax, ymax],'k-', lw=4)
ax.plot([-xmax, xmax], [ ymax, ymax],'k-', lw=4)
ax.plot([ xmax, xmax], [ ymax,-ymax],'k-', lw=4)
ax.plot([ xmax,-xmax], [-ymax,-ymax],'k-', lw=4)
ax.plot([0, 0], [-ymax, ymax],'k-', lw=4)
ax.plot([-0.01*xmax,0.01*xmax], [0,0],'k-', lw=2)

skipquarters = 0

for hh in range(5, 20 + 1):
    xx,yy = get_pos(strt_day, stop_day, (hh-mytz)%24, 00)
    if len(xx)>1 : 
        ax.plot(xx,yy,'k-', lw=2)
        xlab,ylab = get_hr_label_pos(xx, yy)
        ax.text(xlab, ylab, hh, props)

    if skipquarters == 0:
        for mm in range(15, 46, 15):
            xx,yy = get_pos(strt_day, stop_day, (hh-mytz)%24, mm)
            if (len(xx)>1):
                ax.plot(xx,yy,'k--', lw=0.5)

# Message which dial to use at when
msg_txt = 'Use this dial from winter solstice to summer solstice.\n'
msg_txt = msg_txt + 'Subtract one hour if dayligt savings is NOT in effect.'
message_xpos = xlo + 0.25*(xhi-xlo)
message_ypos = ylo + 0.15*(yhi-ylo)
ax.text(message_xpos, message_ypos, msg_txt, props)

loc_print_xpos = xlo + 0.75*(xhi-xlo)
loc_print_ypos = message_ypos
ax.text(loc_print_xpos, loc_print_ypos, loc_print, props)
ax.axis('off')
fig.subplots_adjust(left = 0., bottom = 0., right = 1., top = 1.)


# Finish plotting and write output graphics file
if opfmt == 'svg':
    fig.savefig(opfile, format=opfmt, dpi=DPI)
elif opfmt == 'pdf':
    pdf_pages = PdfPages(opfile)
    pdf_pages.savefig(fig)
    pdf_pages.close()



# Now make the book stand
from matplotlib import patches

print('\n Making the stand')
opname = 'sundial_book_stand'
opfile = opname + '.' + opfmt

spineH = 1                 # The spine of the book is spineHt above ground
standW = 0.8*pgWidth       # The stand's width under each page
x1, y1 = 0, 0              # Coordinates of the bottom-left corner of stand
x2, y2 = 2*standW, 0       # Coordinates of the bottom-right corner
x3, y3 = x2, spineH + standW*np.sin(pa)   # Top-right corner
x4, y4 = standW, spineH                   # Middle of notch
x5, y5 = 0, y3                            # Top-left corner

standX = np.array([x1, x2, x3, x4, x5, x1])
standY = np.array([y1, y2, y3, y4, y5, y1])

f = 0.05                   # Padding to make the PDF slightly larger
xlo, xhi = -f*x2, (1+f)*x2
ylo, yhi = -f*y3, (1+f)*y3
pgXsz, pgYsz = xhi-xlo, yhi-ylo

fig = plt.figure(figsize=(pgXsz, pgYsz), dpi=DPI, clear=True)
ax = fig.add_subplot(111)
ax.set_aspect('equal', adjustable='box')
ax.set_xlim([xlo, xhi])
ax.set_ylim([ylo, yhi])

ax.plot(standX, standY, 'k-', lw=4)
hole1 = patches.Circle((0.25*standW, 0.4*y3), radius=0.25/2, color='black', 
                        fill=False, lw=4)
hole2 = patches.Circle((1.75*standW, 0.4*y3), radius=0.25/2, color='black', 
                        fill=False, lw=4)
ax.add_patch(hole1)
ax.add_patch(hole2)
ax.axis('off')
fig.subplots_adjust(left = 0., bottom = 0., right = 1., top = 1.)


# Finish plotting and write output graphics file
if opfmt == 'svg':
    fig.savefig(opfile, format=opfmt, dpi=DPI)
elif opfmt == 'pdf':
    pdf_pages = PdfPages(opfile)
    pdf_pages.savefig(fig)
    pdf_pages.close()

