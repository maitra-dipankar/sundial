#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 16:43:05 2019

@author: w00322012

All numbers in inches.
Page setup based on:
    https://stackoverflow.com/questions/29400116/using-matplotlib-how-can-i-print-something-actual-size
"""
import math
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

mylat_deg=41.965686       # Wheaton, between library and old science center
lat=mylat_deg*math.pi/180.
sinlat=math.sin(lat)
coslat=math.cos(lat)

ply_thick = 23./32.     # Thickness of the plywood
ft_len = 0.9*ply_thick  # How far the foot will protrude into the baseboard
tab_ht = 1.
tab_wd = 1.
hole_dia = 17./64.      # Clearance hole dia for 1/4-20

# Start at th bottom of the south side
xx = [0.5]; yy = [0.5]

# Go to the top of the south side
xx.append(xx[-1]) ; yy.append(yy[-1] + ply_thick + 2.)

# Make the tab
xx.append(xx[-1] - tab_ht*sinlat) ; yy.append(yy[-1] + tab_ht*coslat)
xx.append(xx[-1] + tab_wd*coslat) ; yy.append(yy[-1] + tab_wd*sinlat)
xx.append(xx[-1] + tab_ht*sinlat) ; yy.append(yy[-1] - tab_ht*coslat)

# Go to the other tab at the top of the north side
ww=10. - 2.*tab_wd
xx.append(xx[-1] + ww*coslat) ; yy.append(yy[-1] + ww*sinlat)

# Make the tab
xx.append(xx[-1] - tab_ht*sinlat) ; yy.append(yy[-1] + tab_ht*coslat)
xx.append(xx[-1] + tab_wd*coslat) ; yy.append(yy[-1] + tab_wd*sinlat)
xx.append(xx[-1] + tab_ht*sinlat) ; yy.append(yy[-1] - tab_ht*coslat)

# Save this x-coordinate to make the foot board
n_footboard_end = xx[-1]
s_footboard_end = xx[0]

# Come down to the bottom of the north side and make the north foot
xx.append(xx[-1]) ; yy.append(yy[0])
xx.append(xx[-1] - tab_wd) ; yy.append(yy[-1])
xx.append(xx[-1]) ; yy.append(yy[-1]+ft_len)

# Go to the south side and make the south foot
xx.append(xx[0] + tab_wd) ; yy.append(yy[-1])
xx.append(xx[-1]) ; yy.append(yy[0])
xx.append(xx[0]) ; yy.append(yy[0])

# Put 3 holes
xh1 = xx[0]+tab_wd   ; yh1 = yy[0]+2.
xh2 = max(xx)-tab_wd ; yh2 = yh1
xh3 = xh2            ; yh3 = yh2+(xh2-xh1)*math.tan(lat)
hole1 = plt.Circle((xh1,yh1), hole_dia/2., color='k', fill=False)
hole2 = plt.Circle((xh2,yh2), hole_dia/2., color='k', fill=False)
hole3 = plt.Circle((xh3,yh3), hole_dia/2., color='k', fill=False)

# Set up the page
xpgsz = 1. + max(xx)-min(xx)
ypgsz = 1. + max(yy)-min(yy)

# The PDF document
pdf_pages = PdfPages('stand.pdf')

# Create a figure instance
fig = plt.figure(figsize=(xpgsz, ypgsz), dpi=600)
ax = fig.add_subplot(111)
ax.set_xlim(min(xx)-0.5, max(xx)+0.5)
ax.set_ylim(min(yy)-0.5, max(yy)+0.5)

ax.plot(xx,yy,'k-')
ax.add_artist(hole1)
ax.add_artist(hole2)
ax.add_artist(hole3)

fig.subplots_adjust(left = 0., bottom = 0., right = 1., top = 1.)

# Done with the page and write PDF to disk
pdf_pages.savefig(fig)
pdf_pages.close()

print(xpgsz, ypgsz)

# Now we make the footboard
xholes_end2end = n_footboard_end-s_footboard_end
xlen_footboard = xholes_end2end + 1.
ylen_footboard = 3. + 23./32. + 1.
xpgsz = xlen_footboard + 0.5
ypgsz = ylen_footboard + 0.5
print(xpgsz, ypgsz)

# Set up the page for the PDF document
pdf_pages = PdfPages('footboard.pdf')

# Create a figure instance
fig = plt.figure(figsize=(xpgsz, ypgsz), dpi=600)
ax = fig.add_subplot(111)
ax.set_xlim(0., xpgsz)
ax.set_ylim(0., ypgsz)

bb = mpatches.Rectangle((0.25,0.25), xlen_footboard, ylen_footboard, 
                        fill=False)
h1 = mpatches.Rectangle((0.75,0.75), tab_wd, ply_thick, fill=False)
h2 = mpatches.Rectangle((0.75+xholes_end2end-tab_wd,0.75), tab_wd, ply_thick, 
                        fill=False)
h3 = mpatches.Rectangle((0.75,0.75+3.), tab_wd, ply_thick, fill=False)
h4 = mpatches.Rectangle((0.75+xholes_end2end-tab_wd,0.75+3.), tab_wd, 
                        ply_thick, fill=False)

ax.add_patch(bb)
ax.add_patch(h1)
ax.add_patch(h2)
ax.add_patch(h3)
ax.add_patch(h4)

hh1 = plt.Circle((1.25,0.25+ylen_footboard/2.), hole_dia/2., fill=False)
hh2 = plt.Circle((1.25+xholes_end2end-tab_wd,0.25+ylen_footboard/2.), hole_dia/2., fill=False)

ax.add_artist(hh1)
ax.add_artist(hh2)


fig.subplots_adjust(left = 0., bottom = 0., right = 1., top = 1.)

# Done with the page and write PDF to disk
pdf_pages.savefig(fig)
pdf_pages.close()

