#!/usr/bin/env python3
# Usage: pipe output to gnuplot ("python snapshot.py | gnuplot"), or save output as file and load it in gnuplot.

import sys; out = sys.stdout # alias
import numpy as np

# Prepare python and gnuplot
metadata=open("datadir/metadata.dat").read().replace(r'"""', '"').replace("false", "0").replace("true", "1")
out.write(metadata) # set parameter values to gnuplot
exec(metadata) # set parameter values to python

##### choose #####
plot_params = r"""
frame = 0 # see datadir files suffix for maximum value; 0=initial condition
max_frame_number = 400
"""
exec(plot_params)
out.write(plot_params)

##### set up gnuplot #####
out.write(r"""
set term pdfcairo size 8.4 cm, 5.6 cm font "Times New Roman,11" fontscale 0.7
set output "trajectories.pdf"
set encoding utf8
set minussign

set xrange [0:t_end]
set yrange [0:s_max]

set xlabel '𝑡' offset 0,1
set ylabel '𝑠_𝑖(𝑡)' offset 1,0
set xtics out scale 0.5 offset 0,0.5
set ytics out scale 0.5 offset 0.5,0
set cbtics scale 0.5

set palette rgbformulae 22,13,10 # blue (<0) - green (=0) - red (>0)
# set palette negative # invert
red(x) = int(255*(3*x-1 < 0 ? 0 : (3*x-1 < 1 ? 3*x-1 : 1)))
green(x) = int(255*(sin(pi*x)))
blue(x) = int(255*(cos(pi*x/2.0)))
rgb(x) = 65536 * red(x) + 256 * green(x) + blue(x)

set cbrange [-0.3:0.3]
set cblabel '𝜎_{𝜏𝜏}/𝐸'
set pm3d border retrace # avoid anti-aliasing artefacts (gnuplot v5.5)

# data: y1, y2, colour
plot \
""")


##### plot #####
# plot stress as colour between node n and n-1.
# PBC assumed between 0 and M-1: 2 regions to colour instead of 1, so add 1 extra
# Data format for boxxyerror: x y xmin xmax ymin ymax color (x and y not really used)
# closed_curve: M+1 regions = M-1 inner regions + [0,s_0] + [s_{M-1},s_max]; open_curve: M regions (springs)
M_rect = M if open_curve else M+1
for n in range(M_rect):
    out.write(r"""'-' using ($0*dt_frame):1:($0*dt_frame):(($0+1)*dt_frame):1:2:(-$3) with boxxyerror lc palette notitle,\
""")
    
# plot node arclenth
for n in range(M_nodes):
    if n % m == 0:
        out.write(r"""'-' using ($0*dt_frame):1 with lines lc 'black' lw 2 notitle,\
""")
    else:
        out.write(r"""'-' using ($0*dt_frame):1 with lines lc 'gray' lw 1 notitle,\
""")
out.write(r"""1/0 notitle
""")

# read arclength and stress data for each node
s = np.zeros([M_nodes,max_frame_number+1])
stress = np.zeros([M+1,max_frame_number+1])
for f in range(max_frame_number+1):
    s[:,f] = np.genfromtxt(f"datadir/springs_{f}.dat", usecols=1, max_rows=M_nodes) # (closed_curve: M_nodes=M, skip last repeated row
    stress[:,f] = np.genfromtxt(f"datadir/springs_{f}.dat", usecols=7)


# feed stress data between node n and n-1:
for n in range(1,M_nodes):
    for f in range(max_frame_number+1):
        out.write("%s\t%s\t%s\n" % (s[n-1,f], s[n,f], stress[n,f]))
    out.write("e\n")

if closed_curve:
    # feed stress data between node 0 and bottom x axis
    for f in range(max_frame_number+1):
        out.write("%s\t%s\t%s\n" % (0, s[0,f], stress[0,f]))
    out.write("e\n")
    # feed stress data  between node M-1 and top x axis
    for f in range(max_frame_number+1):
        out.write("%s\t%s\t%s\n" % (s[M_nodes-1,f], s_max, stress[0,f]))
    out.write("e\n")

# feed node arclength data
for n in range(M_nodes):
    for s_ in s[n,:]:
        out.write("%s\n" % s_)
    out.write("e\n")
