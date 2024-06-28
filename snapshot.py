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

plot_interface = 1
plot_coils = 1
plot_colored_coils = 1
plot_normal_stress_interface = 0 # 1 for colouring normal stress as a coloured offset interface (Model4)
plot_cells = 0 # 1 for plotting a substrate with cell bodies
plot_spring_background = 0
plot_elliptic_spring_nodes = 0
plot_elliptic_cell_nodes = 0

lw_scale = 0.75
interface_lc = 'black'
interface_lw = 1*2*lw_scale
spring_nodes_lc = 'white'
spring_nodes_lw = 2*lw_scale
cell_nodes_lc = 'black'
cell_nodes_lw = 1.5*lw_scale
cell_halfheight = 0.17
spring_background_lw = 1
spring_coil_lw = 1.5*lw_scale
spring_coil_lc = 'grey50' # unless plot_colored_coils: -> stress-based
cellnodesize = 0.04
springnodesize = 0.04 #*1e-10 # to make invisible

"""
exec(plot_params)
out.write(plot_params)

out.write(r"""
if (plot_colored_coils) {spring_coil_lw = 4*lw_scale}

if (plot_coils) {
n = 1; # number of coil periods
lmax = 1.53*s_max/M; # max length of a spring
f(s,l) = sqrt(lmax>l? lmax**2-l**2 : 0)*sin(2*pi*n*s/l)/(4.*n);
coilsamples = n*20;
spring_nodes_lw = 1.*lw_scale;
}
""")

##### set up gnuplot #####
out.write(r"""
set size ratio -1
set style fill solid border lc 'black'
set xrange [-1.15*R:1.15*R]
set yrange [-1.15*R:1.15*R]
set tics front
set colorbox
set cbrange [-0.3:0.3]; # used for coil colours only

if (plot_normal_stress_interface) {
# we use another colorbar to plot normal stress (different units)
# Here we set the min and max value of that other colorbar
cb2_min = -0.8;
cb2_max = 0.8;
}

# clamping function, rescales val in [valmin:valmax] to grey in [0:1]
grey(val,valmin,valmax) = val < valmax ? (val > valmin ? (val-valmin)/(valmax-valmin) : 0) : 1

set palette rgbformulae 22,13,10 # blue (<0) - green (=0) - red (>0)
# set palette negative # invert
red(x) = int(255*(3*x-1 < 0 ? 0 : (3*x-1 < 1 ? 3*x-1 : 1)))
green(x) = int(255*(sin(pi*x)))
blue(x) = int(255*(cos(pi*x/2.0)))
rgb(x) = 65536 * red(x) + 256 * green(x) + blue(x)

# interface parametrisations
if (model eq "Model0") {# circle
	x(t)=R*cos(t);
	y(t)=R*sin(t);
	xp(t) = -R*sin(t); # x'(t)
	yp(t) = R*cos(t); # y'(t)
	W=200;
	H=200;
} else if (model[1:7] eq "Model2") {# smooth cross
	r_theta(t)=R*( cos(t)**4 + sin(t)**4 );
	x(t)=r_theta(t)*cos(t);
	y(t)=r_theta(t)*sin(t);
	r_theta_p(t) = 4*R*(-sin(t)*cos(t)**3 + cos(t)*sin(t)**3);
	xp(t) = r_theta_p(t)*cos(t) - r_theta(t)*sin(t);
	yp(t) = r_theta_p(t)*sin(t) + r_theta(t)*cos(t);
	W=200;
	H=200;
} else if (model eq "Model4") { # sine curve
	x(t) = t;
	y(t) = R*sin(t);
	xp(t)=1;
	yp(t)=R*cos(t);
	kappa(x) = -R*sin(x)/(1+(R*cos(x))**2)**1.5 # curvature
	set xrange [-0.05*2*pi:1.05*2*pi];
    if (plot_cells) {
		set yrange [-1.3*R:1.3*R];
		W=320;
		H=110;
	} else {
		set yrange [-1.2*R:1.2*R];
		W=320;
		H=100;
	}
	set ytics 0.4
}

# tangent vector
tau1(t) = xp(t)/sqrt(xp(t)**2+yp(t)**2)
tau2(t) = yp(t)/sqrt(xp(t)**2+yp(t)**2)

    
# plot parametric interface
if (plot_interface) {
	set parametric;
	set table $interface;
	set samples 200;
	plot [0:u_max] x(t), y(t);
	unset table;
	unset parametric;
}

set term pdfcairo font "Times-New-Roman,12" size 0.3*W/10. cm, 0.3*H/10. cm color
set output 'snapshot.pdf'
set encoding utf8
set minussign
set cbtics format '%g'
set format '%g'

# time label
if (model eq "Model4") { # sine curve
    set label 1 '𝑡 = '.sprintf('%4.2f',frame*dt_frame) at graph 0.435,1 left front textcolor rgb 'black' offset 1,-1
} else {
    set label 1 '𝑡 = '.sprintf('%4.2f',frame*dt_frame) at graph 0,1 left front textcolor rgb 'black' offset 1,-1
}
""")


##### plot #####

if straight_springs:
    # plot interface, nodes as circles, straight springs as stress-coloured lines
    if plot_coils:
        # plot cell nodes as black circles, inner nodes as white circles, springs as stress-coloured coils
        utl = np.genfromtxt(f"datadir/springs_{frame}.dat", usecols=[0,2,3,4,5,6,7])
        # utl[i] = [u_i, x_i, y_i, tau1_i, tau2_i, length_i, stress], length_i=s_i-s_{i-1}

        # portion of interface plot data
        for i in range(1,M+1): # last row is same as first in u_tau1_tau2_length
            out.write(r"""
set parametric
set table $coil""" + "%d" % i + r"""
set sample coilsamples
plot """ + "[{0}:{1}] {5}+({6}-({5}))*(t-{0})/({1}-{0})+(f({2}*(t-{0})/({1}-{0}),{2}))*({3}), {7}+({8}-({7}))*(t-{0})/({1}-{0})+(f({2}*(t-{0})/({1}-{0}),{2}))*({4})".format((utl[i-1,0]-u_max) if (i==M and closed_curve) else utl[i-1,0], utl[i,0], utl[i,5], -(utl[i,2]-utl[i-1,2])/np.sqrt((utl[i,1]-utl[i-1,1])**2+(utl[i,2]-utl[i-1,2])**2), (utl[i,1]-utl[i-1,1])/np.sqrt((utl[i,1]-utl[i-1,1])**2+(utl[i,2]-utl[i-1,2])**2),utl[i-1,1],utl[i,1],utl[i-1,2],utl[i,2]) + r""" #{0}=u_{i-1}, {1}=u_i, {2}=length_i, {3}=n1=-tau2, {4}=n2=tau1, {5}=x_{i-1}, {6}=x_i, {7}=y_{i-1}, {8}=y_i
unset table
unset parametric
""")
                
    # plot cell nodes as black circles, inner nodes as white circles, springs as stress-coloured segments
    out.write(r"""
plot \
""")
    if plot_interface:    
        out.write(r"""$interface with lines lc ''.interface_lc lw interface_lw notitle,\
""")
    if plot_spring_background:
        out.write(f"""'datadir/springs_{frame}.dat' using 3:4:($8) with lines lc palette lw spring_background_lw notitle,\
""")
    if plot_coils:
        for i in range(1,M+1):
            if plot_colored_coils:
                out.write("$coil%d with lines lw spring_coil_lw lc palette cb %f" % (i,utl[i,6]) + r""" notitle,\
""")
            else:
                out.write("$coil%d with lines lw spring_coil_lw lc ''.spring_coil_lc" % i + r""" notitle,\
""")
        if plot_elliptic_spring_nodes:
            out.write(f"""'datadir/springs_{frame}.dat' using 3:4:(2*springnodesize):(2*springnodesize):(atan2($6,$5)*180/pi) with ellipses lc ''.spring_nodes_lc lw spring_nodes_lw notitle,\
""")
        else:
            out.write(f"""'datadir/springs_{frame}.dat' using 3:4:(springnodesize) with circles lc ''.spring_nodes_lc lw spring_nodes_lw notitle,\
""")
                
    else:
        out.write(f"""'datadir/springs_{frame}.dat' using 3:4:(springnodesize) with circles lc ''.spring_nodes_lc lw spring_nodes_lw notitle,\
""")
    out.write(f"""'datadir/springs_{frame}.dat' every m using 3:4:(cellnodesize) with circles lc ''.cell_nodes_lc lw cell_nodes_lw notitle
""")

            
else: # curved springs        
    if plot_coils:
        # plot cell nodes as white circles, inner nodes as black circles, springs as stress-coloured coils
        utl = np.genfromtxt(f"datadir/springs_{frame}.dat", usecols=[0,4,5,6])
        # utl[i] = [u_i, tau1_i, tau2_i, length_i], length_i=s_i-s_{i-1}

        # portion of interface plot data
        for i in range(1,M+1): # last row is same as first in u_tau1_tau2_length
            out.write(r"""
set parametric
set table $coil""" + "%d" % i + r"""
set sample coilsamples
plot """ + "[{0}:{1}] x(t)+(f({2}*(t-{0})/({1}-{0}),{2}))*(-tau2(t)), y(t)+(f({2}*(t-{0})/({1}-{0}),{2}))*(tau1(t))".format((utl[i-1,0]-u_max) if (i==M and closed_curve) else utl[i-1,0], utl[i,0], utl[i,3], -0.5*(utl[i,2]+utl[i-1,2]), 0.5*(utl[i,1]+utl[i-1,1])) + r""" #{0}=u_{i-1}, {1}=u_i, {2}=length_i, {3}=n1=-tau2, {4}=n2=tau1
unset table
unset parametric
""")
            
    # plot cell nodes as white circles, inner nodes as black circles, springs as portions of the parametric surface
    s = np.genfromtxt(f"datadir/springs_{frame}.dat", usecols=[0,7])
    kappa = np.genfromtxt(f"datadir/springs_{frame}.dat", usecols=[0,8])

    # portion of interface plot data
    for i in range(1,M+1): # last row is same as first in s for closed_curve
        out.write(r"""
set parametric
set format # reset otherwise with latex, points are outputted wrapped in latex commands!            
set table $interface""" + "%d" % i + r"""
set samples 20

plot """ + "[%f:%f] x(t), y(t)" % ((s[i-1,0]-u_max) if (i==M and closed_curve) else s[i-1,0], s[i,0]) + r"""

unset table
unset parametric
""")

    out.write(r"""
plot \
""")
    if plot_interface:
        if plot_cells:
            out.write(r"""$interface using 1:($2-cell_halfheight) with filledcurves x1 fc "#dddddd" fs solid 1 notitle,\
$interface using 1:($2-cell_halfheight):($2+cell_halfheight) with filledcurves fc "#f2f2f2" fs solid 1 notitle ,\
$interface using 1:($2-cell_halfheight) with lines lc "grey70" lw 1 notitle,\
$interface using 1:($2+cell_halfheight) with lines lc "grey70" lw 1 notitle,\
$interface with lines lc ''.interface_lc lw interface_lw notitle,\
""")
        else:
            out.write(r"""$interface with lines lc ''.interface_lc lw interface_lw notitle,\
""")

    if plot_spring_background:
        for i in range(1,M+1):# gnuplot joins i to i-1 
            out.write("$interface%d with points lw spring_background_lw lc palette cb %f" % (i,-s[i,1]) + r""" notitle,\
""")
    if plot_normal_stress_interface:
        for i in range(1,M+1):# gnuplot joins i to i-1
            out.write("$interface%d using 1:($2-cell_halfheight-0.04):(rgb(grey(%f*kappa($1),cb2_min,cb2_max))) with lines lw 2*spring_coil_lw lc rgb var" % (i,-s[i,1]) + r""" notitle,\
$interface using 1:($2-cell_halfheight) with lines lc "grey70" lw 1 notitle,\
$interface using 1:($2+cell_halfheight) with lines lc "grey70" lw 1 notitle,\
""")            
    if plot_coils:
        for i in range(1,M+1):
            if plot_colored_coils:
                out.write("$coil%d with lines lw spring_coil_lw lc palette cb %f" % (i,-s[i,1]) + r""" notitle,\
""")
            else:
                out.write("$coil%d with lines lw spring_coil_lw lc ''.spring_coil_lc " % i + r""" notitle,\
""")
        if plot_elliptic_spring_nodes:
            out.write(f"""'datadir/springs_{frame}.dat' using 3:4:(2*springnodesize):(2*springnodesize):(atan2($6,$5)*180/pi) with ellipses lc ''.spring_nodes_lc lw spring_nodes_lw notitle,\
""")
        else:
            out.write(f"""'datadir/springs_{frame}.dat' using 3:4:(springnodesize) with circles lc ''.spring_nodes_lc lw spring_nodes_lw notitle,\
""")
    else:                
        out.write(f"""'datadir/springs_{frame}.dat' using 3:4:(springnodesize) with circles lc ''.spring_nodes_lc lw spring_nodes_lw notitle,\
""")
    if plot_elliptic_cell_nodes:        
        out.write(f"""'datadir/springs_{frame}.dat' every m using 3:4:(0.1*cellnodesize):(7*cellnodesize):(atan2($6,$5)*180/pi) with ellipses lc ''.cell_nodes_lc lw cell_nodes_lw notitle,\
""")
    out.write(f"""'datadir/springs_{frame}.dat' every m using 3:4:(cellnodesize) with circles lc ''.cell_nodes_lc lw cell_nodes_lw notitle
""")

out.write(r"""
unset output
""")
