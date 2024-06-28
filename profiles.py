#!/usr/bin/env python3
# Usage: pipe output to gnuplot ("python snapshot.py | gnuplot"), or save output as file and load it in gnuplot.

import sys; out = sys.stdout # alias
import numpy as np

# Prepare python and gnuplot
metadata=open("datadir/metadata.dat").read().replace(r'"""', '"').replace("false", "0").replace("true", "1")
out.write(metadata) # set parameter values to gnuplot
exec(metadata) # set parameter values to python

metadata_cont=open("datadir-continuum/metadata.dat").read().replace(r'"""', '"').replace("false", "0").replace("true", "1")

##### choose #####
plot_params = r"""
frame = 20 # see datadir files suffix for maximum value; 0=initial condition

"""
exec(plot_params)
out.write(plot_params)

##### set up gnuplot #####
out.write(r"""
lw_scale = 0.7
W=8
H=4.5

# restoring forces
if (hookean_force) {
    f(l) = k*(l-a)
} else if (lineardiff_force) {
    f(l) = k*a**2*(1.0/a - 1.0/l)
} else if (porousdiff_force) {
    f(l) = 0.5*k*a**3*(1.0/a**2 - 1.0/l**2)
}

set term pdfcairo size 2*W cm, H cm color font "Times New Roman,11" fontscale 0.7
set encoding utf8
set minussign
set key reverse Left samplen 2.5 # at graph 0,0.68 left center

set output 'profiles.pdf'
set multiplot layout 1,2 margins char 6.5,0.5,2.5,0.2 spacing char 10
""")


##### plot #####
out.write(r"""
# time label
set label 1 '𝑡 = '.sprintf('%4.2f', frame*dt_frame) at graph 0,1 left front textcolor rgb 'black' offset 1,-1

# density
set xrange [0:s_max]
set yrange [0.8:2.6]
set ytics 0.5
set xlabel '𝑠' offset 0,0.5
set ylabel '𝑞(𝑠, 𝑡)' offset 1,0
plot '< sed -n ''2,2p;3q'' """ + f"datadir/springs_{frame}.dat" + r""" | awk ''BEGIN { FS=OFS="\t" } { $2 = ''0.0'' } 1'' """ + " && head -n -1 " + f"datadir/springs_{frame}.dat" + " && tail -n 1 " + f"datadir/springs_{frame}.dat" + r""" | awk ''BEGIN { FS=OFS="\t" } { $2 = ''""" + "%s" % s_max + r"""'' } 1'' ' using 2:(1/(m*$7)) with fsteps lw 1.5 title """ + f"""'𝑚 = {m}',\
"""
+ f""" 'datadir-continuum/data_{frame}.dat' using 2:3 with lines lc 'black' lw 1.5 title 'cont.'
"""
+ r"""

# stress
set xrange [0:s_max]
set yrange [-0.6:0.6]
set ytics 0.25
set xlabel '𝑠'
set ylabel '𝜎_{𝜏𝜏}(𝑠, 𝑡)/𝐸' offset 1,0
plot '< sed -n ''2,2p;3q'' """ + f"datadir/springs_{frame}.dat" + r""" | awk ''BEGIN { FS=OFS="\t" } { $2 = ''0.0'' } 1'' """ + " && head -n -1 " + f"datadir/springs_{frame}.dat" + " && tail -n 1 " + f"datadir/springs_{frame}.dat" + r""" | awk ''BEGIN { FS=OFS="\t" } { $2 = ''""" + "%s" % s_max + r"""'' } 1'' ' using 2:(-$8) with fsteps lw 1.5 notitle,\
"""
+
f""" 'datadir-continuum/data_{frame}.dat' using 2:(-f(1./(m*$3))/(k*a)) with lines lc 'black' lw 1.5 notitle

unset multiplot
unset output
"""
)
