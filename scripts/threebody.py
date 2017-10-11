#! /usr/bin/env python

from math import *
import numpy as np
import time
import argparse
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description="Compute 3body histograms and plot them")

parser.add_argument("-i",
                    action="store", nargs='?', default="threebody.log",
                    required=False, dest="inputfilename",
                    help="the name of input file")
parser.add_argument("-o",
                    action="store", nargs='?', default="histogram.out",
                    required=False, dest="outputfilename",
                    help="the name of output file")
parser.add_argument("-s",
                    action="store",nargs='?',
                    required=True, dest="system_type",
                    help="System type: isotropic or anisotropic")
parser.add_argument("-rl",
                    action="store",nargs='?', type = float,
                    required=False, dest="rl", default= float(3),
                    help="r low")
parser.add_argument("-rh",
                    action="store",nargs='?', type = float,
                    required=False, dest="rh", default= float(6),
                    help="r high")
parser.add_argument("-rs",
                    action="store",nargs='?', type = float,
                    required=False, dest="rstep", default= float(0.05),
                    help="r step")
parser.add_argument("-zl",
                    action="store",nargs='?', type = float,
                    required=False, dest="zl", default= float(3),
                    help="z low")
parser.add_argument("-zh",
                    action="store",nargs='?', type = float,
                    required=False, dest="zh", default= float(6),
                    help="z high")
parser.add_argument("-zs",
                    action="store",nargs='?', type = float,
                    required=False, dest="zstep", default= float(0.05),
                    help="z step")

args = parser.parse_args()
inputfilename = args.inputfilename
outputfilename = args.outputfilename
system_type = args.system_type
rl = args.rl
rh = args.rh
rstep = args.rstep
zl = args.zl
zh = args.zh
zstep = args.zstep

time1 = time.time()

# Read the input file
def build_isotropic_histogram(inputfilename, r_bins, r_values, r_counter):

    with open(inputfilename, 'r') as inputfile:
        print "reading the input file"
        for line in inputfile:
            # ignore the comments
            line = line.split('#', 1)[0]
            line = line.strip()
            if line != "":
                columns = line.split()
                try:
                    rik = float(columns[0])
                    energy = float(columns[4])
                    r_values, r_counter = evolve_1d_histogram(r_bins, r_values, rik, energy, r_counter)
                except ValueError:
                    print "something went wrong, check your input file. Exiting..."
                    exit(1)

    return r_values, r_counter

def build_anisotropic_histogram(inputfilename, r_bins, z_bins, values, counter):

    with open(inputfilename, 'r') as inputfile:
        print "reading the input file"
        for line in inputfile:
            # ignore the comments
            line = line.split('#', 1)[0]
            line = line.strip()
            if line != "":
                columns = line.split()
                try:
                    xik = float(columns[1])
                    yik = float(columns[2])
                    rik = np.sqrt(np.square(xik) + np.square(yik))
                    zik = float(columns[3])
                    energy = float(columns[4])
                    values, counter = evolve_2d_histogram(r_bins, z_bins, values, counter, rik, zik, energy)
                except ValueError:
                    print "something went wrong, check your input file. Exiting..."
                    exit(1)

    return values, counter

def evolve_1d_histogram(bins, values, xvalue, yvalue, counter):
    if xvalue > bins[0] and xvalue < bins[-1]:
        bin = int((xvalue - bins[0])/(bins[1] - bins[0]))
        values[bin] += yvalue
        counter[bin] += 1
    return values, counter

def evolve_2d_histogram(r_bins, z_bins, values, counter, rik, zik, energy):
    if rik > r_bins[0] and rik < r_bins[-1] and zik > z_bins[0] and zik < z_bins[-1]:
        r_bin = int((rik - r_bins[0])/(r_bins[1] - r_bins[0]))
        z_bin = int((zik - z_bins[0])/(z_bins[1] - z_bins[0]))
        values[r_bin, z_bin] += energy
        counter[r_bin, z_bin] += 1
    return values, counter

def write_isotropic_output(outputfilename, r_centers, r_values, r_counter):
    with open(outputfilename, 'w') as outputfile:
        outputfile.write('#binCenter    average_energy    counts\n')
        np.savetxt(outputfile, np.c_[r_centers, r_values, r_counter])

def write_anisotropic_output(outputfilename, r_centers, z_centers, values, counter):
    with open(outputfilename, 'w') as outputfile:
        outputfile.write('#binCenter_r  binCenter_z  average_energy    counts\n')
        for i, z in enumerate(z_centers):
            for j, r in enumerate(r_centers):
                outputfile.write("%5.3f %5.3f %5.8f %i\n" % (z, r, values[i, j], counter[i, j]))

def plot_1d(centers, values, plotname):
    ax = plt.subplot()
    plt.plot(centers, values, linewidth = 1.5)
    plt.xlabel(r'$\mathrm{d_{ij}}$', fontsize = 18)
    plt.ylabel(r'$\mathrm{Three-body\/energy}$', fontsize = 18)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    plt.tick_params(which='both', width=2)

    plt.savefig(plotname, dpi=1000, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='letter', format='png',
        transparent=True, bbox_inches='tight', pad_inches=0.1,
        frameon=None)
    plt.clf()

def plot_2d(r_centers, z_centers, values, plotname):
    ax = plt.subplot()
    plt.pcolormesh(r_centers, z_centers, values)
    plt.colorbar()
    CP = plt.contour(r_centers, z_centers, values, 15, linewidths=2, colors='white', extend='neither')
    plt.clabel(CP, inline=1, fmt='%1.6f', fontsize=14, colors='white')
    plt.xlabel(r'$\mathrm{r_{ij}}$', fontsize = 18)
    plt.ylabel(r'$\mathrm{z_{ij}}$', fontsize = 18)

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    plt.tick_params(which='both', width=2)

    plt.savefig(plotname, dpi=1000, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='letter', format='png',
        transparent=True, bbox_inches='tight', pad_inches=0.1,
        frameon=None)
    plt.clf()

def plot_projection(r_centers, z_centers, n, values, plotname):
    ax = plt.subplot()
    plt.xlabel(r'$\mathrm{r_{ij}}$', fontsize = 18)
    plt.ylabel(r'$\mathrm{Three-body\/energy}$', fontsize = 18)
    jump = int(len(z_centers)/n)
    for i in range(n):
        plt.plot(r_centers, values[i*jump,:], label = r'$\mathrm{z_{ij}}$'+ ' = ' + str(z_centers[i*jump]))
        plt.legend()

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    plt.tick_params(which='both', width=2)

    plt.savefig(plotname, dpi=1000, facecolor='w', edgecolor='w',
        orientation='portrait', papertype='letter', format='png',
        transparent=True, bbox_inches='tight', pad_inches=0.1,
        frameon=None)
    plt.clf()

def run_isotropic(rl, rh, rstep, inputfilename, outputfilename):
    r_bins = np.arange(rl, rh + rstep, rstep)
    r_values = np.zeros(len(r_bins)-1, dtype=float)
    r_counter = np.zeros(len(r_bins)-1)
    r_values, r_counter = build_isotropic_histogram(inputfilename, r_bins, r_values, r_counter)
    r_values = r_values / r_counter
    # In case there is nan in the array, replace it with 0
    r_values = np.nan_to_num(r_values)

    r_centers = (r_bins[:-1] + r_bins[1:])/2
    write_isotropic_output(outputfilename, r_centers, r_values, r_counter)
    return r_values, r_counter, r_centers

def run_anisotropic(rl, rh, rstep, zl, zh, zstep, inputfilename, outputfilename):
    #rstep = 1; zstep = 1
    r_bins = np.arange(rl, rh + rstep, rstep)
    z_bins = np.arange(zl, zh + zstep, zstep)
    values = np.zeros((len(r_bins)-1, len(z_bins)-1), dtype=float)
    counter = np.zeros((len(r_bins)-1, len(z_bins)-1))
    r_centers = (r_bins[:-1] + r_bins[1:])/2
    z_centers = (z_bins[:-1] + z_bins[1:])/2


    values, counter = build_anisotropic_histogram(inputfilename, r_bins, z_bins, values, counter)
    values = values / counter
    # In case there is nan in the array, replace it with 0
    values = np.nan_to_num(values)

    return values, counter, r_centers, z_centers

if system_type == 'isotropic':
    r_values, r_counter, r_centers = run_isotropic(rl, rh, rstep, inputfilename, outputfilename)
    plot_1d(r_centers, r_values, '3body_1d.png')

if system_type == 'anisotropic':
    values, counter, r_centers, z_centers = run_anisotropic(rl, rh, rstep, zl, zh, zstep, inputfilename, outputfilename)
    write_anisotropic_output(outputfilename, r_centers, z_centers, values, counter)
    plot_2d(r_centers, z_centers, values, '3body_2d.png')
    plot_projection(r_centers, z_centers, 5, values, '3body_2d_projection.png')

time2 = time.time()
calculation_time = time2-time1
print "Calculation time = " + str(calculation_time) + "seconds"
