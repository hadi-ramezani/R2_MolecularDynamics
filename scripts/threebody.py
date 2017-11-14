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
parser.add_argument("-thetal",
                    action="store",nargs='?', type = float,
                    required=False, dest="thetal", default= float(0),
                    help="theta low")
parser.add_argument("-thetah",
                    action="store",nargs='?', type = float,
                    required=False, dest="thetah", default= float(180),
                    help="theta high")
parser.add_argument("-thetas",
                    action="store",nargs='?', type = float,
                    required=False, dest="thetastep", default= float(5),
                    help="theta step")
parser.add_argument("-thetaZl",
                    action="store",nargs='?', type = float,
                    required=False, dest="thetaZl", default= float(0),
                    help="thetaZ low")
parser.add_argument("-thetaZh",
                    action="store",nargs='?', type = float,
                    required=False, dest="thetaZh", default= float(90),
                    help="thetaZ high")
parser.add_argument("-thetaZs",
                    action="store",nargs='?', type = float,
                    required=False, dest="thetaZstep", default= float(5),
                    help="thetaZ step")

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
thetal = args.thetal
thetah = args.thetah
thetastep = args.thetastep
thetaZl = args.thetaZl
thetaZh = args.thetaZh
thetaZstep = args.thetaZstep

time1 = time.time()

# Read the input file
def build_isotropic_histogram(inputfilename, r_bins, r_values, r_counter, theta_bins, theta_counter, thetaZ_bins, thetaZ_counter):

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
                    thetakij = float(columns[4])
                    thetazij = float(columns[5])
                    energy = float(columns[6])
                    r_values, r_counter, theta_counter, thetaZ_counter = evolve_1d_histogram(r_bins, r_values, rik, energy, r_counter,
                                                                             theta_bins, theta_counter, thetakij,
                                                                             thetaZ_bins, thetaZ_counter, thetazij)
                except ValueError:
                    print "something went wrong, check your input file. Exiting..."
                    exit(1)

    return r_values, r_counter, theta_counter, thetaZ_counter

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
                    energy = float(columns[6])
                    values, counter = evolve_2d_histogram(r_bins, z_bins, values, counter, rik, zik, energy)
                except ValueError:
                    print "something went wrong, check your input file. Exiting..."
                    exit(1)

    return values, counter

def evolve_1d_histogram(bins, values, xvalue, yvalue, counter, theta_bins, theta_counter, thetakij, thetaZ_bins, thetaZ_counter, thetazij):
    if (xvalue > bins[0] and xvalue < bins[-1]) and (thetakij > theta_bins[0] and thetakij < theta_bins[-1]) and (thetazij > thetaZ_bins[0] and thetazij < thetaZ_bins[-1]):
        bin = int((xvalue - bins[0])/(bins[1] - bins[0]))
        values[bin] += yvalue
        counter[bin] += 1
        theta_bin = int((thetakij - theta_bins[0])/(theta_bins[1] - theta_bins[0]))
        theta_counter[theta_bin, bin] +=1
        thetaZ_bin = int((thetazij - thetaZ_bins[0])/(thetaZ_bins[1] - thetaZ_bins[0]))
        thetaZ_counter[thetaZ_bin] +=1

    return values, counter, theta_counter, thetaZ_counter

def evolve_2d_histogram(r_bins, z_bins, values, counter, rik, zik, energy):
    if (rik > r_bins[0] and rik < r_bins[-1]) and (zik > z_bins[0] and zik < z_bins[-1]):
        r_bin = int((rik - r_bins[0])/(r_bins[1] - r_bins[0]))
        z_bin = int((zik - z_bins[0])/(z_bins[1] - z_bins[0]))
        values[z_bin, r_bin] += energy
        counter[z_bin, r_bin] += 1
    return values, counter

def write_isotropic_output(outputfilename, r_centers, r_values, r_counter):
    with open(outputfilename, 'w') as outputfile:
        outputfile.write('#binCenter    average    counts\n')
        np.savetxt(outputfile, np.c_[r_centers, r_values, r_counter])

def write_anisotropic_output(outputfilename, r_centers, z_centers, values, counter):
    with open(outputfilename, 'w') as outputfile:
        outputfile.write('#binCenter_r  binCenter_z  average_energy    counts\n')
        for j, z in enumerate(z_centers):
            for i, r in enumerate(r_centers):
                outputfile.write("%5.3f %5.3f %5.8f %i\n" % (z, r, values[j, i], counter[j, i]))

def write_theta_output(outputfilename, r_centers, theta_centers, theta_counter):
    with open(outputfilename, 'w') as outputfile:
        outputfile.write('#binCenter_r  binCenter_theta  counts\n')
        for j, theta in enumerate(theta_centers):
            for i, r in enumerate(r_centers):
                outputfile.write("%5.3f %5.3f %i\n" % (theta, r, theta_counter[j, i]))

def write_thetaZ_output(outputfilename, centers, counter):
    with open(outputfilename, 'w') as outputfile:
        outputfile.write('#binCenter    counts\n')
        np.savetxt(outputfile, np.c_[centers, counter])

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

def plot_thetaZ(thetaZ_centers, thetaZ_counter, plotname):
    ax = plt.subplot()
    thetaZ_counter = thetaZ_counter / float(np.sum(thetaZ_counter))

    plt.plot(thetaZ_centers, thetaZ_counter, linewidth = 1.5)
    plt.xlabel(r'$\mathrm{\theta_{zik}}$', fontsize = 18)
    plt.ylabel(r'$\mathrm{P(\theta_{zik})}$', fontsize = 18)
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

def plot_theta(r_centers, theta_centers, theta_counter, plotname):
    ax = plt.subplot()
    #convert counts to probability. This gives the probability of theta at any r. Sum of the probabilities at any r
    # adds to 1
    for col in range(theta_counter.shape[1]):
        theta_counter[:,col] = theta_counter[:,col]/float(np.sum(theta_counter[:,col]))

    plt.pcolormesh(r_centers, theta_centers, theta_counter)
    plt.colorbar()
    #CP = plt.contour(r_centers, theta_centers, theta_counter, 15, linewidths=2, colors='white', extend='neither')
    #plt.clabel(CP, inline=1, fmt='%1.6f', fontsize=14, colors='white')
    plt.xlabel(r'$\mathrm{d_{ij}}$', fontsize = 18)
    plt.ylabel(r'$\mathrm{\theta_{kij}}$', fontsize = 18)

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

def run_isotropic(rl, rh, rstep, thetal, thetah, thetastep, thetaZl, thetaZh, thetaZstep, inputfilename, outputfilename):
    r_bins = np.arange(rl, rh + rstep, rstep)
    r_values = np.zeros(len(r_bins)-1, dtype=float)
    r_counter = np.zeros(len(r_bins)-1)
    theta_bins = np.arange(thetal, thetah + thetastep, thetastep)
    theta_counter = np.zeros((len(theta_bins)-1, len(r_bins)-1))
    thetaZ_bins = np.arange(thetaZl, thetaZh + thetaZstep, thetaZstep)
    thetaZ_counter = np.zeros(len(thetaZ_bins)-1)

    r_values, r_counter, theta_counter, thetaZ_counter = build_isotropic_histogram(inputfilename, r_bins, r_values, r_counter,
                                                                   theta_bins, theta_counter,
                                                                   thetaZ_bins, thetaZ_counter)
    r_values = r_values / r_counter

    # In case there is nan in the array, replace it with 0
    r_values = np.nan_to_num(r_values)
    theta_counter = np.nan_to_num(theta_counter)

    r_centers = (r_bins[:-1] + r_bins[1:])/2
    theta_centers = (theta_bins[:-1] + theta_bins[1:])/2
    thetaZ_centers = (thetaZ_bins[:-1] + thetaZ_bins[1:])/2

    write_isotropic_output(outputfilename, r_centers, r_values, r_counter)
    write_theta_output(outputfilename+'_theta', r_centers, theta_centers, theta_counter)
    write_thetaZ_output(outputfilename+'_thetaZ', thetaZ_centers, thetaZ_counter)
    return r_values, r_counter, r_centers, theta_counter, theta_centers,  thetaZ_counter, thetaZ_centers

def run_anisotropic(rl, rh, rstep, zl, zh, zstep, inputfilename, outputfilename):
    #rstep = 1; zstep = 1
    r_bins = np.arange(rl, rh + rstep, rstep)
    z_bins = np.arange(zl, zh + zstep, zstep)
    values = np.zeros((len(z_bins)-1, len(r_bins)-1), dtype=float)
    counter = np.zeros((len(z_bins)-1, len(r_bins)-1))
    r_centers = (r_bins[:-1] + r_bins[1:])/2
    z_centers = (z_bins[:-1] + z_bins[1:])/2

    values, counter = build_anisotropic_histogram(inputfilename, r_bins, z_bins, values, counter)
    values = values / counter
    # In case there is nan in the array, replace it with 0
    values = np.nan_to_num(values)
    write_anisotropic_output(outputfilename, r_centers, z_centers, values, counter)
    return values, counter, r_centers, z_centers

if system_type == 'isotropic':
    r_values, r_counter, r_centers, theta_counter, theta_centers, thetaZ_counter, thetaZ_centers = run_isotropic(rl, rh, rstep, thetal, thetah, thetastep, thetaZl, thetaZh, thetaZstep, inputfilename, outputfilename)
    plot_1d(r_centers, r_values, '3body_1d.png')
    plot_theta(r_centers, theta_centers, theta_counter, '3body_theta.png')
    plot_thetaZ(thetaZ_centers, thetaZ_counter, '3body_thetaZ.png')

    values, counter, r_centers, z_centers = run_anisotropic(rl, rh, rstep, zl, zh, zstep, inputfilename, outputfilename+'_aniso')
    plot_2d(r_centers, z_centers, values, '3body_2d.png')
    plot_projection(r_centers, z_centers, 5, values, '3body_2d_projection.png')

if system_type == 'anisotropic':
    values, counter, r_centers, z_centers = run_anisotropic(rl, rh, rstep, zl, zh, zstep, inputfilename, outputfilename)
    plot_2d(r_centers, z_centers, values, '3body_2d.png')
    plot_projection(r_centers, z_centers, 5, values, '3body_2d_projection.png')

    r_values, r_counter, r_centers, theta_counter, theta_centers, thetaZ_counter, thetaZ_centers = run_isotropic(rl, rh, rstep, thetal, thetah, thetastep, thetaZl, thetaZh, thetaZstep, inputfilename, outputfilename+'_iso')
    plot_1d(r_centers, r_values, '3body_1d.png')
    plot_theta(r_centers, theta_centers, theta_counter, '3body_theta.png')
    plot_thetaZ(thetaZ_centers, thetaZ_counter, '3body_thetaZ.png')

time2 = time.time()
calculation_time = time2-time1
print "Calculation time = " + str(calculation_time) + "seconds"
