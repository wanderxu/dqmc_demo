#!/usr/bin/dev python
#coding=utf-8
"""
Author:         Xiao-Yan Xu <wanderxu@gmail.com>
Description:
use chi square to do fitting.

Input:
      file "tmppy.dat" with first k columns of independent variables,
                            following columns with observables and errors
      take care of the number of independent variables k and the model to fitting defined in curvefunc
      for any special case, you need change corresponding places

"""
import math
import numpy as np
import scipy.optimize as opt

## define read file func
def file2list(filename):
		fr = open(filename)
		array = fr.readlines()
		num = len(array)
		returnMat=np.zeros((num,4))  # you can change the dimension 
		index = 0
		for line in array:
				line = line.strip()
				linelist = line.split()
				returnMat[index,:]=linelist[0:4] # you can change the dimension
				index +=1
		return returnMat

## define curvefunc for curve_fit
def curvefunc (xv, e0, ji1, ji2, jt1 ):
        return e0 + ji1*xv[0] + ji2*xv[1] + jt1*xv[2]

def chi_square ( xdata, ydata, ydata_sigma, p0 ):
        popt,pcov = opt.curve_fit(curvefunc, xdata, ydata, p0, sigma=ydata_sigma, absolute_sigma=False )
        perr = np.sqrt(np.diag(pcov))
        rchi_sq = np.sum( ( (ydata-curvefunc(xdata,popt[0],popt[1],popt[2],popt[3]) ) / ydata_sigma )**2 )/(len(ydata)-5.0)
        return popt, perr, rchi_sq


## prepare data
indat=file2list("train.dat")
inmat = np.transpose(indat)
xdata = inmat[1:4]
ydata = inmat[0]
ydata_sigma =len(ydata)*[1.0]
p0=4*[1.0]
popt,perr,rchi_sq = chi_square( xdata, ydata, ydata_sigma, p0 )
print " e0 = ", popt[0], "+/-", perr[0]
print " j1 = ", popt[1], "+/-", perr[1]
print " j2 = ", popt[2], "+/-", perr[2]
print " t1 = ", popt[3], "+/-", perr[3]
print " x^2 = ", rchi_sq
