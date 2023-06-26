#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 11:42:17 2022

@author: silvio
"""
#%% Martijn suggestion

import matplotlib.pyplot as plt
import numpy as np

import os
import sys

# Add QDYN source directory to PATH
# Go up in the directory tree
upup = [os.pardir]*2
qdyn_dir = os.path.join(*upup)
# Get QDYN src directory
src_dir = os.path.abspath(
    os.path.join(
        os.path.join(os.path.abspath(""), qdyn_dir), "/home/silvio/QDyn/qdyn-feature-injectionGalis/src")
)
# Append src directory to Python path
sys.path.append(src_dir)
# Get QDYN plotting library directory
plot_dir = os.path.abspath(
    os.path.join(
        os.path.join(os.path.abspath(""), qdyn_dir), "utils", "post_processing")
)
# Append plotting library directory to Python path
sys.path.append(plot_dir)

# Import QDYN wrapper and plotting library
from pyqdyn import qdyn
import plot_functions as qdyn_plot

import matplotlib.pyplot as plt
from scipy.special import j0, j1, y0, y1
from scipy.optimize import fsolve
from scipy import optimize
from scipy.special import logsumexp

p = qdyn()

# Get the settings dict
set_dict = p.set_dict

re = 1000 #set_dict["SET_DICT_INJECTION"]["re"] 
rw = 0.10 #set_dict["SET_DICT_INJECTION"]["rw"]

reD = re / rw

J1 = j1
Y1 = y1 

def get_roots(x):
    return (J1(x*reD)*Y1(x)-J1(x)*Y1(x * reD))
    #return Fn


def derivative(f, x):
    h = 1e-8
    return (f(x + h) - f(x))/h

def solver(f, x0, epsilon, max_iter):
    xn = x0
    for n in range(0,max_iter):
        y = f(xn)
        if abs(y) < epsilon:
            return xn
        slope = derivative(f, xn)
        if (slope == 0) :
            return None
        xn = xn - y / slope
    return None



def loop(f, L_bound, R_bound, increment, nsolutions):
    solutions = []
   # while np.max(len(solutions)) != nsolutions:
    #    R_bound += 1
    while L_bound <= R_bound:
        solution = solver(f, L_bound, 1e-10, 1000)
        if solution is not None:
            solution = round(solution,15)
            if solution not in solutions:
                solutions.append(solution)
            #    if np.max(len(solutions)) == nsolutions:
             #       break
             #   else:
              #      while np.max(len(solutions)) != nsolutions:
               #         print(solutions)
                #        R_bound += 2000
 #   while np.max(len(solutions)) != nsolutions:
  #      R_bound += 1
   #     if np.max(len(solutions)) == nsolutions:
    #        break
        L_bound += increment
    print(sorted(solutions[:10]))
    print("we found " +str(len(solutions)) + " solutions!")
   

equation = ""
def f(x):
    try:
        y = eval(equation)
    except ZeroDivisionError:
        y = 1e-15
    return y
    
import math 
#loop(lambda x: math.sin(x), -10, 10, 0.5)   
#a = loop(get_roots, 0, 10, 0.005)  
equation="x**2-9"
loop(get_roots,0,0.1,0.0001,10)   

#x = np.arange(0.0005,0.007,0.00001)
#sol = []
#s = get_roots(x)
#sol.append(s)