#! /usr/bin/env python

import sys
from hod_functions import *
from scipy.optimize import curve_fit as cf
from scipy.optimize import fmin as simplex

def chi2(ydata,model,errors=None):
    if (errors is None):
        chi2 = sum((ydata-model)**2)
    else:
        if (len(errors) != len(ydata)):
            sys.exit('STOP chi2: err should have the same lenght as ydata')

        chi2 = sum((ydata-model)**2/errors**2)
    return chi2

def chi2_poisson(ydata,model):
    chi2 = sum((ydata-model)**2/model)
    return chi2

def f_chi2(params, x, y, err, func):
    fchi2 = 0.0
    
    model = func(x,*params)
    fchi2 = chi2(y,model,errors=err)
    
    return fchi2


def go_simplex(func,xdata,ydata,lb,ub,errors=None):
    params = -999. ; chi2n = -999.

    # Check that the degrees of freedom is larger than 1
    dof = len(xdata)-len(lb)
    if (dof<=1):
        sys.exit('HOD_FIT.py: The degrees of freedom should be >1')

    # Get random initial values within the boundaries given
    p1 = random.rand(len(lb))
    p1 = lb + p1*(ub-lb) 

    params = simplex(f_chi2,p1,args=(xdata,ydata,errors,func),full_output=0)

    model = func(xdata,*params)
    chi2n  = chi2(ydata,model,errors)/dof

    return params,chi2n

def fitting(func,xdata,ydata,lb,ub,errors=None,Poisson=False,fixed_param=None):
    params = -999. ; chi2n = -999.

    # Check that the degrees of freedom is larger than 1
    dof = len(xdata)-len(lb)
    if (dof<=1):
        sys.exit('HOD_FIT.py: The degrees of freedom should be >1')

    # Find an adequate initial set of parameters
    ntest = 50
    chi0 = 999.
    for i in range(ntest):
        p1 = random.rand(len(lb))
        p1 = lb + p1*(ub-lb) 

        if (fixed_param is None):
            p2,cov = cf(func,xdata,ydata,p0=p1)
            model = func(xdata,*p2)
        else:
            print 'working on fixing parameters'
            #ia = fixed_param ; a = lb[ia]
            #p2,cov = cf(func+bound([a,a],),xdata,ydata,p0=p1)
            #model = func(xdata,*p2)

        if (Poisson):
            chi = abs(chi2_poisson(ydata,model)/dof -1.)
        else:
            chi = abs(chi2(ydata,model,errors)/dof -1.)

        if (chi<=chi0 or i==0): 
            chi0   = chi
            params = p2
            chi2n  = chi2(ydata,model,errors)/dof
        
    return params,chi2n


def brute_fit(func,xdata,ydata,lb,ub,steps,errors=None,Poisson=False,fixed_param=None):
    params = -999. ; chi2n = -999.  ; nchi0 = 999.

    # Check that the degrees of freedom is larger than 1
    dof = len(xdata)-len(lb)
    if (dof<=1):
        sys.exit('HOD_FIT.py: The degrees of freedom should be >1')
    
    # Find an adequate initial set of parameters
    for i in range(steps):
        p1 = random.rand(len(lb))
        p1 = lb + p1*(ub-lb) 
    
        if (fixed_param is not None):
            ia = fixed_param 
            p1[ia] = lb[ia]
    
        model = func(xdata,*p1)
        if (Poisson):
            nchi = abs(chi2_poisson(ydata,model)/dof -1.)
        else:
            nchi = abs(chi2(ydata,model,errors)/dof -1.)
    
        if (nchi<nchi0 or i==0): 
            nchi0 = nchi
            params=p1
            chi2n = chi2(ydata,model,errors)/dof
        
    return params,chi2n
