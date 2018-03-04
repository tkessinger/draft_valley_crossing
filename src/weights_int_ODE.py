'''
weights_ODE.py
Author: Taylor Kessinger
Date: January 24, 2017
Description: Numerically solve weights ODE and plot for communally recombining phi.
'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as int
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--sigma', type=float, default=0.05)
args=parser.parse_args()

sigma= args.sigma #fitness standard deviation
mu = 1e-4 #mutation rate for second mutant
p = 0.01 #fixation probability
stepsize = .01
t_final = 100 #end time for solving ODE
r = 0.00001


numsteps = 1000 #number of fitness intervals
terr = 1e-3

s = 0 #advantage of double mutant
delta = 0 #disadvantage of single mutant

x_vals = np.linspace(-5,5,numsteps+1)

#fit_dist = 1.0/(np.sqrt(2*np.pi)*sigma)*np.exp(-(x_vals)**2/(2*sigma**2))*(x_vals[1]-x_vals[0]) #Gaussian fitness distribution
fit_dist = 1.0/np.sqrt(2*np.pi)*np.exp(-(x_vals)**2/(2))*(x_vals[1]-x_vals[0])

#fixation probability from Good et al. (2012)
def p_x(x, x_c, sigma):
    if x+s < 0:
        return 0
    elif x+s < x_c:
        return x_c*np.exp(((x+s)**2-x_c**2)/(2*sigma**2))
    else:
        return x+s

#implementing the equation for phi from Neher and Shraiman (2011)
# def dphi(t,phi,f_args):
#     mu,x,sigma,p = f_args[0],f_args[1],f_args[2],f_args[3]
#     z = mu*p
#     v = sigma**2
#      
#     theta = x-v*t+z-delta
#     return -(z + theta*phi - phi**2)

#alternate form of dphi
def dphi(t,phi,f_args):
    mu,x,sigma,p = f_args[0],f_args[1],f_args[2],f_args[3]
    z = mu*p
    #theta = x/sigma+sigma*t+z/sigma-delta/sigma
    theta = x/sigma+t*sigma+z/sigma-delta/sigma
    #theta = x+z
    return z/sigma + theta*phi - phi**2

def dPhi(t, Phi, f_args):
    mu,x_vals,sigma,p,r = f_args[0],f_args[1],f_args[2],f_args[3],f_args[4]
    z = mu*p
    #theta = x/sigma+sigma*t+z/sigma-delta/sigma
    #theta = x_vals+t*sigma+z/sigma-delta/sigma
    #theta = x_vals+t*sigma+z/sigma-delta/sigma
    #theta = x_vals+t*sigma+z/sigma-delta/sigma-r/sigma
    #theta = np.exp(x_vals+t*sigma+z/sigma-delta/sigma-r/sigma)
    theta = x_vals + z/sigma - delta/sigma - r/sigma
    #theta = x_vals/sigma+z/sigma-delta/sigma-r/sigma
    
    
    #print theta[0],theta[-1], t
    #theta = x+z
    return z/sigma + r/sigma*np.sum(np.multiply(fit_dist,Phi.T)) + theta*Phi - (1+theta)*Phi**2
    #a + theta*Phi - Phi**2

def get_Phi(r,sigma,t_final):
    soln = int.ode(dPhi).set_integrator('zvode', nsteps = 500000, with_jacobian=False)
    #p = p_x(x,3*sigma,sigma) UNCOMMENT THIS to use p_x above
    soln.set_initial_value(np.zeros(len(x_vals)),0).set_f_params([mu,x_vals,sigma,p,r]) #phi(0) = 0
    tf = t_final
    dt = stepsize
    tmp_soln = []
    while soln.successful() and soln.t <= tf:
        soln.integrate(soln.t+dt)
        tmp_soln.append(soln.y)
        #print x, soln.t, soln.y
        print soln.t
    tmp_soln = np.array(tmp_soln)
    return tmp_soln

'''plt.figure()
plt.plot(x_vals,tmp_phi_vals)
plt.show()
phi_vals.append(tmp_phi_vals)
plt.plot(x_vals, mu*p*np.exp((x_vals - mu*p)**2))
plt.plot(x_vals, sigma*(x_vals + mu*p))
'''

r_vals = [0,0.01]
sigma_vals = [0.01]
solns = []


for r in r_vals:
    tmp_solns = []
    for sigma in sigma_vals:
        tmp_solns.append(get_Phi(r,sigma,t_final))
    tmp_solns = np.array(tmp_solns)
    solns.append(tmp_solns)

solns = np.array(solns)

t_final += 1 #sometimes one has to manually increment t for reasons that are not clear

def plot_phi_t():
    #plot phi(t) versus t for different values of x
    plt_x_vals = [np.int(len(x_vals)/2),np.int(3*len(x_vals)/4),np.int(len(x_vals)-1)]
    plt.figure()
    for ri, r in enumerate(r_vals):
        for si, sigma in enumerate(sigma_vals):
            for xi, x in enumerate(plt_x_vals):
                time_array = 1.0*np.arange(len(solns[ri,si,:,x]))/len(solns[ri,si,:,x])*t_final
                plt.plot(time_array,solns[ri,si,:,x]*sigma,label=r'$x/\sigma = $'+str(x_vals[x])+ r'$, r = $'+str(r))
            plt.xlabel(r'$t$')
            plt.ylabel(r'$\phi (t)$')
            plt.legend(loc=1)
            #ax=plt.gca()
            #ax.set_yscale('log')
            plt.xlim([0,t_final])
            plt.tight_layout()
            plt.show()

def plot_Phi_t():
    plt.figure()
    for ri, r in enumerate(r_vals):
        for si, sigma in enumerate(sigma_vals):
            time_array = 1.0*np.arange(len(solns[ri,si,:,0]))/len(solns[ri,si,:,0])*t_final
            plt.plot(time_array, [np.sum(np.multiply(fit_dist,solns[ri,si,t,:].T))*sigma for t in range(len(time_array))],label=r'$\Phi(t), r = $'+str(r))
            plt.xlabel(r'$t$')
            plt.ylabel(r'$\Phi (t)$')
            plt.legend(loc=1)
            #ax=plt.gca()
            #ax.set_yscale('log')
            plt.xlim([0,t_final])
            plt.tight_layout()
            plt.show()

def plot_phi_x():
    plt.figure()
    for ri, r in enumerate(r_vals):
        plt.hlines(mu*p,min(x_vals),max(x_vals))
        for si, sigma in enumerate(sigma_vals):
            plt.plot(x_vals,solns[ri,si,-1,:]*sigma, label=r'$\sigma = $'+str(sigma) + r'$, r = $' + str(r))
    plt.legend(loc=2)
    plt.show()
    plt.xlabel(r'$x/\sigma$')
    plt.ylabel(r'$\phi(x,t)$')
    #ax=plt.gca()
    #ax.set_yscale('log')
    plt.tight_layout()

def plot_Phi_x():
    plt.figure()
    plt.plot(x_vals,fit_dist, ls = '--')
    for ri, r in enumerate(r_vals):
        for si, sigma in enumerate(sigma_vals):
            plt.plot(x_vals,np.multiply(fit_dist,solns[ri,si,-1,:].T)*sigma, label=r'$\sigma = $'+str(sigma) + r'$, r = $' + str(r))
    plt.legend(loc=2)
    plt.xlabel(r'$x/\sigma$')
    plt.ylabel(r'$\Phi(t)$')
    plt.tight_layout()
    plt.show()

plot_phi_t()
plot_Phi_t()
plot_phi_x()
plot_Phi_x()

#explanation of crossing probabilities:
#in principle phi = 1-\hat{p} = 1-<e^{-zw}>
#which is simply the overall probability that a bubble gives rise to a successful double mutant
#\frac{1}{N\mu p} therefore gives the expected crossing time
for ri, r in enumerate(r_vals):
    for si, sigma in enumerate(sigma_vals):
        crossing_prob = np.sum(np.multiply(fit_dist,solns[ri,si,-1,:].T))*sigma
        print "crossing prob taking final value of phi"
        print crossing_prob, r, sigma

