'''
weights_ODE.py
Author: Taylor Kessinger
Date: January 4, 2017
Description: Numerically solve weights ODE and plot.
'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as int
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--sigma', type=float, default=0.005)
args=parser.parse_args()

sigma= args.sigma #fitness standard deviation
mu = 1e-5 #mutation rate for second mutant
p = 0.001 #fixation probability
stepsize = 1
t_final = 100000 #end time for solving ODE

numsteps = 50 #number of fitness intervals
terr = 1e-3

s = 0 #advantage of double mutant
delta = 0 #disadvantage of single mutant

#fixation probability from Good et al. (2012)
def p_x(x, x_c, sigma):
    if x+s < 0:
        return 0
    elif x+s < x_c:
        return x_c*np.exp(((x+s)**2-x_c**2)/(2*sigma**2))
    else:
        return x+s

#implementing the equation for phi from Neher and Shraiman (2011)
def dphi(t,phi,f_args):
    mu,x,sigma,p = f_args[0],f_args[1],f_args[2],f_args[3]
    z = mu*p
    v = sigma**2
    
    theta = x-v*t+z-delta
    return z + theta*phi - phi**2

#alternate form of dphi
#def dphi(t,phi,f_args):
#    mu,x,sigma,p = f_args[0],f_args[1],f_args[2],f_args[3]
#    z = mu*p
#    theta = x/sigma+sigma*t+z/sigma-delta/sigma
#    #theta = x+z
#    return z/sigma + theta*phi - phi**2



x_vals = np.linspace(-10,10,numsteps+1)

fit_dist = 1.0/(np.sqrt(2*np.pi)*sigma)*np.exp(-(x_vals)**2/(2*sigma**2))*(x_vals[1]-x_vals[0]) #Gaussian fitness distribution

tmp_phi_vals = []
for x in x_vals:
    x = sigma*x #un-normalize x, so sigma matters
    soln = int.ode(dphi).set_integrator('lsoda', method='bdf', with_jacobian=False)
    #p = p_x(x,3*sigma,sigma) UNCOMMENT THIS to use p_x above
    soln.set_initial_value(0,0).set_f_params([mu,x,sigma,p]) #phi(0) = 0
    tf = t_final
    dt = stepsize
    tmp_soln = []
    while soln.successful() and soln.t <= tf:
        soln.integrate(soln.t+dt)
        tmp_soln.append(soln.y)
        print x, soln.t, soln.y
    tmp_phi_vals.append(tmp_soln)

tmp_phi_vals = np.array(tmp_phi_vals)
'''plt.figure()
plt.plot(x_vals,tmp_phi_vals)
plt.show()
phi_vals.append(tmp_phi_vals)
plt.plot(x_vals, mu*p*np.exp((x_vals - mu*p)**2))
plt.plot(x_vals, sigma*(x_vals + mu*p))
'''

t_final += 1 #sometimes one has to manually increment t for reasons that are not clear

#plot phi(t) versus t for different values of x
plt_x_vals = [0,10,20,30,40,50]
plt.figure()
for xi, x in enumerate(plt_x_vals):
    plt.plot(np.arange(0,t_final,stepsize),tmp_phi_vals[x,:],label=r"$x/\sigma = $"+str(x_vals[x]))
plt.xlabel(r'$t$')
plt.ylabel(r'$\phi (t)$')
plt.legend(loc=1)
#ax=plt.gca()
#ax.set_yscale('log')
plt.tight_layout()
plt.show()

#plot max_x(phi(t)) for different values of x
plt.figure()
plt.plot(x_vals,np.max(tmp_phi_vals,axis=1), label=r'$\sigma = $'+str(sigma))
plt.legend(loc=2)
plt.show()
plt.xlabel(r'$x/\sigma$')
plt.ylabel(r'$max_x (\phi(x,t))$')
#ax=plt.gca()
#ax.set_yscale('log')
plt.tight_layout()

#plot phi(t_final) for different values of x
plt.figure()
plt.plot(x_vals,tmp_phi_vals[:,-1], label=r'$\sigma = $'+str(sigma))
plt.legend(loc=2)
plt.show()
plt.xlabel(r'$x/\sigma$')
plt.ylabel(r'$max_t (\phi(x,t))$')
#ax=plt.gca()
#ax.set_yscale('log')
plt.tight_layout()

#plot phi(x) for different values of t
plot_times = [100,200,500,1000,2000,5000,10000,20000,50000,100000]
plt.figure()
for ti, time in enumerate(plot_times):
    plt.plot(x_vals,tmp_phi_vals[:,time], label=r'$\tau = $'+str(time))
plt.hlines(mu*p,-10*sigma,10*sigma,linestyles='dashed')
plt.legend(loc=2)
plt.show()
plt.xlabel(r'$x/\sigma$')
plt.ylabel(r'$\phi(t,x)$')
plt.tight_layout()
#ax=plt.gca()
#ax.set_yscale('log')

#experimental: try integrating over phi(t)dt and see what happens
int_phi = np.sum(tmp_phi_vals, axis=1)*stepsize
int_phi = int_phi.flatten()
plt_x_vals = [0,10,20,30,40,50]
plt.figure()
plt.plot(x_vals,int_phi)
plt.xlabel(r'$x$')
plt.ylabel(r'$\int_0^\infty \phi(t) dt$')
plt.legend(loc=1)
#ax=plt.gca()
#ax.set_yscale('log')
plt.tight_layout()
plt.show()

#explanation of crossing probabilities:
#in principle phi = 1-\hat{p} = 1-<e^{-zw}>
#which is simply the overall probability that a bubble gives rise to a successful double mutant
#\frac{1}{N\mu p} therefore gives the expected crossing time
crossing_prob = np.sum(fit_dist*np.max(tmp_phi_vals,axis=1))
print "crossing prob taking max(phi)"
print crossing_prob, sigma

crossing_prob = np.sum(fit_dist*tmp_phi_vals[:,100])
print "crossing prob taking final value of phi"
print crossing_prob, sigma

crossing_prob = np.sum(fit_dist*int_phi)
print "crossing prob taking integral over phi"
print crossing_prob, sigma

