import numpy as np
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['figure.dpi'] = 300

# Use parameters used in Example 1.20 of Ramo et al., 3rd edition.
book_config = False

V = [80, 60, 100, 20] # Vl, Vb, Vt, Vr
#V = [80, 0, 0, 0] # Vl, Vb, Vt, Vr
xo = 1
yo = 1
Nx = 64
Ny = 64

# Show grid points on plots
show_grid = False

if book_config:
    V = [80, 60, 100, 20] # Vl, Vb, Vt, Vr
    xo = 1 
    yo = 1
    Nx = 4
    Ny = 4
    show_grid = True

# Stop numerical soln when change at all grid points from one step to next
# is small than local_tol.
local_tol = 0.00001 

def log(msg):
    with open('HW2_2.txt','a') as f:
        f.write(msg + "\n")
    print(msg)    

    
def nancorners(M):
    
    # Set un-used parts of grid to nan so values not plotted
    M[0,0] = np.nan
    M[-1,0] = np.nan
    M[-1,-1] = np.nan
    M[0,-1] = np.nan

    return M


def step(Phi):
    Phi = np.copy(Phi)
    for i in range(1, Phi.shape[0] - 1):
        for j in range(1, Phi.shape[1] - 1):
            Phi[i, j] = (1/4)*(Phi[i-1, j] + Phi[i+1, j] + Phi[i, j-1] + Phi[i, j+1])
        
    return Phi


def exact(x, y):
    # Exact solution is in the form of a series sum.
    # Exact solution is approximated by first non-zero Nmax/2 terms in sum.
    Nmax = 66
    pi = np.pi
    sinh = np.sinh
    sin = np.sin
    Phi_e = np.nan*np.ones((len(x),len(y)))
    for i in np.arange(1, x.shape[0] - 1):
        for j in np.arange(1, x.shape[1] - 1):
            P = 0.
            for n in np.arange(1, Nmax+1, 2):
                R1 = (4/pi)*(1/n)*(1/sinh(n*pi*yo/xo));
                R2 = (4/pi)*(1/n)*(1/sinh(n*pi*xo/yo));
                Pl =  V[0]*sin(n*pi*y[i, j]/yo)*sinh(n*pi*(xo-x[i, j])/yo); 
                Pb =  V[1]*sin(n*pi*x[i, j]/xo)*sinh(n*pi*y[i, j]/xo); 
                Pt =  V[2]*sin(n*pi*x[i, j]/xo)*sinh(n*pi*(yo-y[i, j])/xo); 
                Pr =  V[3]*sin(n*pi*y[i, j]/yo)*sinh(n*pi*x[i, j]/yo); 
                P = P + R1*Pl + R1*Pr + R2*Pt + R2*Pb;
            Phi_e[i,j] = P
    return Phi_e


def plot(x, y, Phi, title='', show_grid=False):

    fig, ax = plt.subplots()

    if show_grid:
        # Show grid points
        ax.plot(X.flatten(),Y.flatten(),'k.')

    ax.set_title(title)
    ax.set_aspect(1)
    # Choose number of colors as 8 so colors line up with tick increments on colorbar.
    m = ax.pcolor(x, y, Phi, cmap=plt.get_cmap('viridis', 9))
    m.set_clim(15, 105) # Set colorbar limits
    cb = fig.colorbar(m)
    cb.ax.set_title('$\Phi$ [V]')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    return m, cb

x = np.arange(0, Nx) # [0, 1, ..., Nx-1]
y = np.arange(0, Ny) # [0, 1, ..., Ny-1]

x = xo*x/(Nx-1) # x = [0, xo/(Nx-1), ..., xo]
y = yo*y/(Nx-1) # y = [0, yo/(Ny-1), ..., yo]

X, Y = np.meshgrid(x, y, indexing='ij')

# Set all values to average of boundary values
Phi = np.mean(V)*np.ones((Ny, Nx), dtype='float64')

# Set boundary potentials
Phi[:,0]  = V[0] # Vl
Phi[-1,:] = V[1] # Vb
Phi[0,:]  = V[2] # Vt
Phi[:,-1] = V[3] # Vr

if book_config:
    for s in range(5):
        Phi_new = step(Phi)
        log('Step {}'.format(s+1))
        log('Φ1\t{0:.2f}'.format(Phi_new[1,1]))
        log('Φ2\t{0:.2f}'.format(Phi_new[1,2]))
        log('Φ3\t{0:.2f}'.format(Phi_new[2,1]))
        log('Φ4\t{0:.2f}'.format(Phi_new[2,2]))
        #d = np.max(np.abs(Phi.flatten()-Phi_new.flatten()))
        #log('abs(max) difference between steps {0:.2f}'.format(d))
        Phi = Phi_new
else:
    # d will be max absolute change. Set to value larger than stopping value.
    d = 2*local_tol 
    # s  is step number
    s = 1 
    while d >= local_tol:
        Phi_new = step(Phi)
        d = np.max(np.abs(Phi.flatten()-Phi_new.flatten()))
        if d < local_tol:
            break
        s = s + 1
        Phi = Phi_new
    log('{0:d} steps for max abs change in values to be < {1:f}'.format(s, local_tol))

# Flip solution for plotting.
Phi = np.flipud(Phi)

# Compute exact solution.
Phi_e = exact(X, Y)

# Transpose and flip solution for plotting.
Phi_e = np.flipud(Phi_e.T)

# Compute mean-squared error of numerical solution relative to exact solution.
mse = np.sqrt(np.mean((Phi[1:-1,1:-1].flatten() - Phi_e[1:-1,1:-1].flatten())**2))

# Change x and y grid so plotted colored patches are centered on grid points
dx = x[1]-x[0]
dy = y[1]-y[0]
x = x-dx/2
x = np.append(x, x[-1] + dx)
y = y-dy/2
y = np.append(y, y[-1] + dy)

# Set corners to be nan so they are not plotted. The corners are not a part
# of the solution (and values are ambiguous as they are a part of two boundaries.)
X = nancorners(X)
Y = nancorners(Y)
Phi = nancorners(Phi)

# Set boundary values in exact solution for plotting.
Phi_e[0, :]  = Phi[0, :]
Phi_e[-1, :] = Phi[-1, :]
Phi_e[:, 0]  = Phi[:, 0]
Phi_e[:, -1] = Phi[:, -1]

if book_config:
    log('Series solution'.format(s+1))
    log('Φ1\t{0:.2f}'.format(Phi_e[2,1]))
    log('Φ2\t{0:.2f}'.format(Phi_e[2,2]))
    log('Φ3\t{0:.2f}'.format(Phi_e[1,1]))
    log('Φ4\t{0:.2f}'.format(Phi_e[1,2]))
    
m, cb = plot(x, y, Phi, title='Numerical solution', show_grid=show_grid)
if book_config:
    plt.savefig('HW2_2.pdf', bbox_inches='tight')

plt.savefig('HW2_figures/HW2_2_{0:d}x{1:d}_Numerical.pdf'.format(Nx, Ny), bbox_inches='tight')

m, cb = plot(x, y, Phi, title='Series solution')
plt.savefig('HW2_figures/HW2_2_{0:d}x{1:d}_Series.pdf'.format(Nx, Ny), bbox_inches='tight')

m, cb = plot(x, y, Phi-Phi_e, title='(Numerical-Series); MSE = {0:.2f} [V]'.format(mse))
cb.ax.set_title('$\Delta\Phi$ [V]')
m.set_clim(-2.25, 2.25) # Set colorbar limits
plt.savefig('HW2_figures/HW2_2_{0:d}x{1:d}_Difference.pdf'.format(Nx, Ny), bbox_inches='tight')
