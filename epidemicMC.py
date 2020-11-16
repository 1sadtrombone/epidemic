import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import get_context, Pool
from functools import partial

def MC(N, p=0.5):
    # returns <X>, <Y>, <Z>
    # after N realizations
    
    Lx = 100
    Ly = 100
    tmax = 100
    
    Y = np.zeros((tmax, N, Lx, Ly), dtype=bool) # infected
    Z = np.zeros((tmax, N, Lx, Ly), dtype=bool) # immune
    # everything else is susceptible
    
    # initialize
    Y[0,:,0] = 1
    
    for n in range(N):
        for t in range(tmax-1):
            ps = np.random.uniform(size=(4,Lx,Ly)) < p
            # do we get infected? One for each direction
            
            # from left?
            Y[t+1,n][np.where(ps[0] * np.roll(Y[t,n], 1, axis=1))] = 1
            # from right?
            Y[t+1,n][np.where(ps[1] * np.roll(Y[t,n], -1, axis=1))] = 1
            # from top? (no periodic bounds)
            Y[t+1,n][np.where(ps[2] * np.vstack((
                                    np.zeros(Lx),
                                    np.roll(Y[t,n], 1, axis=0)[1:]
                                    )))] = 1
            # from bottom? (no periodic bounds)
            Y[t+1,n][np.where(ps[3] * np.vstack((
                                    np.roll(Y[t,n], -1, axis=0)[:-1],
                                    np.zeros(Lx)
                                    )))] = 1
            # the above might reinfect sites 
            # which were immune or already infected
            # but that's ok, because
            
            Z[t+1,n] = Y[t,n]
            Z[t+1,n] += Z[t,n]
            Y[t+1,n][np.where(Z[t+1,n])] = 0
            
            if np.sum(Y[t+1,n]) == 0:
                
                # we know Y[t] == 0 forever now,
                # so I just fill out
                
                Z[t+1:tmax,n] = Z[t+1,n]
                continue
                
    return np.mean(Z, axis=1)

if __name__=="__main__":
    n_cores = 4
    tasks_per_core = 50
    
    p = 0.5
    
    with get_context("spawn").Pool(n_cores) as pool:

        # get averaged Z over N realizations
        Z = np.mean(np.array(pool.map(partial(MC, p=p), [tasks_per_core]*n_cores)), axis=0)

    print(Z.shape)
        
    plt.imshow(Z[90], aspect='auto')
    plt.colorbar()
    plt.show()

