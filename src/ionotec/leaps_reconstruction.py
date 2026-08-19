import numpy as np

def fit_lin(t,sig):
    N=len(t)
    if N<2: return None, None
    # Coefficients of paraboloid of error function (N*mse)
    a,b,c,d,e=0,N,0,0,0 # b=N for B^2 coef
    # Iterate over subseries to calculate coefficients
    for i in range(N):
        a+=t[i]**2 # A^2 coef
        c+=2*t[i] # A*B coef
        d-=2*t[i]*sig[i] # A coef
        e-=2*sig[i] # B coef

    if (a==0) and (c==0):
        return float('NaN'),float('NaN'),float('NaN'),float('NaN')
        
    # Forward A and B parameters of linear fit (solve the minimum of mse)
    A=-(2*b*d-c*e)/(4*a*b-c**2)
    B=-(c*d-2*a*e)/(c**2-4*a*b)

    
    #sigma=0
    #sigmas = []
    ## Calculate value of err function
    #for i in range(N):
    #    s = abs(sig[i]-A*t[i]-B)
    #    sigma+=s
    #    sigmas.append(s)
    #max_deviation = max(sigmas)
    #mean_deviation = np.mean(sigmas)

    return A,B#,max_deviation,mean_deviation
