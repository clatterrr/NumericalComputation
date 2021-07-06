import numpy as np

rho_fluid = 1
rho_air = 0
nu_fluid = 1
nu_air = 0


def getRhoAndGradRho(gradrho,rho,gradphi,phi):
    H = 0
    Hp = 0
    eps = 1e-10
    if phi > eps:
        H = 1
        Hp = 0
    elif phi < -eps:
        H = -1
        Hp = 0
    else:
        H = phi / eps
        Hp = 1 / eps
    rho = rho_fluid * (1 + H)/2 + rho_air * (1 - H)/2
    gradrho[0] = (rho_fluid - rho_air)/2 * Hp * gradphi[0]
    gradrho[1] = (rho_fluid - rho_air)/2 * Hp * gradphi[1]

def getRhoAndNu(rho,mu,phi):
    H = 0
    eps = 1e-10
    if phi > eps:
        H = 1
    elif phi < -eps:
        H = -1
    else:
        H = phi / eps
        
    rho = rho_fluid * (1 + H)/2 + rho_air * (1 - H) / 2
    mu = nu_fluid * (1 + H)/2 + nu_air*(1 - H)/2
    