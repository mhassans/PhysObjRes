import sys
import numpy as np

def getEtaPhiBins(tower_eta, tower_iEta, tower_iPhi):
    etaBin = -100
    phiBin = -100
    if (int(np.sign(tower_eta))==1):
        etaBin = 51 + tower_iEta
        phiBin = 1 + tower_iPhi
    elif (int(np.sign(tower_eta))==-1):
        etaBin = 20 - tower_iEta
        phiBin = 1 + tower_iPhi
    else:
        print('ERROR: found eta=0?')
        sys.exit(1)
    return (etaBin, phiBin)

def sumTowers(hist, gen_eta, gen_phi, numNeighbors):
    """
    numNeighbors = 0 means 1 (=1x1) tower which eta-phi is pointing to.
    numNeighbors = 1 means 9 (=3x3) towers with 8 around the center + 1 center.
    numNeighbors = 2 means 25 (=5x5) towers with 24 around the center tower + 1 center.
    etc.
    """
    
    if int(numNeighbors)!=numNeighbors:
        print("ERROR: number of towers must be an integer!")
        sys.exit(1)
    
    genEtaBin = hist.GetXaxis().FindBin(gen_eta)
    genPhiBin = hist.GetYaxis().FindBin(gen_phi)
    
    nBinsPhi = hist.GetNbinsY()
    nBinsEta = hist.GetNbinsX() #Number of eta bins inclusive of both endcaps and the barrel
    if (nBinsPhi!=72):
        print("number of phi bins updated/changed?")
        sys.exit(1)
    if (nBinsEta!=70):
        print("number of eta bins updated/changed?")
        sys.exit(1)

    energy = []
    eta = []
    phi = []
    for ix in range(-numNeighbors, 1+numNeighbors):
        for iy in range(-numNeighbors, 1+numNeighbors):
            tower_eta = ix+genEtaBin
            tower_phi = iy+genPhiBin
            
            if tower_phi > nBinsPhi:
                tower_phi = tower_phi - nBinsPhi
            elif tower_phi < 1:
                tower_phi = tower_phi + nBinsPhi
            
            if (tower_eta >= 1) and (tower_eta <= nBinsEta):
                energy.append(hist.GetBinContent(tower_eta, tower_phi))
                eta.append(tower_eta)
                phi.append(tower_phi)

    energy_sum = sum(energy)
    eta_avg = sum(eta[i]*energy[i]/energy_sum for i in range(len(energy))) if energy_sum!=0 else -999
    phi_avg = sum(phi[i]*energy[i]/energy_sum for i in range(len(energy))) if energy_sum!=0 else -999
    
    return (energy_sum, eta_avg, phi_avg)

