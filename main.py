from array import array
import ROOT
import math
import numpy as np
from funcs import getEtaPhiBins, sumTowers

ROOT.gStyle.SetOptStat(0)

def getNtupleVars(ntupleDir, ntupleName, varsDir):
    fIn = ROOT.TFile.Open(ntupleDir + ntupleName + '.root', "READ")
    treeIn = fIn.Get("hgcalTriggerNtuplizer/HGCalTriggerNtuple")

    varsName = ntupleName.replace("ntuple_", "vars_") #variables
    fOut = ROOT.TFile(varsDir + varsName + ".root", "RECREATE")
    treeOut = ROOT.TTree(varsName, varsName)

    GenEnergy_1 = array('f', [-100.])
    GenEta_1 = array('f', [-100.])
    GenPhi_1 = array('f', [-100.])
    EM1x1_1 = array('f', [-100.])
    EM3x3_1 = array('f', [-100.])
    EM5x5_1 = array('f', [-100.])
    EM7x7_1 = array('f', [-100.])
    EM9x9_1 = array('f', [-100.])
    EM11x11_1 = array('f', [-100.])
    Had1x1_1 = array('f', [-100.])
    Had3x3_1 = array('f', [-100.])
    Had5x5_1 = array('f', [-100.])
    Had7x7_1 = array('f', [-100.])
    Had9x9_1 = array('f', [-100.])
    Had11x11_1 = array('f', [-100.])
    
    GenEnergy_2 = array('f', [-100.])
    GenEta_2 = array('f', [-100.])
    GenPhi_2 = array('f', [-100.])
    EM1x1_2 = array('f', [-100.])
    EM3x3_2 = array('f', [-100.])
    EM5x5_2 = array('f', [-100.])
    EM7x7_2 = array('f', [-100.])
    EM9x9_2 = array('f', [-100.])
    EM11x11_2 = array('f', [-100.])
    Had1x1_2 = array('f', [-100.])
    Had3x3_2 = array('f', [-100.])
    Had5x5_2 = array('f', [-100.])
    Had7x7_2 = array('f', [-100.])
    Had9x9_2 = array('f', [-100.])
    Had11x11_2 = array('f', [-100.])

    treeOut.Branch('GenEnergy_1', GenEnergy_1, 'GenEnergy_1/F')
    treeOut.Branch('GenEta_1', GenEta_1, 'GenEta_1/F')
    treeOut.Branch('GenPhi_1', GenPhi_1, 'GenPhi_1/F')
    treeOut.Branch('EM1x1_1', EM1x1_1, 'EM1x1_1/F')
    treeOut.Branch('EM3x3_1', EM3x3_1, 'EM3x3_1/F')
    treeOut.Branch('EM5x5_1', EM5x5_1, 'EM5x5_1/F')
    treeOut.Branch('EM7x7_1', EM7x7_1, 'EM7x7_1/F')
    treeOut.Branch('EM9x9_1', EM9x9_1, 'EM9x9_1/F')
    treeOut.Branch('EM11x11_1', EM11x11_1, 'EM11x11_1/F')
    treeOut.Branch('Had1x1_1', Had1x1_1, 'Had1x1_1/F')
    treeOut.Branch('Had3x3_1', Had3x3_1, 'Had3x3_1/F')
    treeOut.Branch('Had5x5_1', Had5x5_1, 'Had5x5_1/F')
    treeOut.Branch('Had7x7_1', Had7x7_1, 'Had7x7_1/F')
    treeOut.Branch('Had9x9_1', Had9x9_1, 'Had9x9_1/F')
    treeOut.Branch('Had11x11_1', Had11x11_1, 'Had11x11_1/F')
    
    treeOut.Branch('GenEnergy_2', GenEnergy_2, 'GenEnergy_2/F')
    treeOut.Branch('GenEta_2', GenEta_2, 'GenEta_2/F')
    treeOut.Branch('GenPhi_2', GenPhi_2, 'GenPhi_2/F')
    treeOut.Branch('EM1x1_2', EM1x1_2, 'EM1x1_2/F')
    treeOut.Branch('EM3x3_2', EM3x3_2, 'EM3x3_2/F')
    treeOut.Branch('EM5x5_2', EM5x5_2, 'EM5x5_2/F')
    treeOut.Branch('EM7x7_2', EM7x7_2, 'EM7x7_2/F')
    treeOut.Branch('EM9x9_2', EM9x9_2, 'EM9x9_2/F')
    treeOut.Branch('EM11x11_2', EM11x11_2, 'EM11x11_2/F')
    treeOut.Branch('Had1x1_2', Had1x1_2, 'Had1x1_2/F')
    treeOut.Branch('Had3x3_2', Had3x3_2, 'Had3x3_2/F')
    treeOut.Branch('Had5x5_2', Had5x5_2, 'Had5x5_2/F')
    treeOut.Branch('Had7x7_2', Had7x7_2, 'Had7x7_2/F')
    treeOut.Branch('Had9x9_2', Had9x9_2, 'Had9x9_2/F')
    treeOut.Branch('Had11x11_2', Had11x11_2, 'Had11x11_2/F')

    etaBinStep = 0.0870
    minBinEta = -35
    maxBinEta = 35
    minEta = minBinEta * etaBinStep
    maxEta = maxBinEta * etaBinStep
    nBinsEta = maxBinEta - minBinEta
    
    phiBinStep = 2*math.pi/72
    minBinPhi = -36
    maxBinPhi = 36
    minPhi = minBinPhi * phiBinStep
    maxPhi = maxBinPhi * phiBinStep
    nBinsPhi = maxBinPhi - minBinPhi

    histEM = ROOT.TH2D("histEM","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)
    histHad = ROOT.TH2D("histHad","",nBinsEta,minEta,maxEta, nBinsPhi,minPhi,maxPhi)

    for entryNum in range(0, treeIn.GetEntries()):
        if entryNum%100==0:
            print("entryNum=", entryNum)
        treeIn.GetEntry(entryNum)
        tower_iPhi = getattr(treeIn,"tower_iPhi")
        tower_iEta = getattr(treeIn,"tower_iEta")
        tower_etEm = getattr(treeIn,"tower_etEm")
        tower_etHad = getattr(treeIn,"tower_etHad")
        tower_n = getattr(treeIn,"tower_n")
        tower_eta = getattr(treeIn,"tower_eta")        
        gen_eta = getattr(treeIn,"gen_eta")
        gen_phi = getattr(treeIn,"gen_phi")
        gen_energy = getattr(treeIn,"gen_energy")
        
        histEM.Reset()
        histHad.Reset()

        for towerID in range(tower_n):
            etaPhiBins = getEtaPhiBins(tower_eta[towerID], tower_iEta[towerID], tower_iPhi[towerID])
            histEM.SetBinContent(etaPhiBins[0], etaPhiBins[1], tower_etEm[towerID])
            histHad.SetBinContent(etaPhiBins[0], etaPhiBins[1], tower_etHad[towerID])

        index_1 = 0
        GenEnergy_1[0] = gen_energy[index_1]/np.cosh(gen_eta[index_1])
        GenEta_1[0] = gen_eta[index_1]
        GenPhi_1[0] = gen_phi[index_1]

        EM1x1_1[0] = sumTowers(histEM, gen_eta[index_1], gen_phi[index_1], numNeighbors=0)
        EM3x3_1[0] = sumTowers(histEM, gen_eta[index_1], gen_phi[index_1], numNeighbors=1)
        EM5x5_1[0] = sumTowers(histEM, gen_eta[index_1], gen_phi[index_1], numNeighbors=2)
        EM7x7_1[0] = sumTowers(histEM, gen_eta[index_1], gen_phi[index_1], numNeighbors=3)
        EM9x9_1[0] = sumTowers(histEM, gen_eta[index_1], gen_phi[index_1], numNeighbors=4)
        EM11x11_1[0] = sumTowers(histEM, gen_eta[index_1], gen_phi[index_1], numNeighbors=5)

        Had1x1_1[0] = sumTowers(histHad, gen_eta[index_1], gen_phi[index_1], numNeighbors=0)
        Had3x3_1[0] = sumTowers(histHad, gen_eta[index_1], gen_phi[index_1], numNeighbors=1)
        Had5x5_1[0] = sumTowers(histHad, gen_eta[index_1], gen_phi[index_1], numNeighbors=2)
        Had7x7_1[0] = sumTowers(histHad, gen_eta[index_1], gen_phi[index_1], numNeighbors=3)
        Had9x9_1[0] = sumTowers(histHad, gen_eta[index_1], gen_phi[index_1], numNeighbors=4)
        Had11x11_1[0] = sumTowers(histHad, gen_eta[index_1], gen_phi[index_1], numNeighbors=5)
        
        index_2 = 1
        GenEnergy_2[0] = gen_energy[index_2]/np.cosh(gen_eta[index_2])
        GenEta_2[0] = gen_eta[index_2]
        GenPhi_2[0] = gen_phi[index_2]

        EM1x1_2[0] = sumTowers(histEM, gen_eta[index_2], gen_phi[index_2], numNeighbors=0)
        EM3x3_2[0] = sumTowers(histEM, gen_eta[index_2], gen_phi[index_2], numNeighbors=1)
        EM5x5_2[0] = sumTowers(histEM, gen_eta[index_2], gen_phi[index_2], numNeighbors=2)
        EM7x7_2[0] = sumTowers(histEM, gen_eta[index_2], gen_phi[index_2], numNeighbors=3)
        EM9x9_2[0] = sumTowers(histEM, gen_eta[index_2], gen_phi[index_2], numNeighbors=4)
        EM11x11_2[0] = sumTowers(histEM, gen_eta[index_2], gen_phi[index_2], numNeighbors=5)

        Had1x1_2[0] = sumTowers(histHad, gen_eta[index_2], gen_phi[index_2], numNeighbors=0)
        Had3x3_2[0] = sumTowers(histHad, gen_eta[index_2], gen_phi[index_2], numNeighbors=1)
        Had5x5_2[0] = sumTowers(histHad, gen_eta[index_2], gen_phi[index_2], numNeighbors=2)
        Had7x7_2[0] = sumTowers(histHad, gen_eta[index_2], gen_phi[index_2], numNeighbors=3)
        Had9x9_2[0] = sumTowers(histHad, gen_eta[index_2], gen_phi[index_2], numNeighbors=4)
        Had11x11_2[0] = sumTowers(histHad, gen_eta[index_2], gen_phi[index_2], numNeighbors=5)
        
        treeOut.Fill()
    
    fOut.Write()
    fOut.Close()



def main():
    ntupleDir = "inputNtuples/"
    ntupleName = "ntuple_VBFHToInv_NoPU_WithEnergySplit"
    varsDir = "varsDir/"
    getNtupleVars(ntupleDir, ntupleName, varsDir)
    
if __name__=='__main__':
    main()
