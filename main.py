from array import array
import sys
import ROOT
import math
import numpy as np
from funcs import getEtaPhiBins, sumTowers
#from optparse import OptionParser
import argparse


ROOT.gStyle.SetOptStat(0)
    
parser = argparse.ArgumentParser()

parser.add_argument('--jet', action='store_true')
parser.add_argument('--nojet', action='store_true')
parser.add_argument('--input', type=str, required=True)

args = parser.parse_args()

if(bool(args.nojet)==bool(args.jet)):
    print("ERROR: either '--jet' or '--nojet' must be used as input argument")
    sys.exit(1)

ntupleFile = args.input
outDir = "varsDir/"
isjet = args.jet

fIn = ROOT.TFile.Open(ntupleFile, "READ")
treeIn = fIn.Get("hgcalTriggerNtuplizer/HGCalTriggerNtuple")

varsFile = ntupleFile.split('/')[-1].replace("ntuple_", "vars_").replace(".root","") #variables
fOut = ROOT.TFile(outDir + varsFile + ".root", "RECREATE")
treeOut = ROOT.TTree(varsFile, varsFile)

GenEnergy_1 = array('f', [-100.])
GenEta_1 = array('f', [-100.])
GenPhi_1 = array('f', [-100.])
GenEnergy_2 = array('f', [-100.])
GenEta_2 = array('f', [-100.])
GenPhi_2 = array('f', [-100.])
treeOut.Branch('GenEnergy_1', GenEnergy_1, 'GenEnergy_1/F')
treeOut.Branch('GenEta_1', GenEta_1, 'GenEta_1/F')
treeOut.Branch('GenPhi_1', GenPhi_1, 'GenPhi_1/F')
treeOut.Branch('GenEnergy_2', GenEnergy_2, 'GenEnergy_2/F')
treeOut.Branch('GenEta_2', GenEta_2, 'GenEta_2/F')
treeOut.Branch('GenPhi_2', GenPhi_2, 'GenPhi_2/F')

mtxSizes = ['1x1', '3x3', '5x5', '7x7', '9x9', '11x11']
detectorVariables = ['EtEM', 'EtHad', 'EtTot', 'EtaEM', 'EtaHad', 'EtaTot', 'PhiEM', 'PhiHad', 'PhiTot']
energyRank = ['_1', '_2']

for rank in energyRank:
    for det in detectorVariables:
        for size in mtxSizes:
            var = det+size+rank
            exec(var+" = array('f', [-100.])")
            exec("treeOut.Branch('"+var+"', "+var+", '"+var+"/F')")

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
    if(isjet):
        jetflag = 'jet'
    else:
        jetflag = ''
    gen_eta = getattr(treeIn,"gen" + jetflag + "_eta")
    gen_phi = getattr(treeIn,"gen" + jetflag + "_phi")
    gen_energy = getattr(treeIn,"gen" + jetflag + "_energy")
    
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
    
    index_2 = 1
    only1particle = True if len(gen_eta)==1 else False
    GenEnergy_2[0] = gen_energy[index_2]/np.cosh(gen_eta[index_2]) if not only1particle else -999
    GenEta_2[0] = gen_eta[index_2] if not only1particle else -999
    GenPhi_2[0] = gen_phi[index_2] if not only1particle else -999

    for det in ['EM', 'Had']:
        for count, size in enumerate(mtxSizes):
            for rank in energyRank:
                var = det+size+rank
                if((rank=='_1') or (not only1particle)):
                    exec("Et"+var+"[0], Eta"+var+"[0], Phi"+var+"[0] = sumTowers(hist"+det+\
                        ", gen_eta[index"+rank+"], gen_phi[index"+rank+"], numNeighbors="+str(count)+")")
                else:
                    exec("Et"+var+"[0], Eta"+var+"[0], Phi"+var+"[0] = (-999, -999, -999)")

    for size in mtxSizes:
        for rank in energyRank:
            var = size+rank
            if((rank=='_1') or (not only1particle)):
                exec("EtTot"+var+"[0] = EtEM"+var+"[0] + EtHad"+var+"[0]")
                exec("EtaTot"+var+"[0] = (EtaEM"+var+"[0]*EtEM"+var+"[0] + EtaHad"
                    +var+"[0]*EtHad"+var+"[0])/(EtTot"+var+"[0]) if (EtTot"+var+"[0])!=0 else -999")
                exec("PhiTot"+var+"[0] = (PhiEM"+var+"[0]*EtEM"+var+"[0] + PhiHad"
                    +var+"[0]*EtHad"+var+"[0])/(EtTot"+var+"[0]) if (EtTot"+var+"[0])!=0 else -999")
            else:
                exec("EtTot"+var+"[0] = -999")
                exec("EtaTot"+var+"[0] = -999")
                exec("PhiTot"+var+"[0] = -999")
                
    treeOut.Fill()

fOut.Write()
fOut.Close()




