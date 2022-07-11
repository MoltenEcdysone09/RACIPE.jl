#Defining custom Shifted Hill Equations
#Shifted Hill Equation for -ve regulation
Hsn(thr,N,lbd,cop) = ((1/lbd) + (1-(1/lbd))*(1/(1 + (N/thr)^cop)))
#Shifted Hill Equation for +ve regulation
Hsp(thr,N,lbd,cop) = ((lbd) + (1-lbd)*(1/(1 + (N/thr)^cop)))/lbd

#Defining the reaction network with ModellingToolkit.jl macro
function rn!(du, u, p, t)
	MIR200, TGFB, MIR34A, VIME, OVOL2, ECAL, ZEB1, SNAIL1 = u
	gMIR200, gTGFB, gMIR34A, gVIME, gOVOL2, gECAL, gZEB1, gSNAIL1 ,kMIR200, kTGFB, kMIR34A, kVIME, kOVOL2, kECAL, kZEB1, kSNAIL1 ,tZEB1MIR200, fZEB1MIR200, nZEB1MIR200, tSNAIL1MIR200, fSNAIL1MIR200, nSNAIL1MIR200, tOVOL2TGFB, fOVOL2TGFB, nOVOL2TGFB, tMIR200TGFB, fMIR200TGFB, nMIR200TGFB, tZEB1MIR34A, fZEB1MIR34A, nZEB1MIR34A, tSNAIL1MIR34A, fSNAIL1MIR34A, nSNAIL1MIR34A, tZEB1VIME, fZEB1VIME, nZEB1VIME, tSNAIL1VIME, fSNAIL1VIME, nSNAIL1VIME, tOVOL2VIME, fOVOL2VIME, nOVOL2VIME, tZEB1OVOL2, fZEB1OVOL2, nZEB1OVOL2, tZEB1ECAL, fZEB1ECAL, nZEB1ECAL, tSNAIL1ECAL, fSNAIL1ECAL, nSNAIL1ECAL, tMIR200ZEB1, fMIR200ZEB1, nMIR200ZEB1, tOVOL2ZEB1, fOVOL2ZEB1, nOVOL2ZEB1, tSNAIL1ZEB1, fSNAIL1ZEB1, nSNAIL1ZEB1, tTGFBSNAIL1, fTGFBSNAIL1, nTGFBSNAIL1, tMIR34ASNAIL1, fMIR34ASNAIL1, nMIR34ASNAIL1, tSNAIL1SNAIL1, fSNAIL1SNAIL1, nSNAIL1SNAIL1, tOVOL2SNAIL1, fOVOL2SNAIL1, nOVOL2SNAIL1 = p
	du[1] = gMIR200*Hsn(tZEB1MIR200,ZEB1,fZEB1MIR200,nZEB1MIR200)*Hsn(tSNAIL1MIR200,SNAIL1,fSNAIL1MIR200,nSNAIL1MIR200) - kMIR200*MIR200
	du[2] = gTGFB*Hsn(tOVOL2TGFB,OVOL2,fOVOL2TGFB,nOVOL2TGFB)*Hsn(tMIR200TGFB,MIR200,fMIR200TGFB,nMIR200TGFB) - kTGFB*TGFB
	du[3] = gMIR34A*Hsn(tZEB1MIR34A,ZEB1,fZEB1MIR34A,nZEB1MIR34A)*Hsn(tSNAIL1MIR34A,SNAIL1,fSNAIL1MIR34A,nSNAIL1MIR34A) - kMIR34A*MIR34A
	du[4] = gVIME*Hsp(tZEB1VIME,ZEB1,fZEB1VIME,nZEB1VIME)*Hsp(tSNAIL1VIME,SNAIL1,fSNAIL1VIME,nSNAIL1VIME)*Hsn(tOVOL2VIME,OVOL2,fOVOL2VIME,nOVOL2VIME) - kVIME*VIME
	du[5] = gOVOL2*Hsn(tZEB1OVOL2,ZEB1,fZEB1OVOL2,nZEB1OVOL2) - kOVOL2*OVOL2
	du[6] = gECAL*Hsn(tZEB1ECAL,ZEB1,fZEB1ECAL,nZEB1ECAL)*Hsn(tSNAIL1ECAL,SNAIL1,fSNAIL1ECAL,nSNAIL1ECAL) - kECAL*ECAL
	du[7] = gZEB1*Hsn(tMIR200ZEB1,MIR200,fMIR200ZEB1,nMIR200ZEB1)*Hsn(tOVOL2ZEB1,OVOL2,fOVOL2ZEB1,nOVOL2ZEB1)*Hsp(tSNAIL1ZEB1,SNAIL1,fSNAIL1ZEB1,nSNAIL1ZEB1) - kZEB1*ZEB1
	du[8] = gSNAIL1*Hsp(tTGFBSNAIL1,TGFB,fTGFBSNAIL1,nTGFBSNAIL1)*Hsn(tMIR34ASNAIL1,MIR34A,fMIR34ASNAIL1,nMIR34ASNAIL1)*Hsn(tSNAIL1SNAIL1,SNAIL1,fSNAIL1SNAIL1,nSNAIL1SNAIL1)*Hsn(tOVOL2SNAIL1,OVOL2,fOVOL2SNAIL1,nOVOL2SNAIL1) - kSNAIL1*SNAIL1
end