#Defining custom Shifted Hill Equations
#Shifted Hill Equation for -ve regulation
Hsn(Thr,N,Fld,Hill) = ((Fld) + (1-Fld)*(1/(1 + (N/Thr)^Hill)))
#Shifted Hill Equation for +ve regulation
Hsp(Thr,N,Fld,Hill) = ((Fld) + (1-Fld)*(1/(1 + (N/Thr)^Hill)))/Fld

#Defining the reaction network with ModellingToolkit.jl macro
function rn!(du, u, p, t)
	MIR200, TGFB, MIR34A, VIME, OVOL2, ECAL, ZEB1, SNAIL1 = u
	Prod_MIR200, Prod_TGFB, Prod_MIR34A, Prod_VIME, Prod_OVOL2, Prod_ECAL, Prod_ZEB1, Prod_SNAIL1 ,Deg_MIR200, Deg_TGFB, Deg_MIR34A, Deg_VIME, Deg_OVOL2, Deg_ECAL, Deg_ZEB1, Deg_SNAIL1 ,InhThr_SNAIL1_MIR200, InhFld_SNAIL1_MIR200, Hill_SNAIL1_MIR200, InhThr_ZEB1_MIR200, InhFld_ZEB1_MIR200, Hill_ZEB1_MIR200, InhThr_MIR200_TGFB, InhFld_MIR200_TGFB, Hill_MIR200_TGFB, InhThr_OVOL2_TGFB, InhFld_OVOL2_TGFB, Hill_OVOL2_TGFB, InhThr_SNAIL1_MIR34A, InhFld_SNAIL1_MIR34A, Hill_SNAIL1_MIR34A, InhThr_ZEB1_MIR34A, InhFld_ZEB1_MIR34A, Hill_ZEB1_MIR34A, InhThr_OVOL2_VIME, InhFld_OVOL2_VIME, Hill_OVOL2_VIME, ActThr_SNAIL1_VIME, ActFld_SNAIL1_VIME, Hill_SNAIL1_VIME, ActThr_ZEB1_VIME, ActFld_ZEB1_VIME, Hill_ZEB1_VIME, InhThr_ZEB1_OVOL2, InhFld_ZEB1_OVOL2, Hill_ZEB1_OVOL2, InhThr_SNAIL1_ECAL, InhFld_SNAIL1_ECAL, Hill_SNAIL1_ECAL, InhThr_ZEB1_ECAL, InhFld_ZEB1_ECAL, Hill_ZEB1_ECAL, InhThr_MIR200_ZEB1, InhFld_MIR200_ZEB1, Hill_MIR200_ZEB1, InhThr_OVOL2_ZEB1, InhFld_OVOL2_ZEB1, Hill_OVOL2_ZEB1, ActThr_SNAIL1_ZEB1, ActFld_SNAIL1_ZEB1, Hill_SNAIL1_ZEB1, InhThr_MIR34A_SNAIL1, InhFld_MIR34A_SNAIL1, Hill_MIR34A_SNAIL1, InhThr_OVOL2_SNAIL1, InhFld_OVOL2_SNAIL1, Hill_OVOL2_SNAIL1, InhThr_SNAIL1_SNAIL1, InhFld_SNAIL1_SNAIL1, Hill_SNAIL1_SNAIL1, ActThr_TGFB_SNAIL1, ActFld_TGFB_SNAIL1, Hill_TGFB_SNAIL1 = p
	du[1] = Prod_MIR200*Hsn(InhThr_SNAIL1_MIR200,SNAIL1,InhFld_SNAIL1_MIR200,Hill_SNAIL1_MIR200)*Hsn(InhThr_ZEB1_MIR200,ZEB1,InhFld_ZEB1_MIR200,Hill_ZEB1_MIR200) - Deg_MIR200*MIR200
	du[2] = Prod_TGFB*Hsn(InhThr_MIR200_TGFB,MIR200,InhFld_MIR200_TGFB,Hill_MIR200_TGFB)*Hsn(InhThr_OVOL2_TGFB,OVOL2,InhFld_OVOL2_TGFB,Hill_OVOL2_TGFB) - Deg_TGFB*TGFB
	du[3] = Prod_MIR34A*Hsn(InhThr_SNAIL1_MIR34A,SNAIL1,InhFld_SNAIL1_MIR34A,Hill_SNAIL1_MIR34A)*Hsn(InhThr_ZEB1_MIR34A,ZEB1,InhFld_ZEB1_MIR34A,Hill_ZEB1_MIR34A) - Deg_MIR34A*MIR34A
	du[4] = Prod_VIME*Hsn(InhThr_OVOL2_VIME,OVOL2,InhFld_OVOL2_VIME,Hill_OVOL2_VIME)*Hsp(ActThr_SNAIL1_VIME,SNAIL1,ActFld_SNAIL1_VIME,Hill_SNAIL1_VIME)*Hsp(ActThr_ZEB1_VIME,ZEB1,ActFld_ZEB1_VIME,Hill_ZEB1_VIME) - Deg_VIME*VIME
	du[5] = Prod_OVOL2*Hsn(InhThr_ZEB1_OVOL2,ZEB1,InhFld_ZEB1_OVOL2,Hill_ZEB1_OVOL2) - Deg_OVOL2*OVOL2
	du[6] = Prod_ECAL*Hsn(InhThr_SNAIL1_ECAL,SNAIL1,InhFld_SNAIL1_ECAL,Hill_SNAIL1_ECAL)*Hsn(InhThr_ZEB1_ECAL,ZEB1,InhFld_ZEB1_ECAL,Hill_ZEB1_ECAL) - Deg_ECAL*ECAL
	du[7] = Prod_ZEB1*Hsn(InhThr_MIR200_ZEB1,MIR200,InhFld_MIR200_ZEB1,Hill_MIR200_ZEB1)*Hsn(InhThr_OVOL2_ZEB1,OVOL2,InhFld_OVOL2_ZEB1,Hill_OVOL2_ZEB1)*Hsp(ActThr_SNAIL1_ZEB1,SNAIL1,ActFld_SNAIL1_ZEB1,Hill_SNAIL1_ZEB1) - Deg_ZEB1*ZEB1
	du[8] = Prod_SNAIL1*Hsn(InhThr_MIR34A_SNAIL1,MIR34A,InhFld_MIR34A_SNAIL1,Hill_MIR34A_SNAIL1)*Hsn(InhThr_OVOL2_SNAIL1,OVOL2,InhFld_OVOL2_SNAIL1,Hill_OVOL2_SNAIL1)*Hsn(InhThr_SNAIL1_SNAIL1,SNAIL1,InhFld_SNAIL1_SNAIL1,Hill_SNAIL1_SNAIL1)*Hsp(ActThr_TGFB_SNAIL1,TGFB,ActFld_TGFB_SNAIL1,Hill_TGFB_SNAIL1) - Deg_SNAIL1*SNAIL1
end
