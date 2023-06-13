#Defining custom Shifted Hill Equations
#Shifted Hill Equation for -ve regulation
Hsn(Thr,N,Fld,Hill) = ((Fld) + (1-Fld)*(1/(1 + (N/Thr)^Hill)))
#Shifted Hill Equation for +ve regulation
Hsp(Thr,N,Fld,Hill) = ((Fld) + (1-Fld)*(1/(1 + (N/Thr)^Hill)))/Fld

#Defining the reaction network with ModellingToolkit.jl macro
function rn!(du, u, p, t)
	B, A = u
	Prod_B, Prod_A ,Deg_B, Deg_A ,InhThr_A_B, InhFld_A_B, Hill_A_B, InhThr_B_A, InhFld_B_A, Hill_B_A = p
	du[1] = Prod_B*Hsn(InhThr_A_B,A,InhFld_A_B,Hill_A_B) - Deg_B*B
	du[2] = Prod_A*Hsn(InhThr_B_A,B,InhFld_B_A,Hill_B_A) - Deg_A*A
end
