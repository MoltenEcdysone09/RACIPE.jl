#Defining custom Shifted Hill Equations
#Shifted Hill Equation for -ve regulation
Hsn(thr,N,lbd,cop) = ((1/lbd) + (1-(1/lbd))*(1/(1 + (N/thr)^cop)))
#Shifted Hill Equation for +ve regulation
Hsp(thr,N,lbd,cop) = ((lbd) + (1-lbd)*(1/(1 + (N/thr)^cop)))/lbd

#Defining the reaction network with ModellingToolkit.jl macro
function rn!(du, u, p, t)
	B, A = u
	gB, gA ,kB, kA ,tAB, fAB, nAB, tBA, fBA, nBA = p
	du[1] = gB*Hsn(tAB,A,fAB,nAB) - kB*B
	du[2] = gA*Hsn(tBA,B,fBA,nBA) - kA*A
end
