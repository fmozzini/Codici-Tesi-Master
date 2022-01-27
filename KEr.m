function KEr=KEr(W,Ht,x,y,z)
% Ht: weigth
% W: weigth 
% Piano saggitale(y) -> Imx = 3.44*Ht.^2 + 0.144*W-8.04;
% Piano frontale/coronale(z) -> Imy =  3.52*Ht.^2 + 0.125*W-7.78;
% Piano trasversale (?) 
% Imx = Izz
% Imy = Ixx, Iyy 
% Ixx =  3.52*Ht.^2 + 0.125*W-7.78;
% Iyy =  3.52*Ht.^2 + 0.125*W-7.78;
% Izz = 3.44*Ht.^2 + 0.144*W-8.04;
 %Nuovo
Ixx = 3.44*Ht.^2 + 0.144*W-8.04;
Iyy = 3.52*Ht.^2 + 0.125*W-7.78;
Izz = 3.44*Ht.^2 + 0.144*W-8.04;
KEr = 0.5*(Ixx*x.^2 + Iyy*y.^2 + Izz*z.^2);
end 
