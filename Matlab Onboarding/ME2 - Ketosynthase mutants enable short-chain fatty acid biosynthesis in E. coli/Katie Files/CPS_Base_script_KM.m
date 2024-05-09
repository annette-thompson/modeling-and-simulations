%All units are uM and sec.

%Specify initial conditions
%vector of initial conditions (units of uM)
init_cond = zeros(336,1);

init_cond(3) = 500;%s3 (Acetyl-CoA)
init_cond(4) = 10;%s6 (holo ACP)
init_cond(5) = 1000;%s7 (NADPH)
init_cond(6) = 1000;%s8 (NADH)
init_cond(8) = 500;%p2 (malonyl-CoA)

% init_cond(3) = 300;%s3 (Acetyl-CoA)
% init_cond(4) = 10;%s6 (holo ACP)
% init_cond(5) = 1300;%s7 (NADPH)
% init_cond(6) = 1300;%s8 (NADH)
% init_cond(8) = 1500;%p2 (malonyl-CoA)

compart=1E-3*6.7e-16;
init_cond(309) = 9.61452e-19/compart;%s3 (FtsH)
init_cond(310) = 6.39307e-19/compart;%s6 (LpxC)
init_cond(311) = 2.54062e-19/compart;%s8 (KdtA)

%Specify enzyme concentrations
%vector of enzyme concentrations (units of uM)
%ORDER: (  ACC FabD FabH FabG FabZ FabI TesA FabF FabA FabB)

%Specify time range to solve system (seconds)
time_range = [0 720];

conc=xlsread('FabB and FabF modelling conc')

FabB=conc(:,1);
FabF=conc(:,2);
%%

%Specify enzyme concentrations
%vector of enzyme concentrations (units of uM)
%ORDER: (  ACC FabD FabH FabG FabZ FabI TesA FabF FabA FabB)
enz_conc = [0  1    1    1    1    1    10    1    1    1];
%enz_conc = [0  1    1    1    2.59018648976186	1   10   0.497550254250183	0.0394246097521289	1.52877593179746];



conc_matrix=zeros(23,10);
for i=1:23
    for j=1:10
        conc_matrix(i,j)=enz_conc(1,j);
    end
end

%fill in FabB and FabF values
for i=1:23
        
        conc_matrix(i,10)=FabB(i);
        conc_matrix(i,8)=FabF(i);
end

%%
%Specify time range to solve system (seconds)
Frac_CL=zeros(23,14);
C16_tot=zeros(23,1);
FA_tot=zeros(23,1);
chain=zeros(23,14);
% RateC16=zeros(17,1);
% RateFA=zeros(17,1);

CL=[4,6,8,10,12,12,14,14,16,16,18,18,20,20];

for i=1:23
   
        [total_FA,F_weighted,frac_FA] = CPS_Base_solv(init_cond,conc_matrix(i,:),time_range);
%   F_weighted=palm_equiv: vector of palmitic acid equivalents at each time point
%   in time vals.
%   total_FA = sum(F_raw_all(end,:));%total fatty acid (concentration)
%frac_FA = F_raw_all(end,:)./total_FA;%fraction of each chain length fatty acid.
        
% for k=1:length(time_vals)
%             if time_vals(k)>=148 && time_vals(k-1)<148
%               index=k;
%             end
%         end
      

        C16_tot(i,1)=F_weighted(end)  %Find C16 equivalents at 12 min
        FA_tot(i,1)=total_FA; %Calculate total uM FA at 12 min
        Frac_CL(i,:)=frac_FA;
         for j=1:14
        chain(i,j) = frac_FA(j)*CL(j);
         end
         avg_CL(i,1)=sum(chain(i,:));
         
 
end
