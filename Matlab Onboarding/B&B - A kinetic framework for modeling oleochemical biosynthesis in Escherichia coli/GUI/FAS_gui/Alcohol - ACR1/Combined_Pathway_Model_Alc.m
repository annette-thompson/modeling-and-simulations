function dy = Combined_Pathway_Model_Alc(t,y,param_list,opt_struct,p_vec_k)
%Contains all the differential equations and enzyme balances that define the FAS model
%   Input:
%       t: time values (required as input for MATLAB ODE solver, sec)
%       y: concentration values (all components and intermediates, uM)
%       param_list: structure containing all kinetic parameters
%       opt_struct: structure containining additional parameters of
%       model(such as initial enzyme concentration)
%   Output:
%       dy: values of differential equations for given 
%       concentrations and kinetic parameters

%Parameters (defined from values in param_list)

%ACC (not used)
k1_1f = param_list.k1_1f;
k1_1r = param_list.k1_1r;
k1_2f = param_list.k1_2f;
k1_2r = param_list.k1_2r;
kcat1_1 = param_list.kcat1_1;
k1_3f = param_list.k1_3f;
k1_3r = param_list.k1_3r;
kcat1_2 = param_list.kcat1_2;
e1tot = opt_struct.enzyme_conc(1);

%FabD
k2_1f = param_list.k2_1f;
k2_1r = param_list.k2_1r;
k2_2f = param_list.k2_2f;
k2_2r = param_list.k2_2r;
k2_3f = param_list.k2_3f;
k2_3r = param_list.k2_3r;
k2_4f = param_list.k2_4f;
k2_4r = param_list.k2_4r;
e2tot = opt_struct.enzyme_conc(2);

%FabH
k3_1f = param_list.k3_1f;
k3_1r = param_list.k3_1r;
k3_2f = param_list.k3_2f;
k3_2r = param_list.k3_2r;
k3_3f = param_list.k3_3f;
k3_3r = param_list.k3_3r;
k3_4f = param_list.k3_4f;
k3_4r = param_list.k3_4r;
k3_5f = param_list.k3_5f;
k3_5r = param_list.k3_5r;
kcat3 = param_list.kcat3;
k3_inh_f = opt_struct.ACP_inh(1);
k3_inh_r = opt_struct.ACP_inh(2);
e3tot = opt_struct.enzyme_conc(3);

%FabG
k4_1f = param_list.k4_1f;
k4_1r = param_list.k4_1r;
k4_2f = param_list.k4_2f;
k4_2r = param_list.k4_2r;
kcat4 = param_list.kcat4;
k4_inh_f = opt_struct.ACP_inh(3);
k4_inh_r = opt_struct.ACP_inh(4);
e4tot = opt_struct.enzyme_conc(4);

%FabZ
k5_1f = param_list.k5_1f;
k5_1r = param_list.k5_1r;
kcat5 = param_list.kcat5;
k5_2r = param_list.k5_2r;
k5_3f = k5_1r;
k5_3r = k5_1f;
kcat9_est = param_list.kcat9;
kcat9_un_est = kcat9_est.*[1,1,1,opt_struct.scaling_factor_FabA_unsat,1,1,1,1,1];
ratio_ZA =  kcat9_un_est(4)/kcat5(4);
rel_ZA = 0.17946003;
kcat5_un = kcat5.*[1,1,1,rel_ZA*ratio_ZA,1,1,1,1,1];
k5_2r_un = k5_2r.*[1,1,1,rel_ZA*ratio_ZA,1,1,1,1,1];
k5_3f_un = k5_3f;
k5_3r_un = k5_3r;
k5_inh_f = opt_struct.ACP_inh(5);
k5_inh_r = opt_struct.ACP_inh(6);
e5tot = opt_struct.enzyme_conc(5);

%FabI
k6_1f = param_list.k6_1f;
k6_1r = param_list.k6_1r;
k6_2f = param_list.k6_2f;
k6_2r = param_list.k6_2r;
kcat6 = param_list.kcat6;
%For turning off C20 production (not seen in vitro or in vivo)
kcat6(9) = 0;
k6_inh_f = opt_struct.ACP_inh(7);
k6_inh_r = opt_struct.ACP_inh(8);
e6tot = opt_struct.enzyme_conc(6);

%TesA
thio_scaling = [0.454318854716137	2.19560691520507	68.4348965082940	12.9343678377699];
if strcmp(opt_struct.TesA_fitting_source,'TesA') || strcmp(opt_struct.TesA_fitting_source,'Pf')
    k7_1f = param_list.k7_1f;
elseif strcmp(opt_struct.TesA_fitting_source,'CpFatB')
    k7_1f = thio_scaling(3).*param_list.k7_1f;
elseif strcmp(opt_struct.TesA_fitting_source,'BTE')
    k7_1f = thio_scaling(4).*param_list.k7_1f;
else
    k7_1f = 0;
end
k7_1r = param_list.k7_1r;
if strcmp(opt_struct.TesA_fitting_source,'TesA') || strcmp(opt_struct.TesA_fitting_source,'Pf')
    kcat7 = param_list.kcat7;
elseif strcmp(opt_struct.TesA_fitting_source,'CpFatB')
    kcat7 = thio_scaling(1).*param_list.kcat7;

elseif strcmp(opt_struct.TesA_fitting_source, 'BTE')
    kcat7 = thio_scaling(2).*param_list.kcat7;
else
    kcat7 = 0;
end 
k7_inh_f = opt_struct.ACP_inh(9);
k7_inh_r = opt_struct.ACP_inh(10);
e7tot = opt_struct.enzyme_conc(7);

%FabF
k8_1f = param_list.k8_1f;
k8_1r = param_list.k8_1r;
k8_2f = param_list.k8_2f;
k8_2r = param_list.k8_2r;
k8_3f = param_list.k8_3f;
k8_3r = param_list.k8_3r;
k8_4f = param_list.k8_4f;
k8_4r = param_list.k8_4r;
k8_5f = param_list.k8_5f;
k8_5r = param_list.k8_5r;
k8_6f = param_list.k8_6f;
k8_6r = param_list.k8_6r;
k8_7f = param_list.k8_7f;
k8_7r = param_list.k8_7r;
k8_8f = param_list.k8_8f;
k8_8r = param_list.k8_8r;
k8_9f = param_list.k8_9f;
k8_9r = param_list.k8_9r;
kcat8 = param_list.kcat8;
kcat8_H = param_list.kcat8_H;
kcat8_CO2 = param_list.kcat8_CO2;
kcat10_est = param_list.kcat10;
ratio_FB = kcat10_est(5)/kcat8(4);
rel_FB = 0.123650181;
%rel_FB = opt_struct.rel_FB;
%kcat8_un = kcat8(4).*[1,1,1,rel_FB*ratio_FB,0.9,0.289,0.34,0.34,0.34];%specificity of reaction with unsaturated acyl chains
kcat8_un = kcat8(4).*[1,1,1,rel_FB*ratio_FB,0.9,0.289,0.34,0,0];%specificity of reaction with unsaturated acyl chains, no C20
k8_inh_f = opt_struct.ACP_inh(11);
k8_inh_r = opt_struct.ACP_inh(12);
e8tot = opt_struct.enzyme_conc(8);

%FabA
k9_1f = param_list.k9_1f;
k9_1r = param_list.k9_1r;
kcat9 = param_list.kcat9;
k9_2r = param_list.k9_2r;
k9_3f = k9_1r;
k9_3r = k9_1f;
k9_1f_un = k9_1f.*[1,1,1,1,0.0358,0.0358,0.0358,0.0358,0.0358];%specificity of reaction with unsaturated acyl chains
k9_1r_un = k9_1r;
kcat9_un = kcat9.*[1,1,1,opt_struct.scaling_factor_FabA_unsat,1,1,1,1,1];
k9_2r_un = k9_2r;
k9_3f_un = k9_3f;
k9_3r_un = k9_3r;
k9_inh_f = opt_struct.ACP_inh(13);
k9_inh_r = opt_struct.ACP_inh(14);
e9tot = opt_struct.enzyme_conc(9);

%FabB
k10_1f = param_list.k10_1f;
k10_1r = param_list.k10_1r;
k10_2f = param_list.k10_2f;
k10_2r = param_list.k10_2r;
k10_3f = param_list.k10_3f;
k10_3r = param_list.k10_3r;
k10_4f = param_list.k10_4f;
k10_4r = param_list.k10_4r;
k10_5f = param_list.k10_5f;
k10_5r = param_list.k10_5r;
k10_6f = param_list.k10_6f;
k10_6r = param_list.k10_6r;
k10_7f = param_list.k10_7f;
k10_7r = param_list.k10_7r;
k10_8f = param_list.k10_8f;
k10_8r = param_list.k10_8r;
k10_9f = param_list.k10_9f;
k10_9r = param_list.k10_9r;
kcat10 = param_list.kcat10;
kcat10_H = param_list.kcat10_H;
kcat10_CO2 = param_list.kcat10_CO2;
%kcat10_un = kcat10(5).*[1,1,1,1,1,0.125,0.0229,0.0229,0.0229];%specificity of reaction with unsaturated acyl chains
kcat10_un = kcat10(5).*[1,1,1,1,1,0.125,0.0229,0,0];%specificity of reaction with unsaturated acyl chains, no C20
k10_inh_f = opt_struct.ACP_inh(15);
k10_inh_r = opt_struct.ACP_inh(16);
e10tot = opt_struct.enzyme_conc(10);


%reaction 'LpxA_forward':  kinetic parameter 'Kma'
Kma=0.82*1E3;
%reaction 'LpxA_forward':  kinetic parameter 'Kmb'
Kmb=0.0016*1E3;
%reaction 'LpxA_forward':  kinetic parameter 'A_kcat'
A_kcat=0.007941921;%7.17;
%reaction 'LpxA_reverse':  kinetic parameter 'Kma'
Kma1=0.82*1E3;
%reaction 'LpxA_reverse':  kinetic parameter 'Kmb'
Kmb1=0.0016*1E3;
%reaction 'LpxA_reverse':  kinetic parameter 'A_kcat'
A_kcat1=717;
%reaction 'y(2,:)':  kinetic parameter 'A_kcat'
A_kcat2=3.3;
%reaction 'y(2,:)':  kinetic parameter 'Km'
Km=0.00019*1E3;
%reaction 'LpxC_translate':  kinetic parameter 'k1'
k1=0.148;
%reaction 'LpxC_degrade':  kinetic parameter 'k1'
k11=9.62e-05;
%reaction 'LpxH':  kinetic parameter 'Km'
Km1=0.0617*1E3;
%reaction 'LpxH':  kinetic parameter 'Ki'
Ki=0.015*1E3;
%reaction 'LpxH':  kinetic parameter 'A_kcat'
A_kcat3=47;
%reaction 'LpxB':  kinetic parameter 'Kma'
Kma2=0.287*1E3;
%reaction 'LpxB':  kinetic parameter 'Kmb'
Kmb2=0.381*1E3;
%reaction 'LpxB':  kinetic parameter 'A_kcat'
A_kcat4=129;
%reaction 'LpxK':  kinetic parameter 'A_kcat'
A_kcat5=2.1;
%reaction 'LpxK':  kinetic parameter 'Km'
Km2=0.04*1E3;
%reaction 'KdtA_translate':  kinetic parameter 'k1'
k12=0.176;
%reaction 'KdtA_degrade':  kinetic parameter 'k1'
k13=0.00018;
%reaction 'y(4,:)':  kinetic parameter 'A_kcat'
A_kcat6=16.7;
%reaction 'y(4,:)':  kinetic parameter 'Kma'
Kma3=0.088*1E3;
%reaction 'y(4,:)':  kinetic parameter 'Kmb'
Kmb3=0.052*1E3;
%reaction 'y(4,:)':  kinetic parameter 'ki'
ki1=0.0317*1E3;
%reaction 'LpxL':  kinetic parameter 'A_kcat'
A_kcat7=131;
%reaction 'LpxL':  kinetic parameter 'Km'
Km3=0.015*1E3;
%reaction 'LpxM':  kinetic parameter 'A_kcat'
A_kcat8=0.6;
%reaction 'LpxM':  kinetic parameter 'Km'
Km4=0.00275*1E3;
%reaction 'MsbA':  kinetic parameter 'A_kcat'
A_kcat9=166;
%reaction 'MsbA':  kinetic parameter 'Km'
Km5=0.021*1E3;
%reaction 'LpxC_proteolysis':  kinetic parameter 'k1'
k14=2.015*1E-3;
%reaction 'FtsH_activate_LpxC':  kinetic parameter 'k1'
k15=0.14*1E-3;
%reaction 'FtsH_activate_KdtA':  kinetic parameter 'k1'
k16=32.3*1E-3;
%reaction 'KdtA_proteolysis':  kinetic parameter 'k1'
k17=6.8*1E-3;
%reaction 'LpxD':  kinetic parameter 'A_kcat'
A_kcat10=23;
%reaction 'LpxD':  kinetic parameter 'Kma'
Kma4=0.0025*1E3;
%reaction 'LpxD':  kinetic parameter 'Kmb'
Kmb4=0.0032*1E3;
%reaction 'LpxD':  kinetic parameter 'ki'
ki2=0.0094*1E3;
%reaction 'KdtA_alternate':  kinetic parameter 'A_kcat'
A_kcat11=1.9;
%reaction 'KdtA_alternate':  kinetic parameter 'Kma'
Kma5=0.088*1E3;
%reaction 'KdtA_alternate':  kinetic parameter 'Kmb'
Kmb5=0.052*1E3;
%reaction 'FtsH_inactive_LpxC':  kinetic parameter 'k1'
k18=0.1;
%reaction 'FtsH_inactive_KdtA':  kinetic parameter 'k1'
k19=0.1;

km_laur = 0.0032*1E3;
km_myr = 0.0032*1E3;
Laur_acp = 2E-2;
Myr_acp = 2E-2;

%LipidA_multiplier (concentration)
lipidA_mult = 0;

compart=1E-3*6.7e-16;

% Initial values:
%metabolite 'UDP-GlcNac': ode
%UDP_Glc=3.32108e-15/compart;
UDP_Glc = 2.881E-16/compart;%430uM
%metabolite 'Beta-hydroxymyristoylACP': ode
Beta_hy=3.32108e-15/compart;
%metabolite 'ACP': ode
ACP=1.70039e-18/compart;
%metabolite 'CMP-KDO': ode
CMP_KDO=3.32108e-15/compart;
%metabolite 'LpxA': reactions
LpxA=lipidA_mult*1.1026e-18/compart;
%metabolite 'promoter-LpxC': reactions
promote=1.66054e-21/compart;
%metabolite 'LpxD': reactions
LpxD=lipidA_mult*7.52224e-19/compart;
%metabolite 'LpxH': reactions
LpxH=lipidA_mult*2.93915e-19/compart;
%metabolite 'LpxB': reactions
LpxB=lipidA_mult*6.37647e-19/compart;
%metabolite 'LpxK': reactions
LpxK=lipidA_mult*7.17353e-19/compart;
%metabolite 'promoter-KdtA': reactions
promot=1.66054e-21/compart;
%metabolite 'LpxL': reactions
LpxL=lipidA_mult*1.54098e-18/compart;
%metabolite 'LpxM': reactions
LpxM=lipidA_mult*6.1772e-18/compart;
%metabolite 'MsbA': reactions
MsbA=lipidA_mult*3.42071e-19/compart;

FtsH_init = lipidA_mult*9.6145195545121833e-19/compart;

%Phospholipid_multiplier (concentration)
phos_mult = 0;

%Phospholipid Constants
kcatP1 = 3874.446011;%10;%PlsB
kcatP2 = 3874.446011;%10;%PlsC
KmP1 = 15;%PlsB C16-ACP
KmP2 = 140;%PlsB G3P
KmP3 = 25;%PlsB C18:1-ACP
KmP4 = 15;%PlsC C16:1-ACP
KmP5 = 140;%PlsC LysP
KmP6 = 25;%PlsC C18:1-ACP

G3P = 0.18E3;
PlsB = phos_mult*0.3569;
PlsC = phos_mult*0.3569;


%Alcohol production kinetic constants
%FadD
Km_ATP = 36; %doi.org/10.1016/0167-4838(85)90269-9
Km_CoA = 40.67; %doi.org/10.1016/0167-4838(85)90269-9
Km_FA_AMP = [11.2  40.8  4.3 1.6  2.4  2.4  8.9  1.6  2.4  2.4  8.5]; %doi.org/10.1016/S0021-9258(19)69262-8
Km_FA = [11.2  40.8  4.3 1.6  2.4  2.4  8.9  1.6  2.4  2.4  8.5]; %doi.org/10.1016/S0021-9258(19)69262-8
kcatD1 = 11.7210624038738*0.029.*[0.103 0.132 0.299 1.000 0.798 0.704 0.434 1.000 0.798 0.704 0.519]; %10.1042/0264-6021:3600699 AND Promiscuous Fatty Acyl-CoA Ligases Produce Acyl-CoA and Acyl-SNAC Precursors for Polyketide Biosynthesis
kcatD2 = 60.0916342073819*10*0.029; %10.1042/0264-6021:3600699

%ACR1
Km_ACR1 = 32; %https://www.brenda-enzymes.org/enzyme.php?ecno=1.2.1.B25#TURNOVER%20NUMBER%20[1/s]
kcatACR1 = 0.006*[0  0  1*p_vec_k(2)  1*p_vec_k(2)  1*p_vec_k(3)  1*p_vec_k(4)  0  1*p_vec_k(2)  1*p_vec_k(3)  1*p_vec_k(4)  0]; %https://www.brenda-enzymes.org/enzyme.php?ecno=1.2.1.B25#TURNOVER%20NUMBER%20[1/s]46

%Ahr
kcatAhr = [464.4 516.5948127 246.2524496 5.353314121 8.029971182 8.029971182 8.029971182 8.029971182 8.029971182 8.029971182 8.029971182 8.029971182]; %DOI 10.1007/s00253-012-4474-5
Km_Ahr_NADPH = 60;
Km_FA_Ald = [2100 340 340 340 340 340 340 340 340 340 340 340];

%Enzyme concentrations (uM)
FadD = opt_struct.aux_enz(1);
ACR1 = opt_struct.aux_enz(2);
Ahr = opt_struct.aux_enz(3);


%Enzyme concentration balance equations
e1 = e1tot - y(83,:) - y(84,:) - y(85,:) - y(86,:);%ACC
e2 = e2tot - y(87,:) - y(88,:) - y(89,:);%FabD
e3 = e3tot - y(90,:) - y(91,:) - y(92,:) - y(95,:) - y(110,:) - y(125,:) - y(140,:) - y(155,:) - y(170,:) - y(185,:) - y(200,:) - y(215,:)...
    - y(96,:) - y(111,:) - y(126,:) - y(141,:) - y(156,:) - y(171,:) - y(186,:) - y(201,:) - y(216,:) - y(293,:);%FabH
e4 = e4tot - y(93,:) - y(97,:) - y(112,:) - y(127,:) - y(142,:) - y(157,:) - y(172,:) - y(187,:)...
    - y(202,:) - y(217,:) - y(226,:) - y(241,:) - y(256,:) - y(271,:) - y(286,:) - y(294,:);%FabG
e5 = e5tot - y(98,:) - y(113,:) - y(128,:) - y(143,:) - y(158,:) - y(173,:) - y(188,:) - y(203,:) - y(218,:)...
    - y(227,:) - y(242,:) - y(257,:) - y(272,:) - y(287,:)...
    - y(99,:) - y(114,:) - y(129,:) - y(144,:) - y(159,:) - y(174,:) - y(189,:) - y(204,:) - y(219,:)...
    - y(228,:) - y(243,:) - y(258,:) - y(273,:) - y(288,:)...
    - y(295,:) - y(305,:);%FabZ
e6 = e6tot - y(94,:) - y(100,:) - y(115,:) - y(130,:) - y(145,:) - y(160,:) - y(175,:) - y(190,:)...
    - y(205,:) - y(220,:) - y(229,:) - y(244,:) - y(259,:) - y(274,:) - y(289,:) - y(296,:);%FabI
e7 = e7tot - y(101,:) - y(116,:) - y(131,:) - y(146,:) - y(161,:) - y(176,:) - y(191,:) - y(206,:) - y(221,:)...
    - y(230,:) - y(245,:) - y(260,:) - y(275,:) - y(290,:) - y(297,:);%TesA
e8 = e8tot - y(102,:) - y(103,:) - y(104,:) - y(117,:) - y(118,:) - y(119,:) - y(132,:) - y(133,:) - y(134,:)...
     - y(147,:) - y(148,:) - y(149,:) - y(162,:) - y(231,:) - y(163,:) - y(232,:) - y(164,:) - y(233,:)...
     - y(177,:) - y(246,:) - y(178,:) - y(247,:) - y(179,:) - y(248,:)...
     - y(192,:) - y(261,:) - y(193,:) - y(262,:) - y(194,:) - y(263,:)...
     - y(207,:) - y(276,:) - y(208,:) - y(277,:) - y(209,:) - y(278,:) - y(298,:)...
     - y(306,:) - y(307,:) - y(308,:)...
     - y(327,:) - y(328,:) - y(329,:) - y(330,:) - y(331,:);%FabF
e9 = e9tot - y(105,:) - y(120,:) - y(135,:) - y(150,:) - y(165,:) - y(180,:) - y(195,:) - y(210,:) - y(222,:)...
    - y(234,:) - y(249,:) - y(264,:) - y(279,:) - y(291,:)...
    - y(106,:) - y(121,:) - y(136,:) - y(151,:) - y(166,:) - y(181,:) - y(196,:) - y(211,:) - y(223,:)...
    - y(235,:) - y(250,:) - y(265,:) - y(280,:) - y(292,:) - y(301,:)...
    - y(299,:);%FabA
e10 = e10tot - y(107,:) - y(108,:) - y(109,:) - y(122,:) - y(123,:) - y(124,:) - y(137,:) - y(138,:) - y(139,:)...
     - y(152,:) - y(153,:) - y(154,:) - y(167,:) - y(236,:) - y(168,:) - y(237,:) - y(169,:) - y(238,:)...
     - y(182,:) - y(251,:) - y(183,:) - y(252,:) - y(184,:) - y(253,:)...
     - y(197,:) - y(266,:) - y(198,:) - y(267,:) - y(199,:) - y(268,:)...
     - y(212,:) - y(281,:) - y(213,:) - y(282,:) - y(214,:) - y(283,:)...
     - y(302,:) - y(303,:) - y(304,:) - y(300,:)...
     - y(332,:) - y(333,:) - y(334,:) - y(335,:) -y(336,:);%FabB


dy = zeros(380,size(y,2));

%set of differential equations
dy(1,:) = k1_1r.*y(83,:) - k1_1f.*e1.*y(1,:)...
    - kcatD1(1).*FadD.*y(1,:).*y(21,:)./((Km_ATP+y(1,:)).*(Km_FA(1)+y(21,:)))...
    - kcatD1(2).*FadD.*y(1,:).*y(26,:)./((Km_ATP+y(1,:)).*(Km_FA(2)+y(26,:))) ...
    - kcatD1(3).*FadD.*y(1,:).*y(31,:)./((Km_ATP+y(1,:)).*(Km_FA(3)+y(31,:)))...
    - kcatD1(4).*FadD.*y(1,:).*y(37,:)./((Km_ATP+y(1,:)).*(Km_FA(4)+y(37,:)))...
    - kcatD1(5).*FadD.*y(1,:).*y(47,:)./((Km_ATP+y(1,:)).*(Km_FA(5)+y(47,:)))...
    - kcatD1(6).*FadD.*y(1,:).*y(57,:)./((Km_ATP+y(1,:)).*(Km_FA(6)+y(57,:)))...
    - kcatD1(7).*FadD.*y(1,:).*y(67,:)./((Km_ATP+y(1,:)).*(Km_FA(7)+y(67,:)))...
    - kcatD1(8).*FadD.*y(1,:).*y(42,:)./((Km_ATP+y(1,:)).*(Km_FA(8)+y(42,:)))...
    - kcatD1(9).*FadD.*y(1,:).*y(52,:)./((Km_ATP+y(1,:)).*(Km_FA(9)+y(52,:)))...
    - kcatD1(10).*FadD.*y(1,:).*y(62,:)./((Km_ATP+y(1,:)).*(Km_FA(10)+y(62,:)))...
    - kcatD1(11).*FadD.*y(1,:).*y(72,:)./((Km_ATP+y(1,:)).*(Km_FA(11)+y(72,:)));

dy(2,:) = k1_2r.*y(84,:) - k1_2f.*y(83,:).*y(2,:);
dy(3,:) = k1_3r.*y(86,:) + k3_1r.*y(90,:) - k1_3f.*y(85,:).*y(3,:) - k3_1f.*e3.*y(3,:)...
    - k8_4f.*e8.*y(3,:) + k8_4r.*y(330,:) - k10_4f.*e10.*y(3,:) + k10_4r.*y(335,:);
dy(4,:) = k2_3r.*y(89,:) - k2_3f.*y(88,:).*y(4,:)...
    + kcat7(1).*y(101,:)...
    + kcat7(2).*y(116,:)...
    + kcat7(3).*y(131,:)...
    + kcat7(4).*y(146,:)...
    + kcat7(5).*y(161,:)...
    + kcat7(6).*y(176,:)...
    + kcat7(7).*y(191,:)...
    + kcat7(8).*y(206,:)...
    + kcat7(9).*y(221,:)...
    + kcat7(5).*y(230,:)...
    + kcat7(6).*y(245,:)...
    + kcat7(7).*y(260,:)...
    + kcat7(8).*y(275,:)...
    + kcat7(9).*y(290,:)...
    + k8_2f(1).*y(102,:) - k8_2r(1).*y(103,:).*y(4,:)...
    + k8_2f(2).*y(117,:) - k8_2r(2).*y(118,:).*y(4,:)...
    + k8_2f(3).*y(132,:) - k8_2r(3).*y(133,:).*y(4,:)...
    + k8_2f(4).*y(147,:)- k8_2r(4).*y(148,:).*y(4,:)...
    + k8_2f(5).*y(162,:)- k8_2r(5).*y(163,:).*y(4,:)...
    + k8_2f(6).*y(177,:)- k8_2r(6).*y(178,:).*y(4,:)...
    + k8_2f(7).*y(192,:)- k8_2r(7).*y(193,:).*y(4,:)...
    + k8_2f(8).*y(207,:)- k8_2r(8).*y(208,:).*y(4,:)...
    + k8_2f(5).*y(231,:)- k8_2r(5).*y(232,:).*y(4,:)...
    + k8_2f(6).*y(246,:)- k8_2r(6).*y(247,:).*y(4,:)...
    + k8_2f(7).*y(261,:)- k8_2r(7).*y(262,:).*y(4,:)...
    + k8_2f(8).*y(276,:)- k8_2r(8).*y(277,:).*y(4,:)...
    + k10_2f(1).*y(107,:) - k10_2r(1).*y(108,:).*y(4,:)...
    + k10_2f(2).*y(122,:) - k10_2r(2).*y(123,:).*y(4,:)...
    + k10_2f(3).*y(137,:) - k10_2r(3).*y(138,:).*y(4,:)...
    + k10_2f(4).*y(152,:)- k10_2r(4).*y(153,:).*y(4,:)...
    + k10_2f(5).*y(167,:)- k10_2r(5).*y(168,:).*y(4,:)...
    + k10_2f(6).*y(182,:)- k10_2r(6).*y(183,:).*y(4,:)...
    + k10_2f(7).*y(197,:)- k10_2r(7).*y(198,:).*y(4,:)...
    + k10_2f(8).*y(212,:)- k10_2r(8).*y(213,:).*y(4,:)...
    + k10_2f(5).*y(236,:)- k10_2r(5).*y(237,:).*y(4,:)...
    + k10_2f(6).*y(251,:)- k10_2r(6).*y(252,:).*y(4,:)...
    + k10_2f(7).*y(266,:)- k10_2r(7).*y(267,:).*y(4,:)...
    + k10_2f(8).*y(281,:)- k10_2r(8).*y(282,:).*y(4,:)...
    + k10_2f(4).*y(302,:)- k10_2r(4).*y(303,:).*y(4,:)...
    + k8_2f(4).*y(306,:) - k8_2r(4).*y(307,:).*y(4,:)...
    + k3_inh_r.*y(293,:)  - k3_inh_f.*e3.*y(4,:)...
    + k4_inh_r.*y(294,:)  - k4_inh_f.*e4.*y(4,:)...
    + k5_inh_r.*y(295,:)  - k5_inh_f.*e5.*y(4,:)...
    + k6_inh_r.*y(296,:)  - k6_inh_f.*e6.*y(4,:)...
    + k7_inh_r.*y(297,:)  - k7_inh_f.*e7.*y(4,:)...
    + k8_inh_r.*y(298,:)  - k8_inh_f.*e8.*y(4,:)...
    + k9_inh_r.*y(299,:)  - k9_inh_f.*e9.*y(4,:)...
    + k10_inh_r.*y(300,:)  - k10_inh_f.*e10.*y(4,:)...
    + LpxD.*A_kcat10./(1+y(311,:)./ki2).*y(318,:).*y(44,:)./(Kma4.*Kmb4+y(318,:).*Kmb4+y(44,:).*Kma4+y(318,:).*y(44,:))...
    + LpxA.*A_kcat.*UDP_Glc.*y(44,:)./(Kma.*Kmb+UDP_Glc.*Kmb+y(44,:).*Kma+UDP_Glc.*y(44,:))...
    + A_kcat7.*LpxL.*y(320,:).*y(36,:)./((Km3+y(320,:).*(km_laur+y(36,:))))...
    + A_kcat8.*LpxM.*y(317,:).*y(46,:)./((Km4+y(317,:).*(km_myr+y(46,:))))...
    - LpxA.*A_kcat1.*y(314,:).*y(4,:)./(Kma1.*Kmb1+y(314,:).*Kmb1+y(4,:).*Kma1+y(314,:).*y(4,:))...
    + kcatP1.*PlsB.*G3P.*y(56,:)./((KmP1+G3P).*(KmP2+y(56,:)))...
    + kcatP2*PlsC.*y(324,:).*y(61,:)./((KmP4 + y(324,:)).*(KmP5 + y(61,:)))...
    + kcatP1.*PlsB.*G3P.*y(71,:)./((KmP1+G3P).*(KmP3+y(71,:)))...
    + kcatP2*PlsC.*y(324,:).*y(71,:)./((KmP4 + y(324,:)).*(KmP6 + y(71,:)))...
    + k8_9f.*y(328,:) - k8_9r.*y(329,:).*y(4,:) + k10_9f.*y(333,:) - k10_9r.*y(334,:).*y(4,:)...
    ;
dy(5,:) = k4_1r(1).*y(93,:) - k4_1f(1).*e4.*y(5,:);
dy(6,:) = k6_1r(1).*y(94,:) - k6_1f(1).*e6.*y(6,:)...
     - kcatAhr(2).*Ahr.*y(6,:).*y(359,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(2)+y(359,:)))...
     - kcatAhr(3).*Ahr.*y(6,:).*y(360,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(3)+y(360,:)))...
     - kcatAhr(4).*Ahr.*y(6,:).*y(361,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(4)+y(361,:)))...
     - kcatAhr(5).*Ahr.*y(6,:).*y(362,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(5)+y(362,:)))...
     - kcatAhr(6).*Ahr.*y(6,:).*y(363,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(6)+y(363,:)))...
     - kcatAhr(7).*Ahr.*y(6,:).*y(364,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(7)+y(364,:)))...
     - kcatAhr(8).*Ahr.*y(6,:).*y(365,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(8)+y(365,:)))...
     - kcatAhr(9).*Ahr.*y(6,:).*y(366,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(9)+y(366,:)))...
     - kcatAhr(10).*Ahr.*y(6,:).*y(367,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(10)+y(367,:)))...
     - kcatAhr(11).*Ahr.*y(6,:).*y(368,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(11)+y(368,:)))...
     - kcatAhr(12).*Ahr.*y(6,:).*y(369,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(12)+y(369,:)));
dy(7,:) = kcat1_2.*y(86,:);
dy(8,:) = kcat1_2.*y(86,:) + k2_1r.*y(87,:) - k2_1f.*e2.*y(8,:);
dy(9,:) = k2_2f.*y(87,:) + k3_2f.*y(90,:) - k2_2r.*y(88,:).*y(9,:) - k3_2r.*y(91,:).*y(9,:)...
    + k8_5f.*y(330,:) - k8_5r.*y(329,:).*y(9,:) + k10_5f.*y(335,:) - k10_5r.*y(334,:).*y(9,:)...
    - kcatD2.*FadD.*y(9,:).*y(337,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(1)+y(337,:)))...
    - kcatD2.*FadD.*y(9,:).*y(338,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(2)+y(338,:)))...
    - kcatD2.*FadD.*y(9,:).*y(339,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(3)+y(339,:)))...
    - kcatD2.*FadD.*y(9,:).*y(340,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(4)+y(340,:)))...
    - kcatD2.*FadD.*y(9,:).*y(341,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(5)+y(341,:)))...
    - kcatD2.*FadD.*y(9,:).*y(342,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(6)+y(342,:)))...
    - kcatD2.*FadD.*y(9,:).*y(343,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(7)+y(343,:)))...
    - kcatD2.*FadD.*y(9,:).*y(344,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(8)+y(344,:)))...
    - kcatD2.*FadD.*y(9,:).*y(345,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(9)+y(345,:)))...
    - kcatD2.*FadD.*y(9,:).*y(346,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(10)+y(346,:)))...
    - kcatD2.*FadD.*y(9,:).*y(347,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(11)+y(347,:)));

dy(10,:) = k2_4f.*y(89,:) + k3_3r.*y(92,:) - k2_4r.*e2.*y(10,:) - k3_3f.*y(91,:).*y(10,:)...
    + k8_3r(1).*y(104,:) - k8_3f(1).*y(103,:).*y(10,:)...
    + k8_3r(2).*y(119,:) - k8_3f(2).*y(118,:).*y(10,:)...
    + k8_3r(3).*y(134,:) - k8_3f(3).*y(133,:).*y(10,:)...
    + k8_3r(4).*y(149,:)- k8_3f(4).*y(148,:).*y(10,:)...
    + k8_3r(5).*y(164,:)- k8_3f(5).*y(163,:).*y(10,:)...
    + k8_3r(6).*y(179,:)- k8_3f(6).*y(178,:).*y(10,:)...
    + k8_3r(7).*y(194,:)- k8_3f(7).*y(193,:).*y(10,:)...
    + k8_3r(8).*y(209,:)- k8_3f(8).*y(208,:).*y(10,:)...
    + k8_3r(5).*y(233,:)- k8_3f(5).*y(232,:).*y(10,:)...
    + k8_3r(6).*y(248,:)- k8_3f(6).*y(247,:).*y(10,:)...
    + k8_3r(7).*y(263,:)- k8_3f(7).*y(262,:).*y(10,:)...
    + k8_3r(8).*y(278,:)- k8_3f(8).*y(277,:).*y(10,:)...
    + k10_3r(1).*y(109,:) - k10_3f(1).*y(108,:).*y(10,:)...
    + k10_3r(2).*y(124,:) - k10_3f(2).*y(123,:).*y(10,:)...
    + k10_3r(3).*y(139,:) - k10_3f(3).*y(138,:).*y(10,:)...
    + k10_3r(4).*y(154,:)- k10_3f(4).*y(153,:).*y(10,:)...
    + k10_3r(5).*y(169,:)- k10_3f(5).*y(168,:).*y(10,:)...
    + k10_3r(6).*y(184,:)- k10_3f(6).*y(183,:).*y(10,:)...
    + k10_3r(7).*y(199,:)- k10_3f(7).*y(198,:).*y(10,:)...
    + k10_3r(8).*y(214,:)- k10_3f(8).*y(213,:).*y(10,:)...
    + k10_3r(5).*y(238,:)- k10_3f(5).*y(237,:).*y(10,:)...
    + k10_3r(6).*y(253,:)- k10_3f(6).*y(252,:).*y(10,:)...
    + k10_3r(7).*y(268,:)- k10_3f(7).*y(267,:).*y(10,:)...
    + k10_3r(8).*y(283,:)- k10_3f(8).*y(282,:).*y(10,:)...
    + k10_3r(4).*y(304,:)- k10_3f(4).*y(303,:).*y(10,:)...
    + k8_3r(4).*y(308,:) - k8_3f(4).*y(307,:).*y(10,:)...
    - k8_7f.*e8.*y(10,:) + k8_7r.*y(327,:) - k8_6f.*y(329,:).*y(10,:) + k8_6r.*y(331,:)...
    - k10_7f.*e10.*y(10,:) + k10_7r.*y(332,:) - k10_6f.*y(334,:).*y(10,:) + k10_6r.*y(336,:)...
    ;
dy(11,:) = kcat3.*y(92,:) + kcat8(1).*y(104,:) + kcat8(2).*y(119,:) + kcat8(3).*y(134,:) + kcat8(4).*y(149,:)...
    + kcat8(5).*y(164,:) + kcat8(6).*y(179,:) + kcat8(7).*y(194,:) + kcat8(8).*y(209,:)...
    + kcat8_un(5).*y(233,:) + kcat8_un(6).*y(248,:) + kcat8_un(7).*y(263,:) + kcat8_un(8).*y(278,:)...
    + kcat10(1).*y(109,:) + kcat10(2).*y(124,:) + kcat10(3).*y(139,:) + kcat10(4).*y(154,:)...
    + kcat10(5).*y(169,:) + kcat10(6).*y(184,:) + kcat10(7).*y(199,:) + kcat10(8).*y(214,:)...
    + kcat10_un(4).*y(304,:) + kcat10_un(5).*y(238,:) + kcat10_un(6).*y(253,:) + kcat10_un(7).*y(268,:) + kcat10_un(8).*y(283,:)...
    + kcat8_CO2.*y(327,:) + kcat8_H.*y(331,:) + kcat10_CO2.*y(332,:) + kcat10_H.*y(336,:)...
    + kcat8_un(4).*y(308,:)...
;

dy(12,:) = kcat3.*y(92,:) + k4_2r(1).*y(97,:) - k4_2f(1).*y(93,:).*y(12,:)...
    + kcat8_H.*y(331,:) + kcat10_H.*y(336,:);

dy(17,:) = kcat8(1).*y(104,:)+ kcat10(1).*y(109,:) + k4_2r(2).*y(112,:)- k4_2f(2).*y(93,:).*y(17,:);
dy(22,:) = kcat8(2).*y(119,:)+ kcat10(2).*y(124,:) + k4_2r(3).*y(127,:)- k4_2f(3).*y(93,:).*y(22,:);
dy(27,:) = kcat8(3).*y(134,:)+ kcat10(3).*y(139,:) + k4_2r(4).*y(142,:)- k4_2f(4).*y(93,:).*y(27,:);
dy(33,:) = kcat8(4).*y(149,:)+ kcat10(4).*y(154,:) + k4_2r(5).*y(157,:)- k4_2f(5).*y(93,:).*y(33,:);
dy(43,:) = kcat8(5).*y(164,:)+ kcat10(5).*y(169,:) + k4_2r(6).*y(172,:)- k4_2f(6).*y(93,:).*y(43,:);
dy(53,:) = kcat8(6).*y(179,:)+ kcat10(6).*y(184,:) + k4_2r(7).*y(187,:)- k4_2f(7).*y(93,:).*y(53,:);
dy(63,:) = kcat8(7).*y(194,:)+ kcat10(7).*y(199,:) + k4_2r(8).*y(202,:)- k4_2f(8).*y(93,:).*y(63,:);
dy(73,:) = kcat8(8).*y(209,:)+ kcat10(8).*y(214,:) + k4_2r(9).*y(217,:)- k4_2f(9).*y(93,:).*y(73,:);



dy(38,:) = kcat8_un(4).*y(308,:) + kcat10_un(4).*y(304,:) + k4_2r(5).*y(226,:)- k4_2f(5).*y(93,:).*y(38,:);
dy(48,:) = kcat8_un(5).*y(233,:) + kcat10_un(5).*y(238,:) + k4_2r(6).*y(241,:)- k4_2f(6).*y(93,:).*y(48,:);
dy(58,:) = kcat8_un(6).*y(248,:) + kcat10_un(6).*y(253,:) + k4_2r(7).*y(256,:)- k4_2f(7).*y(93,:).*y(58,:);
dy(68,:) = kcat8_un(7).*y(263,:) + kcat10_un(7).*y(268,:) + k4_2r(8).*y(271,:)- k4_2f(8).*y(93,:).*y(68,:);
dy(78,:) = kcat8_un(8).*y(278,:) + kcat10_un(8).*y(283,:) + k4_2r(9).*y(286,:)- k4_2f(9).*y(93,:).*y(78,:);

dy(13,:) =  kcat4(1).*y(97,:) + k5_1r(1).*y(98,:) - k5_1f(1).*e5.*y(13,:) + k9_1r(1).*y(105,:) - k9_1f(1).*e9.*y(13,:);
dy(18,:) =  kcat4(2).*y(112,:) + k5_1r(2).*y(113,:) - k5_1f(2).*e5.*y(18,:) + k9_1r(2).*y(120,:) - k9_1f(2).*e9.*y(18,:);
dy(23,:) =  kcat4(3).*y(127,:) + k5_1r(3).*y(128,:) - k5_1f(3).*e5.*y(23,:) + k9_1r(3).*y(135,:) - k9_1f(3).*e9.*y(23,:);
dy(28,:) = kcat4(4).*y(142,:)+ k5_1r(4).*y(143,:)- k5_1f(4).*e5.*y(28,:)+ k9_1r(4).*y(150,:)- k9_1f(4).*e9.*y(28,:);
dy(34,:) = kcat4(5).*y(157,:)+ k5_1r(5).*y(158,:)- k5_1f(5).*e5.*y(34,:)+ k9_1r(5).*y(165,:)- k9_1f(5).*e9.*y(34,:);
dy(44,:) = kcat4(6).*y(172,:)+ k5_1r(6).*y(173,:)- k5_1f(6).*e5.*y(44,:)+ k9_1r(6).*y(180,:)- k9_1f(6).*e9.*y(44,:)...
            +LpxA.*A_kcat1.*y(314,:).*y(4,:)./(Kma1.*Kmb1+y(314,:).*Kmb1+y(4,:).*Kma1+y(314,:).*y(4,:))...
            -LpxD.*A_kcat10./(1+y(311,:)./ki2).*y(318,:).*y(44,:)./(Kma4.*Kmb4+y(318,:).*Kmb4+y(44,:).*Kma4+y(318,:).*y(44,:))...
            -LpxA.*A_kcat.*UDP_Glc.*y(44,:)./(Kma.*Kmb+UDP_Glc.*Kmb+y(44,:).*Kma+UDP_Glc.*y(44,:));
dy(54,:) = kcat4(7).*y(187,:)+ k5_1r(7).*y(188,:)- k5_1f(7).*e5.*y(54,:)+ k9_1r(7).*y(195,:)- k9_1f(7).*e9.*y(54,:);
dy(64,:) = kcat4(8).*y(202,:)+ k5_1r(8).*y(203,:)- k5_1f(8).*e5.*y(64,:)+ k9_1r(8).*y(210,:)- k9_1f(8).*e9.*y(64,:);
dy(74,:) = kcat4(9).*y(217,:)+ k5_1r(9).*y(218,:)- k5_1f(9).*e5.*y(74,:)+ k9_1r(9).*y(222,:)- k9_1f(9).*e9.*y(74,:);

dy(39,:) = kcat4(5).*y(226,:)+ k5_1r(5).*y(227,:)- k5_1f(5).*e5.*y(39,:)+ k9_1r_un(5).*y(234,:)- k9_1f_un(5).*e9.*y(39,:);
dy(49,:) = kcat4(6).*y(241,:)+ k5_1r(6).*y(242,:)- k5_1f(6).*e5.*y(49,:)+ k9_1r_un(6).*y(249,:)- k9_1f_un(6).*e9.*y(49,:);
dy(59,:) = kcat4(7).*y(256,:)+ k5_1r(7).*y(257,:)- k5_1f(7).*e5.*y(59,:)+ k9_1r_un(7).*y(264,:)- k9_1f_un(7).*e9.*y(59,:);
dy(69,:) = kcat4(8).*y(271,:)+ k5_1r(8).*y(272,:)- k5_1f(8).*e5.*y(69,:)+ k9_1r_un(8).*y(279,:)- k9_1f_un(8).*e9.*y(69,:);
dy(79,:) = kcat4(9).*y(286,:)+ k5_1r(9).*y(287,:)- k5_1f(9).*e5.*y(79,:)+ k9_1r_un(9).*y(291,:)- k9_1f_un(9).*e9.*y(79,:);



dy(14,:) =  k5_3f(1).*y(99,:) - k5_3r(1).*e5.*y(14,:) +k9_3f(1).*y(106,:) - k9_3r(1).*e9.*y(14,:) + k6_2r(1).*y(100,:) - k6_2f(1).*y(94,:).*y(14,:);
dy(19,:) =  k5_3f(2).*y(114,:) - k5_3r(2).*e5.*y(19,:) +k9_3f(2).*y(121,:) - k9_3r(2).*e9.*y(19,:) + k6_2r(2).*y(115,:) - k6_2f(2).*y(94,:).*y(19,:);
dy(24,:) =  k5_3f(3).*y(129,:) - k5_3r(3).*e5.*y(24,:) +k9_3f(3).*y(136,:) - k9_3r(3).*e9.*y(24,:) + k6_2r(3).*y(130,:) - k6_2f(3).*y(94,:).*y(24,:);
dy(29,:) = k5_3f(4).*y(144,:) -k5_3r(4).*e5.*y(29,:)+k9_3f(4).*y(151,:) -k9_3r(4).*e9.*y(29,:)+ k6_2r(4).*y(145,:)- k6_2f(4).*y(94,:).*y(29,:);
dy(35,:) = k5_3f(5).*y(159,:) -k5_3r(5).*e5.*y(35,:)+k9_3f(5).*y(166,:) -k9_3r(5).*e9.*y(35,:)+ k6_2r(5).*y(160,:)- k6_2f(5).*y(94,:).*y(35,:);
dy(45,:) = k5_3f(6).*y(174,:) -k5_3r(6).*e5.*y(45,:)+k9_3f(6).*y(181,:) -k9_3r(6).*e9.*y(45,:)+ k6_2r(6).*y(175,:)- k6_2f(6).*y(94,:).*y(45,:);
dy(55,:) = k5_3f(7).*y(189,:) -k5_3r(7).*e5.*y(55,:)+k9_3f(7).*y(196,:) -k9_3r(7).*e9.*y(55,:)+ k6_2r(7).*y(190,:)- k6_2f(7).*y(94,:).*y(55,:);
dy(65,:) = k5_3f(8).*y(204,:) -k5_3r(8).*e5.*y(65,:)+k9_3f(8).*y(211,:) -k9_3r(8).*e9.*y(65,:)+ k6_2r(8).*y(205,:)- k6_2f(8).*y(94,:).*y(65,:);
dy(75,:) = k5_3f(9).*y(219,:) -k5_3r(9).*e5.*y(75,:)+k9_3f(9).*y(223,:) -k9_3r(9).*e9.*y(75,:)+ k6_2r(9).*y(220,:)- k6_2f(9).*y(94,:).*y(75,:);

dy(32,:) = k9_3f_un(4).*y(301,:) - k9_3r_un(4).*e9.*y(32,:) + k10_1r(4).*y(302,:) - k10_1f(4).*e10.*y(32,:)...
          +k5_3f_un(4).*y(305,:) - k5_3r_un(4).*e5.*y(32,:) + k8_1r(4).*y(306,:) - k8_1f(4).*e8.*y(32,:);
dy(40,:) = k5_3f(5).*y(228,:) -k5_3r(5).*e5.*y(40,:)+k9_3f(5).*y(235,:) -k9_3r(5).*e9.*y(40,:)+ k6_2r(5).*y(229,:)- k6_2f(5).*y(94,:).*y(40,:);
dy(50,:) = k5_3f(6).*y(243,:) -k5_3r(6).*e5.*y(50,:)+k9_3f(6).*y(250,:) -k9_3r(6).*e9.*y(50,:)+ k6_2r(6).*y(244,:)- k6_2f(6).*y(94,:).*y(50,:);
dy(60,:) = k5_3f(7).*y(258,:) -k5_3r(7).*e5.*y(60,:)+k9_3f(7).*y(265,:) -k9_3r(7).*e9.*y(60,:)+ k6_2r(7).*y(259,:)- k6_2f(7).*y(94,:).*y(60,:);
dy(70,:) = k5_3f(8).*y(273,:) -k5_3r(8).*e5.*y(70,:)+k9_3f(8).*y(280,:) -k9_3r(8).*e9.*y(70,:)+ k6_2r(8).*y(274,:)- k6_2f(8).*y(94,:).*y(70,:);
dy(80,:) = k5_3f(9).*y(288,:) -k5_3r(9).*e5.*y(80,:)+k9_3f(9).*y(292,:) -k9_3r(9).*e9.*y(80,:)+ k6_2r(9).*y(289,:)- k6_2f(9).*y(94,:).*y(80,:);

dy(15,:) =   kcat6(1).*y(100,:) + k7_1r(1).*y(101,:) + k8_1r(1).*y(102,:)+  k10_1r(1).*y(107,:)  - k7_1f(1).*e7.*y(15,:) - k8_1f(1).*e8.*y(15,:) - k10_1f(1).*e10.*y(15,:) - k3_4f(1).*e3.*y(15,:) + k3_4r(1).*y(95,:) - k3_5f(1).*y(91,:).*y(15,:) + k3_5r(1).*y(96,:);
dy(20,:) =   kcat6(2).*y(115,:) + k7_1r(2).*y(116,:) + k8_1r(2).*y(117,:)+  k10_1r(2).*y(122,:)  - k7_1f(2).*e7.*y(20,:) - k8_1f(2).*e8.*y(20,:) - k10_1f(2).*e10.*y(20,:) - k3_4f(2).*e3.*y(20,:) + k3_4r(2).*y(110,:) - k3_5f(2).*y(91,:).*y(20,:) + k3_5r(2).*y(111,:);
dy(25,:) =   kcat6(3).*y(130,:) + k7_1r(3).*y(131,:) + k8_1r(3).*y(132,:)+  k10_1r(3).*y(137,:)  - k7_1f(3).*e7.*y(25,:) - k8_1f(3).*e8.*y(25,:) - k10_1f(3).*e10.*y(25,:) - k3_4f(3).*e3.*y(25,:) + k3_4r(3).*y(125,:) - k3_5f(3).*y(91,:).*y(25,:) + k3_5r(3).*y(126,:);
dy(30,:) =  kcat6(4).*y(145,:)+ k7_1r(4).*y(146,:)+ k8_1r(4).*y(147,:)+ k10_1r(4).*y(152,:) - k7_1f(4).*e7.*y(30,:)- k8_1f(4).*e8.*y(30,:)- k10_1f(4).*e10.*y(30,:)- k3_4f(4).*e3.*y(30,:)+ k3_4r(4).*y(140,:)- k3_5f(4).*y(91,:).*y(30,:)+ k3_5r(4).*y(141,:);
dy(36,:) =  kcat6(5).*y(160,:)+ k7_1r(5).*y(161,:)+ k8_1r(5).*y(162,:)+ k10_1r(5).*y(167,:) - k7_1f(5).*e7.*y(36,:)- k8_1f(5).*e8.*y(36,:)- k10_1f(5).*e10.*y(36,:)- k3_4f(5).*e3.*y(36,:)+ k3_4r(5).*y(155,:)- k3_5f(5).*y(91,:).*y(36,:)+ k3_5r(5).*y(156,:)...
            -A_kcat7.*LpxL.*y(320,:).*y(36,:)./((Km3+y(320,:).*(km_laur+y(36,:))));
dy(46,:) =  kcat6(6).*y(175,:)+ k7_1r(6).*y(176,:)+ k8_1r(6).*y(177,:)+ k10_1r(6).*y(182,:) - k7_1f(6).*e7.*y(46,:)- k8_1f(6).*e8.*y(46,:)- k10_1f(6).*e10.*y(46,:)- k3_4f(6).*e3.*y(46,:)+ k3_4r(6).*y(170,:)- k3_5f(6).*y(91,:).*y(46,:)+ k3_5r(6).*y(171,:)...
            -A_kcat8.*LpxM.*y(317,:).*y(46,:)./((Km4+y(317,:).*(km_myr+y(46,:))));
dy(56,:) =  kcat6(7).*y(190,:)+ k7_1r(7).*y(191,:)+ k8_1r(7).*y(192,:)+ k10_1r(7).*y(197,:) - k7_1f(7).*e7.*y(56,:)- k8_1f(7).*e8.*y(56,:)- k10_1f(7).*e10.*y(56,:)- k3_4f(7).*e3.*y(56,:)+ k3_4r(7).*y(185,:)- k3_5f(7).*y(91,:).*y(56,:)+ k3_5r(7).*y(186,:)...
            -kcatP1.*PlsB.*G3P.*y(56,:)./((KmP1+G3P).*(KmP2+y(56,:)));
dy(66,:) =  kcat6(8).*y(205,:)+ k7_1r(8).*y(206,:)+ k8_1r(8).*y(207,:)+ k10_1r(8).*y(212,:) - k7_1f(8).*e7.*y(66,:)- k8_1f(8).*e8.*y(66,:)- k10_1f(8).*e10.*y(66,:)- k3_4f(8).*e3.*y(66,:)+ k3_4r(8).*y(200,:)- k3_5f(8).*y(91,:).*y(66,:)+ k3_5r(8).*y(201,:);
dy(76,:) =  kcat6(9).*y(220,:)+ k7_1r(9).*y(221,:)                  - k7_1f(9).*e7.*y(76,:)                  - k3_4f(9).*e3.*y(76,:)+ k3_4r(9).*y(215,:)- k3_5f(9).*y(91,:).*y(76,:)+ k3_5r(9).*y(216,:);


dy(41,:) =  kcat6(5).*y(229,:)+ k7_1r(5).*y(230,:)+ k8_1r(5).*y(231,:)+ k10_1r(5).*y(236,:) - k7_1f(5).*e7.*y(41,:)- k8_1f(5).*e8.*y(41,:)- k10_1f(5).*e10.*y(41,:)- k3_4f(5).*e3.*y(41,:)+ k3_4r(5).*y(224,:)- k3_5f(5).*y(91,:).*y(41,:)+ k3_5r(5).*y(225,:);
dy(51,:) =  kcat6(6).*y(244,:)+ k7_1r(6).*y(245,:)+ k8_1r(6).*y(246,:)+ k10_1r(6).*y(251,:) - k7_1f(6).*e7.*y(51,:)- k8_1f(6).*e8.*y(51,:)- k10_1f(6).*e10.*y(51,:)- k3_4f(6).*e3.*y(51,:)+ k3_4r(6).*y(239,:)- k3_5f(6).*y(91,:).*y(51,:)+ k3_5r(6).*y(240,:);
dy(61,:) =  kcat6(7).*y(259,:)+ k7_1r(7).*y(260,:)+ k8_1r(7).*y(261,:)+ k10_1r(7).*y(266,:) - k7_1f(7).*e7.*y(61,:)- k8_1f(7).*e8.*y(61,:)- k10_1f(7).*e10.*y(61,:)- k3_4f(7).*e3.*y(61,:)+ k3_4r(7).*y(254,:)- k3_5f(7).*y(91,:).*y(61,:)+ k3_5r(7).*y(255,:)...
            -kcatP2*PlsC.*y(324,:).*y(61,:)./((KmP4 + y(324,:)).*(KmP5 + y(61,:)));
dy(71,:) =  kcat6(8).*y(274,:)+ k7_1r(8).*y(275,:)+ k8_1r(8).*y(276,:)+ k10_1r(8).*y(281,:) - k7_1f(8).*e7.*y(71,:)- k8_1f(8).*e8.*y(71,:)- k10_1f(8).*e10.*y(71,:)- k3_4f(8).*e3.*y(71,:)+ k3_4r(8).*y(269,:)- k3_5f(8).*y(91,:).*y(71,:)+ k3_5r(8).*y(270,:)...
            -kcatP1.*PlsB.*G3P.*y(71,:)./((KmP1+G3P).*(KmP3+y(71,:))) - kcatP2*PlsC.*y(324,:).*y(71,:)./((KmP4 + y(324,:)).*(KmP6 + y(71,:)));
dy(81,:) =  kcat6(9).*y(289,:)+ k7_1r(9).*y(290,:)                  - k7_1f(9).*e7.*y(81,:)                  - k3_4f(9).*e3.*y(81,:)+ k3_4r(9).*y(284,:)- k3_5f(9).*y(91,:).*y(81,:)+ k3_5r(9).*y(285,:);


%FFA Production. 
dy(16,:) = kcat7(1).*y(101,:);
dy(21,:) = kcat7(2).*y(116,:) - kcatD1(1).*FadD.*y(1,:).*y(21,:)./((Km_ATP+y(1,:)).*(Km_FA(1)+y(21,:)));
dy(26,:) = kcat7(3).*y(131,:) - kcatD1(2).*FadD.*y(1,:).*y(26,:)./((Km_ATP+y(1,:)).*(Km_FA(2)+y(26,:)));
dy(31,:) = kcat7(4).*y(146,:) - kcatD1(3).*FadD.*y(1,:).*y(31,:)./((Km_ATP+y(1,:)).*(Km_FA(3)+y(31,:)));
dy(37,:) = kcat7(5).*y(161,:) - kcatD1(4).*FadD.*y(1,:).*y(37,:)./((Km_ATP+y(1,:)).*(Km_FA(4)+y(37,:)));
dy(47,:) = kcat7(6).*y(176,:) - kcatD1(5).*FadD.*y(1,:).*y(47,:)./((Km_ATP+y(1,:)).*(Km_FA(5)+y(47,:)));
dy(57,:) = kcat7(7).*y(191,:) - kcatD1(6).*FadD.*y(1,:).*y(57,:)./((Km_ATP+y(1,:)).*(Km_FA(6)+y(57,:)));
dy(67,:) = kcat7(8).*y(206,:) - kcatD1(7).*FadD.*y(1,:).*y(67,:)./((Km_ATP+y(1,:)).*(Km_FA(7)+y(67,:)));
dy(77,:) = kcat7(9).*y(221,:);

dy(42,:) = kcat7(5).*y(230,:) - kcatD1(8).*FadD.*y(1,:).*y(42,:)./((Km_ATP+y(1,:)).*(Km_FA(8)+y(42,:)));
dy(52,:) = kcat7(6).*y(245,:) - kcatD1(9).*FadD.*y(1,:).*y(52,:)./((Km_ATP+y(1,:)).*(Km_FA(9)+y(52,:)));
dy(62,:) = kcat7(7).*y(260,:) - kcatD1(10).*FadD.*y(1,:).*y(62,:)./((Km_ATP+y(1,:)).*(Km_FA(10)+y(62,:)));
dy(72,:) = kcat7(8).*y(275,:) - kcatD1(11).*FadD.*y(1,:).*y(72,:)./((Km_ATP+y(1,:)).*(Km_FA(11)+y(72,:)));
dy(82,:) = kcat7(9).*y(290,:);

dy(83,:) = k1_1f.*e1.*y(1,:) + k1_2r.*y(84,:) - k1_1r.*y(83,:) - k1_2f.*y(83,:).*y(2,:);
dy(84,:) = k1_2f.*y(83,:).*y(2,:) - k1_2r.*y(84,:) - kcat1_1.*y(84,:);
dy(85,:) = kcat1_1.*y(84,:) + k1_3r.*y(86,:) - k1_3f.*y(85,:).*y(3,:);
dy(86,:) = k1_3f.*y(85,:).*y(3,:) - k1_3r.*y(86,:) - kcat1_2.*y(86,:);

dy(87,:) = k2_1f.*e2.*y(8,:) + k2_2r.*y(88,:).*y(9,:) - k2_1r.*y(87,:) - k2_2f.*y(87,:);
dy(88,:) = k2_2f.*y(87,:) + k2_3r.*y(89,:) - k2_2r.*y(88,:).*y(9,:) - k2_3f.*y(88,:).*y(4,:);
dy(89,:) = k2_3f.*y(88,:).*y(4,:) + k2_4r.*e2.*y(10,:) - k2_3r.*y(89,:) - k2_4f.*y(89,:);

dy(90,:) = k3_1f.*e3.*y(3,:) + k3_2r.*y(91,:).*y(9,:) - k3_1r.*y(90,:) - k3_2f.*y(90,:);
dy(91,:) = k3_2f.*y(90,:) + k3_3r.*y(92,:) - k3_2r.*y(91,:).*y(9,:) - k3_3f.*y(91,:).*y(10,:)...
    + k3_5r(1).*y(96,:)  - k3_5f(1).*y(91,:).*y(15,:)...
    + k3_5r(2).*y(111,:)  - k3_5f(2).*y(91,:).*y(20,:)...
    + k3_5r(3).*y(126,:)  - k3_5f(3).*y(91,:).*y(25,:)...
    + k3_5r(4).*y(141,:) - k3_5f(4).*y(91,:).*y(30,:)...
    + k3_5r(5).*y(156,:) - k3_5f(5).*y(91,:).*y(36,:)...
    + k3_5r(6).*y(171,:) - k3_5f(6).*y(91,:).*y(46,:)...
    + k3_5r(7).*y(186,:) - k3_5f(7).*y(91,:).*y(56,:)...
    + k3_5r(8).*y(201,:) - k3_5f(8).*y(91,:).*y(66,:)...
    + k3_5r(9).*y(216,:) - k3_5f(9).*y(91,:).*y(76,:)...
    + k3_5r(5).*y(225,:) - k3_5f(5).*y(91,:).*y(41,:)...
    + k3_5r(6).*y(240,:) - k3_5f(6).*y(91,:).*y(51,:)...
    + k3_5r(7).*y(255,:) - k3_5f(7).*y(91,:).*y(61,:)...
    + k3_5r(8).*y(270,:) - k3_5f(8).*y(91,:).*y(71,:)...
    + k3_5r(9).*y(285,:) - k3_5f(9).*y(91,:).*y(81,:)...
    ;
dy(92,:) = k3_3f.*y(91,:).*y(10,:) - k3_3r.*y(92,:) - kcat3.*y(92,:);
dy(93,:) = k4_1f(1).*e4.*y(5,:) - k4_1r(1).*y(93,:)...
    + k4_2r(1).*y(97,:) - k4_2f(1).*y(93,:).*y(12,:)...
    + k4_2r(2).*y(112,:) - k4_2f(2).*y(93,:).*y(17,:)...
    + k4_2r(3).*y(127,:) - k4_2f(3).*y(93,:).*y(22,:)...
    + k4_2r(4).*y(142,:)- k4_2f(4).*y(93,:).*y(27,:)...
    + k4_2r(5).*y(157,:)- k4_2f(5).*y(93,:).*y(33,:)...
    + k4_2r(6).*y(172,:)- k4_2f(6).*y(93,:).*y(43,:)...
    + k4_2r(7).*y(187,:)- k4_2f(7).*y(93,:).*y(53,:)...
    + k4_2r(8).*y(202,:)- k4_2f(8).*y(93,:).*y(63,:)...
    + k4_2r(9).*y(217,:)- k4_2f(9).*y(93,:).*y(73,:)...
    + k4_2r(5).*y(226,:)- k4_2f(5).*y(93,:).*y(38,:)...
    + k4_2r(6).*y(241,:)- k4_2f(6).*y(93,:).*y(48,:)...
    + k4_2r(7).*y(256,:)- k4_2f(7).*y(93,:).*y(58,:)...
    + k4_2r(8).*y(271,:)- k4_2f(8).*y(93,:).*y(68,:)...
    + k4_2r(9).*y(286,:)- k4_2f(9).*y(93,:).*y(78,:)...
;
dy(97,:) =  k4_2f(1).*y(93,:).*y(12,:) - k4_2r(1).*y(97,:) - kcat4(1).*y(97,:);
dy(112,:) =  k4_2f(2).*y(93,:).*y(17,:) - k4_2r(2).*y(112,:) - kcat4(2).*y(112,:);
dy(127,:) =  k4_2f(3).*y(93,:).*y(22,:) - k4_2r(3).*y(127,:) - kcat4(3).*y(127,:);
dy(142,:) = k4_2f(4).*y(93,:).*y(27,:)- k4_2r(4).*y(142,:)- kcat4(4).*y(142,:);
dy(157,:) = k4_2f(5).*y(93,:).*y(33,:)- k4_2r(5).*y(157,:)- kcat4(5).*y(157,:);
dy(172,:) = k4_2f(6).*y(93,:).*y(43,:)- k4_2r(6).*y(172,:)- kcat4(6).*y(172,:);
dy(187,:) = k4_2f(7).*y(93,:).*y(53,:)- k4_2r(7).*y(187,:)- kcat4(7).*y(187,:);
dy(202,:) = k4_2f(8).*y(93,:).*y(63,:)- k4_2r(8).*y(202,:)- kcat4(8).*y(202,:);
dy(217,:) = k4_2f(9).*y(93,:).*y(73,:)- k4_2r(9).*y(217,:)- kcat4(9).*y(217,:);

dy(226,:) = k4_2f(5).*y(93,:).*y(38,:)- k4_2r(5).*y(226,:)- kcat4(5).*y(226,:);
dy(241,:) = k4_2f(6).*y(93,:).*y(48,:)- k4_2r(6).*y(241,:)- kcat4(6).*y(241,:);
dy(256,:) = k4_2f(7).*y(93,:).*y(58,:)- k4_2r(7).*y(256,:)- kcat4(7).*y(256,:);
dy(271,:) = k4_2f(8).*y(93,:).*y(68,:)- k4_2r(8).*y(271,:)- kcat4(8).*y(271,:);
dy(286,:) = k4_2f(9).*y(93,:).*y(78,:)- k4_2r(9).*y(286,:)- kcat4(9).*y(286,:);

dy(98,:) =  k5_1f(1).*e5.*y(13,:) - k5_1r(1).*y(98,:) - kcat5(1).*y(98,:) + k5_2r(1).*y(99,:);
dy(113,:) =  k5_1f(2).*e5.*y(18,:) - k5_1r(2).*y(113,:) - kcat5(2).*y(113,:) + k5_2r(2).*y(114,:);
dy(128,:) =  k5_1f(3).*e5.*y(23,:) - k5_1r(3).*y(128,:) - kcat5(3).*y(128,:) + k5_2r(3).*y(129,:);
dy(143,:) = k5_1f(4).*e5.*y(28,:)- k5_1r(4).*y(143,:)- kcat5(4).*y(143,:) +k5_2r(4).*y(144,:);
dy(158,:) = k5_1f(5).*e5.*y(34,:)- k5_1r(5).*y(158,:)- kcat5(5).*y(158,:) +k5_2r(5).*y(159,:);
dy(173,:) = k5_1f(6).*e5.*y(44,:)- k5_1r(6).*y(173,:)- kcat5(6).*y(173,:) +k5_2r(6).*y(174,:);
dy(188,:) = k5_1f(7).*e5.*y(54,:)- k5_1r(7).*y(188,:)- kcat5(7).*y(188,:) +k5_2r(7).*y(189,:);
dy(203,:) = k5_1f(8).*e5.*y(64,:)- k5_1r(8).*y(203,:)- kcat5(8).*y(203,:) +k5_2r(8).*y(204,:);
dy(218,:) = k5_1f(9).*e5.*y(74,:)- k5_1r(9).*y(218,:)- kcat5(9).*y(218,:) +k5_2r(9).*y(219,:);

dy(227,:) = k5_1f(5).*e5.*y(39,:)- k5_1r(5).*y(227,:)- kcat5(5).*y(227,:) + k5_2r(5).*y(228,:);
dy(242,:) = k5_1f(6).*e5.*y(49,:)- k5_1r(6).*y(242,:)- kcat5(6).*y(242,:) + k5_2r(6).*y(243,:);
dy(257,:) = k5_1f(7).*e5.*y(59,:)- k5_1r(7).*y(257,:)- kcat5(7).*y(257,:) + k5_2r(7).*y(258,:);
dy(272,:) = k5_1f(8).*e5.*y(69,:)- k5_1r(8).*y(272,:)- kcat5(8).*y(272,:) + k5_2r(8).*y(273,:);
dy(287,:) = k5_1f(9).*e5.*y(79,:)- k5_1r(9).*y(287,:)- kcat5(9).*y(287,:) + k5_2r(9).*y(288,:);


dy(99,:) =  kcat5(1).*y(98,:) - k5_2r(1).*y(99,:) - k5_3f(1).*y(99,:) + k5_3r(1).*e5.*y(14,:);
dy(114,:) =  kcat5(2).*y(113,:) - k5_2r(2).*y(114,:) - k5_3f(2).*y(114,:) + k5_3r(2).*e5.*y(19,:);
dy(129,:) =  kcat5(3).*y(128,:) - k5_2r(3).*y(129,:) - k5_3f(3).*y(129,:) + k5_3r(3).*e5.*y(24,:);
dy(144,:) = kcat5(4).*y(143,:)- k5_2r(4).*y(144,:)- k5_3f(4).*y(144,:)+ k5_3r(4).*e5.*y(29,:)- kcat5_un(4).*y(144,:)+ k5_2r_un(4).*y(305,:);
dy(159,:) = kcat5(5).*y(158,:)- k5_2r(5).*y(159,:)- k5_3f(5).*y(159,:)+ k5_3r(5).*e5.*y(35,:);
dy(174,:) = kcat5(6).*y(173,:)- k5_2r(6).*y(174,:)- k5_3f(6).*y(174,:)+ k5_3r(6).*e5.*y(45,:);
dy(189,:) = kcat5(7).*y(188,:)- k5_2r(7).*y(189,:)- k5_3f(7).*y(189,:)+ k5_3r(7).*e5.*y(55,:);
dy(204,:) = kcat5(8).*y(203,:)- k5_2r(8).*y(204,:)- k5_3f(8).*y(204,:)+ k5_3r(8).*e5.*y(65,:);
dy(219,:) = kcat5(9).*y(218,:)- k5_2r(9).*y(219,:)- k5_3f(9).*y(219,:)+ k5_3r(9).*e5.*y(75,:);

dy(228,:) = kcat5(5).*y(227,:) - k5_2r(5).*y(228,:) - k5_3f(5).*y(228,:) + k5_3r(5).*e5.*y(40,:);
dy(243,:) = kcat5(6).*y(242,:) - k5_2r(6).*y(243,:) - k5_3f(6).*y(243,:) + k5_3r(6).*e5.*y(50,:);
dy(258,:) = kcat5(7).*y(257,:) - k5_2r(7).*y(258,:) - k5_3f(7).*y(258,:) + k5_3r(7).*e5.*y(60,:);
dy(273,:) = kcat5(8).*y(272,:) - k5_2r(8).*y(273,:) - k5_3f(8).*y(273,:) + k5_3r(8).*e5.*y(70,:);
dy(288,:) = kcat5(9).*y(287,:) - k5_2r(9).*y(288,:) - k5_3f(9).*y(288,:) + k5_3r(9).*e5.*y(80,:);


dy(105,:) =  k9_1f(1).*e9.*y(13,:) - k9_1r(1).*y(105,:) - kcat9(1).*y(105,:) + k9_2r(1).*y(106,:);
dy(120,:) =  k9_1f(2).*e9.*y(18,:) - k9_1r(2).*y(120,:) - kcat9(2).*y(120,:) + k9_2r(2).*y(121,:);
dy(135,:) =  k9_1f(3).*e9.*y(23,:) - k9_1r(3).*y(135,:) - kcat9(3).*y(135,:) + k9_2r(3).*y(136,:);
dy(150,:) = k9_1f(4).*e9.*y(28,:)- k9_1r(4).*y(150,:)- kcat9(4).*y(150,:) +k9_2r(4).*y(151,:);
dy(165,:) = k9_1f(5).*e9.*y(34,:)- k9_1r(5).*y(165,:)- kcat9(5).*y(165,:) +k9_2r(5).*y(166,:);
dy(180,:) = k9_1f(6).*e9.*y(44,:)- k9_1r(6).*y(180,:)- kcat9(6).*y(180,:) +k9_2r(6).*y(181,:);
dy(195,:) = k9_1f(7).*e9.*y(54,:)- k9_1r(7).*y(195,:)- kcat9(7).*y(195,:) +k9_2r(7).*y(196,:);
dy(210,:) = k9_1f(8).*e9.*y(64,:)- k9_1r(8).*y(210,:)- kcat9(8).*y(210,:) +k9_2r(8).*y(211,:);
dy(222,:) = k9_1f(9).*e9.*y(74,:)- k9_1r(9).*y(222,:)- kcat9(9).*y(222,:) +k9_2r(9).*y(223,:);

dy(234,:) = k9_1f_un(5).*e9.*y(39,:)- k9_1r_un(5).*y(234,:)- kcat9(5).*y(234,:) + k9_2r(5).*y(235,:);
dy(249,:) = k9_1f_un(6).*e9.*y(49,:)- k9_1r_un(6).*y(249,:)- kcat9(6).*y(249,:) + k9_2r(6).*y(250,:);
dy(264,:) = k9_1f_un(7).*e9.*y(59,:)- k9_1r_un(7).*y(264,:)- kcat9(7).*y(264,:) + k9_2r(7).*y(265,:);
dy(279,:) = k9_1f_un(8).*e9.*y(69,:)- k9_1r_un(8).*y(279,:)- kcat9(8).*y(279,:) + k9_2r(8).*y(280,:);
dy(291,:) = k9_1f_un(9).*e9.*y(79,:)- k9_1r_un(9).*y(291,:)- kcat9(9).*y(291,:) + k9_2r(9).*y(292,:);

dy(106,:) =  kcat9(1).*y(105,:) - k9_2r(1).*y(106,:) - k9_3f(1).*y(106,:) + k9_3r(1).*e9.*y(14,:);
dy(121,:) =  kcat9(2).*y(120,:) - k9_2r(2).*y(121,:) - k9_3f(2).*y(121,:) + k9_3r(2).*e9.*y(19,:);
dy(136,:) =  kcat9(3).*y(135,:) - k9_2r(3).*y(136,:) - k9_3f(3).*y(136,:) + k9_3r(3).*e9.*y(24,:);
dy(151,:) = kcat9(4).*y(150,:)- k9_2r(4).*y(151,:)- k9_3f(4).*y(151,:)+ k9_3r(4).*e9.*y(29,:) - kcat9_un(4).*y(151,:) + k9_2r_un(4).*y(301,:);
dy(166,:) = kcat9(5).*y(165,:)- k9_2r(5).*y(166,:)- k9_3f(5).*y(166,:)+ k9_3r(5).*e9.*y(35,:);
dy(181,:) = kcat9(6).*y(180,:)- k9_2r(6).*y(181,:)- k9_3f(6).*y(181,:)+ k9_3r(6).*e9.*y(45,:);
dy(196,:) = kcat9(7).*y(195,:)- k9_2r(7).*y(196,:)- k9_3f(7).*y(196,:)+ k9_3r(7).*e9.*y(55,:);
dy(211,:) = kcat9(8).*y(210,:)- k9_2r(8).*y(211,:)- k9_3f(8).*y(211,:)+ k9_3r(8).*e9.*y(65,:);
dy(223,:) = kcat9(9).*y(222,:)- k9_2r(9).*y(223,:)- k9_3f(9).*y(223,:)+ k9_3r(9).*e9.*y(75,:);

dy(301,:) = kcat9_un(4).*y(151,:) - k9_2r_un(4).*y(301,:) - k9_3f_un(4).*y(301,:) + k9_3r_un(4).*e9.*y(32,:);

dy(235,:) = kcat9(5).*y(234,:) - k9_2r(5).*y(235,:) - k9_3f(5).*y(235,:) + k9_3r(5).*e9.*y(40,:);
dy(250,:) = kcat9(6).*y(249,:) - k9_2r(6).*y(250,:) - k9_3f(6).*y(250,:) + k9_3r(6).*e9.*y(50,:);
dy(265,:) = kcat9(7).*y(264,:) - k9_2r(7).*y(265,:) - k9_3f(7).*y(265,:) + k9_3r(7).*e9.*y(60,:);
dy(280,:) = kcat9(8).*y(279,:) - k9_2r(8).*y(280,:) - k9_3f(8).*y(280,:) + k9_3r(8).*e9.*y(70,:);
dy(292,:) = kcat9(9).*y(291,:) - k9_2r(9).*y(292,:) - k9_3f(9).*y(292,:) + k9_3r(9).*e9.*y(80,:);


dy(94,:) = k6_1f(1).*e6.*y(6,:) - k6_1r(1).*y(94,:)...
    + k6_2r(1).*y(100,:) - k6_2f(1).*y(94,:).*y(14,:)...
    + k6_2r(2).*y(115,:) - k6_2f(2).*y(94,:).*y(19,:)...
    + k6_2r(3).*y(130,:) - k6_2f(3).*y(94,:).*y(24,:)...
    + k6_2r(4).*y(145,:)- k6_2f(4).*y(94,:).*y(29,:)...
    + k6_2r(5).*y(160,:)- k6_2f(5).*y(94,:).*y(35,:)...
    + k6_2r(6).*y(175,:)- k6_2f(6).*y(94,:).*y(45,:)...
    + k6_2r(7).*y(190,:)- k6_2f(7).*y(94,:).*y(55,:)...
    + k6_2r(8).*y(205,:)- k6_2f(8).*y(94,:).*y(65,:)...
    + k6_2r(9).*y(220,:)- k6_2f(9).*y(94,:).*y(75,:)...
    + k6_2r(5).*y(229,:)- k6_2f(5).*y(94,:).*y(40,:)...
    + k6_2r(6).*y(244,:)- k6_2f(6).*y(94,:).*y(50,:)...
    + k6_2r(7).*y(259,:)- k6_2f(7).*y(94,:).*y(60,:)...
    + k6_2r(8).*y(274,:)- k6_2f(8).*y(94,:).*y(70,:)...
    + k6_2r(9).*y(289,:)- k6_2f(9).*y(94,:).*y(80,:)...
;
dy(100,:) =  k6_2f(1).*y(94,:).*y(14,:) - k6_2r(1).*y(100,:) - kcat6(1).*y(100,:);
dy(115,:) =  k6_2f(2).*y(94,:).*y(19,:) - k6_2r(2).*y(115,:) - kcat6(2).*y(115,:);
dy(130,:) =  k6_2f(3).*y(94,:).*y(24,:) - k6_2r(3).*y(130,:) - kcat6(3).*y(130,:);
dy(145,:) = k6_2f(4).*y(94,:).*y(29,:)- k6_2r(4).*y(145,:)- kcat6(4).*y(145,:);
dy(160,:) = k6_2f(5).*y(94,:).*y(35,:)- k6_2r(5).*y(160,:)- kcat6(5).*y(160,:);
dy(175,:) = k6_2f(6).*y(94,:).*y(45,:)- k6_2r(6).*y(175,:)- kcat6(6).*y(175,:);
dy(190,:) = k6_2f(7).*y(94,:).*y(55,:)- k6_2r(7).*y(190,:)- kcat6(7).*y(190,:);
dy(205,:) = k6_2f(8).*y(94,:).*y(65,:)- k6_2r(8).*y(205,:)- kcat6(8).*y(205,:);
dy(220,:) = k6_2f(9).*y(94,:).*y(75,:)- k6_2r(9).*y(220,:)- kcat6(9).*y(220,:);

dy(229,:) = k6_2f(5).*y(94,:).*y(40,:)- k6_2r(5).*y(229,:)- kcat6(5).*y(229,:);
dy(244,:) = k6_2f(6).*y(94,:).*y(50,:)- k6_2r(6).*y(244,:)- kcat6(6).*y(244,:);
dy(259,:) = k6_2f(7).*y(94,:).*y(60,:)- k6_2r(7).*y(259,:)- kcat6(7).*y(259,:);
dy(274,:) = k6_2f(8).*y(94,:).*y(70,:)- k6_2r(8).*y(274,:)- kcat6(8).*y(274,:);
dy(289,:) = k6_2f(9).*y(94,:).*y(80,:)- k6_2r(9).*y(289,:)- kcat6(9).*y(289,:);

dy(101,:) =  k7_1f(1).*e7.*y(15,:) - k7_1r(1).*y(101,:) - kcat7(1).*y(101,:);
dy(116,:) =  k7_1f(2).*e7.*y(20,:) - k7_1r(2).*y(116,:) - kcat7(2).*y(116,:);
dy(131,:) =  k7_1f(3).*e7.*y(25,:) - k7_1r(3).*y(131,:) - kcat7(3).*y(131,:);
dy(146,:) = k7_1f(4).*e7.*y(30,:)- k7_1r(4).*y(146,:)- kcat7(4).*y(146,:);
dy(161,:) = k7_1f(5).*e7.*y(36,:)- k7_1r(5).*y(161,:)- kcat7(5).*y(161,:);
dy(176,:) = k7_1f(6).*e7.*y(46,:)- k7_1r(6).*y(176,:)- kcat7(6).*y(176,:);
dy(191,:) = k7_1f(7).*e7.*y(56,:)- k7_1r(7).*y(191,:)- kcat7(7).*y(191,:);
dy(206,:) = k7_1f(8).*e7.*y(66,:)- k7_1r(8).*y(206,:)- kcat7(8).*y(206,:);
dy(221,:) = k7_1f(9).*e7.*y(76,:)- k7_1r(9).*y(221,:)- kcat7(9).*y(221,:);

dy(230,:) = k7_1f(5).*e7.*y(41,:)- k7_1r(5).*y(230,:)- kcat7(5).*y(230,:);
dy(245,:) = k7_1f(6).*e7.*y(51,:)- k7_1r(6).*y(245,:)- kcat7(6).*y(245,:);
dy(260,:) = k7_1f(7).*e7.*y(61,:)- k7_1r(7).*y(260,:)- kcat7(7).*y(260,:);
dy(275,:) = k7_1f(8).*e7.*y(71,:)- k7_1r(8).*y(275,:)- kcat7(8).*y(275,:);
dy(290,:) = k7_1f(9).*e7.*y(81,:)- k7_1r(9).*y(290,:)- kcat7(9).*y(290,:);

%FabF
dy(102,:) =  k8_1f(1).*e8.*y(15,:) + k8_2r(1).*y(103,:).*y(4,:) - k8_1r(1).*y(102,:) - k8_2f(1).*y(102,:);
dy(117,:) =  k8_1f(2).*e8.*y(20,:) + k8_2r(2).*y(118,:).*y(4,:) - k8_1r(2).*y(117,:) - k8_2f(2).*y(117,:);
dy(132,:) =  k8_1f(3).*e8.*y(25,:) + k8_2r(3).*y(133,:).*y(4,:) - k8_1r(3).*y(132,:) - k8_2f(3).*y(132,:);
dy(147,:) = k8_1f(4).*e8.*y(30,:)+ k8_2r(4).*y(148,:).*y(4,:)- k8_1r(4).*y(147,:)- k8_2f(4).*y(147,:);
dy(162,:) = k8_1f(5).*e8.*y(36,:)+ k8_2r(5).*y(163,:).*y(4,:)- k8_1r(5).*y(162,:)- k8_2f(5).*y(162,:);
dy(177,:) = k8_1f(6).*e8.*y(46,:)+ k8_2r(6).*y(178,:).*y(4,:)- k8_1r(6).*y(177,:)- k8_2f(6).*y(177,:);
dy(192,:) = k8_1f(7).*e8.*y(56,:)+ k8_2r(7).*y(193,:).*y(4,:)- k8_1r(7).*y(192,:)- k8_2f(7).*y(192,:);
dy(207,:) = k8_1f(8).*e8.*y(66,:)+ k8_2r(8).*y(208,:).*y(4,:)- k8_1r(8).*y(207,:)- k8_2f(8).*y(207,:);

dy(231,:) = k8_1f(5).*e8.*y(41,:)+ k8_2r(5).*y(232,:).*y(4,:)- k8_1r(5).*y(231,:)- k8_2f(5).*y(231,:);
dy(246,:) = k8_1f(6).*e8.*y(51,:)+ k8_2r(6).*y(247,:).*y(4,:)- k8_1r(6).*y(246,:)- k8_2f(6).*y(246,:);
dy(261,:) = k8_1f(7).*e8.*y(61,:)+ k8_2r(7).*y(262,:).*y(4,:)- k8_1r(7).*y(261,:)- k8_2f(7).*y(261,:);
dy(276,:) = k8_1f(8).*e8.*y(71,:)+ k8_2r(8).*y(277,:).*y(4,:)- k8_1r(8).*y(276,:)- k8_2f(8).*y(276,:);


dy(103,:) =  k8_2f(1).*y(102,:) + k8_3r(1).*y(104,:) - k8_2r(1).*y(103,:).*y(4,:) - k8_3f(1).*y(103,:).*y(10,:);
dy(118,:) =  k8_2f(2).*y(117,:) + k8_3r(2).*y(119,:) - k8_2r(2).*y(118,:).*y(4,:) - k8_3f(2).*y(118,:).*y(10,:);
dy(133,:) =  k8_2f(3).*y(132,:) + k8_3r(3).*y(134,:) - k8_2r(3).*y(133,:).*y(4,:) - k8_3f(3).*y(133,:).*y(10,:);
dy(148,:) = k8_2f(4).*y(147,:)+ k8_3r(4).*y(149,:)- k8_2r(4).*y(148,:).*y(4,:)- k8_3f(4).*y(148,:).*y(10,:);
dy(163,:) = k8_2f(5).*y(162,:)+ k8_3r(5).*y(164,:)- k8_2r(5).*y(163,:).*y(4,:)- k8_3f(5).*y(163,:).*y(10,:);
dy(178,:) = k8_2f(6).*y(177,:)+ k8_3r(6).*y(179,:)- k8_2r(6).*y(178,:).*y(4,:)- k8_3f(6).*y(178,:).*y(10,:);
dy(193,:) = k8_2f(7).*y(192,:)+ k8_3r(7).*y(194,:)- k8_2r(7).*y(193,:).*y(4,:)- k8_3f(7).*y(193,:).*y(10,:);
dy(208,:) = k8_2f(8).*y(207,:)+ k8_3r(8).*y(209,:)- k8_2r(8).*y(208,:).*y(4,:)- k8_3f(8).*y(208,:).*y(10,:);

dy(232,:) = k8_2f(5).*y(231,:)+ k8_3r(5).*y(233,:)- k8_2r(5).*y(232,:).*y(4,:)- k8_3f(5).*y(232,:).*y(10,:);
dy(247,:) = k8_2f(6).*y(246,:)+ k8_3r(6).*y(248,:)- k8_2r(6).*y(247,:).*y(4,:)- k8_3f(6).*y(247,:).*y(10,:);
dy(262,:) = k8_2f(7).*y(261,:)+ k8_3r(7).*y(263,:)- k8_2r(7).*y(262,:).*y(4,:)- k8_3f(7).*y(262,:).*y(10,:);
dy(277,:) = k8_2f(8).*y(276,:)+ k8_3r(8).*y(278,:)- k8_2r(8).*y(277,:).*y(4,:)- k8_3f(8).*y(277,:).*y(10,:);


dy(104,:) =  k8_3f(1).*y(103,:).*y(10,:) - k8_3r(1).*y(104,:) - kcat8(1).*y(104,:);
dy(119,:) =  k8_3f(2).*y(118,:).*y(10,:) - k8_3r(2).*y(119,:) - kcat8(2).*y(119,:);
dy(134,:) =  k8_3f(3).*y(133,:).*y(10,:) - k8_3r(3).*y(134,:) - kcat8(3).*y(134,:);
dy(149,:) = k8_3f(4).*y(148,:).*y(10,:)- k8_3r(4).*y(149,:)- kcat8(4).*y(149,:);
dy(164,:) = k8_3f(5).*y(163,:).*y(10,:)- k8_3r(5).*y(164,:)- kcat8(5).*y(164,:);
dy(179,:) = k8_3f(6).*y(178,:).*y(10,:)- k8_3r(6).*y(179,:)- kcat8(6).*y(179,:);
dy(194,:) = k8_3f(7).*y(193,:).*y(10,:)- k8_3r(7).*y(194,:)- kcat8(7).*y(194,:);
dy(209,:) = k8_3f(8).*y(208,:).*y(10,:)- k8_3r(8).*y(209,:)- kcat8(8).*y(209,:);

dy(233,:) = k8_3f(5).*y(232,:).*y(10,:)- k8_3r(5).*y(233,:)- kcat8_un(5).*y(233,:);
dy(248,:) = k8_3f(6).*y(247,:).*y(10,:)- k8_3r(6).*y(248,:)- kcat8_un(6).*y(248,:);
dy(263,:) = k8_3f(7).*y(262,:).*y(10,:)- k8_3r(7).*y(263,:)- kcat8_un(7).*y(263,:);
dy(278,:) = k8_3f(8).*y(277,:).*y(10,:)- k8_3r(8).*y(278,:)- kcat8_un(8).*y(278,:);


%FabB
dy(107,:) =  k10_1f(1).*e10.*y(15,:) + k10_2r(1).*y(108,:).*y(4,:) - k10_1r(1).*y(107,:) - k10_2f(1).*y(107,:);
dy(122,:) =  k10_1f(2).*e10.*y(20,:) + k10_2r(2).*y(123,:).*y(4,:) - k10_1r(2).*y(122,:) - k10_2f(2).*y(122,:);
dy(137,:) =  k10_1f(3).*e10.*y(25,:) + k10_2r(3).*y(138,:).*y(4,:) - k10_1r(3).*y(137,:) - k10_2f(3).*y(137,:);
dy(152,:) = k10_1f(4).*e10.*y(30,:)+ k10_2r(4).*y(153,:).*y(4,:)- k10_1r(4).*y(152,:)- k10_2f(4).*y(152,:);
dy(167,:) = k10_1f(5).*e10.*y(36,:)+ k10_2r(5).*y(168,:).*y(4,:)- k10_1r(5).*y(167,:)- k10_2f(5).*y(167,:);
dy(182,:) = k10_1f(6).*e10.*y(46,:)+ k10_2r(6).*y(183,:).*y(4,:)- k10_1r(6).*y(182,:)- k10_2f(6).*y(182,:);
dy(197,:) = k10_1f(7).*e10.*y(56,:)+ k10_2r(7).*y(198,:).*y(4,:)- k10_1r(7).*y(197,:)- k10_2f(7).*y(197,:);
dy(212,:) = k10_1f(8).*e10.*y(66,:)+ k10_2r(8).*y(213,:).*y(4,:)- k10_1r(8).*y(212,:)- k10_2f(8).*y(212,:);

dy(236,:) = k10_1f(5).*e10.*y(41,:)+ k10_2r(5).*y(237,:).*y(4,:)- k10_1r(5).*y(236,:)- k10_2f(5).*y(236,:);
dy(251,:) = k10_1f(6).*e10.*y(51,:)+ k10_2r(6).*y(252,:).*y(4,:)- k10_1r(6).*y(251,:)- k10_2f(6).*y(251,:);
dy(266,:) = k10_1f(7).*e10.*y(61,:)+ k10_2r(7).*y(267,:).*y(4,:)- k10_1r(7).*y(266,:)- k10_2f(7).*y(266,:);
dy(281,:) = k10_1f(8).*e10.*y(71,:)+ k10_2r(8).*y(282,:).*y(4,:)- k10_1r(8).*y(281,:)- k10_2f(8).*y(281,:);


dy(108,:) =  k10_2f(1).*y(107,:) + k10_3r(1).*y(109,:) - k10_2r(1).*y(108,:).*y(4,:) - k10_3f(1).*y(108,:).*y(10,:);
dy(123,:) =  k10_2f(2).*y(122,:) + k10_3r(2).*y(124,:) - k10_2r(2).*y(123,:).*y(4,:) - k10_3f(2).*y(123,:).*y(10,:);
dy(138,:) =  k10_2f(3).*y(137,:) + k10_3r(3).*y(139,:) - k10_2r(3).*y(138,:).*y(4,:) - k10_3f(3).*y(138,:).*y(10,:);
dy(153,:) = k10_2f(4).*y(152,:)+ k10_3r(4).*y(154,:)- k10_2r(4).*y(153,:).*y(4,:)- k10_3f(4).*y(153,:).*y(10,:);
dy(168,:) = k10_2f(5).*y(167,:)+ k10_3r(5).*y(169,:)- k10_2r(5).*y(168,:).*y(4,:)- k10_3f(5).*y(168,:).*y(10,:);
dy(183,:) = k10_2f(6).*y(182,:)+ k10_3r(6).*y(184,:)- k10_2r(6).*y(183,:).*y(4,:)- k10_3f(6).*y(183,:).*y(10,:);
dy(198,:) = k10_2f(7).*y(197,:)+ k10_3r(7).*y(199,:)- k10_2r(7).*y(198,:).*y(4,:)- k10_3f(7).*y(198,:).*y(10,:);
dy(213,:) = k10_2f(8).*y(212,:)+ k10_3r(8).*y(214,:)- k10_2r(8).*y(213,:).*y(4,:)- k10_3f(8).*y(213,:).*y(10,:);

dy(237,:) = k10_2f(5).*y(236,:)+ k10_3r(5).*y(238,:)- k10_2r(5).*y(237,:).*y(4,:)- k10_3f(5).*y(237,:).*y(10,:);
dy(252,:) = k10_2f(6).*y(251,:)+ k10_3r(6).*y(253,:)- k10_2r(6).*y(252,:).*y(4,:)- k10_3f(6).*y(252,:).*y(10,:);
dy(267,:) = k10_2f(7).*y(266,:)+ k10_3r(7).*y(268,:)- k10_2r(7).*y(267,:).*y(4,:)- k10_3f(7).*y(267,:).*y(10,:);
dy(282,:) = k10_2f(8).*y(281,:)+ k10_3r(8).*y(283,:)- k10_2r(8).*y(282,:).*y(4,:)- k10_3f(8).*y(282,:).*y(10,:);


dy(109,:) =  k10_3f(1).*y(108,:).*y(10,:) - k10_3r(1).*y(109,:) - kcat10(1).*y(109,:);
dy(124,:) =  k10_3f(2).*y(123,:).*y(10,:) - k10_3r(2).*y(124,:) - kcat10(2).*y(124,:);
dy(139,:) =  k10_3f(3).*y(138,:).*y(10,:) - k10_3r(3).*y(139,:) - kcat10(3).*y(139,:);
dy(154,:) = k10_3f(4).*y(153,:).*y(10,:)- k10_3r(4).*y(154,:)- kcat10(4).*y(154,:);
dy(169,:) = k10_3f(5).*y(168,:).*y(10,:)- k10_3r(5).*y(169,:)- kcat10(5).*y(169,:);
dy(184,:) = k10_3f(6).*y(183,:).*y(10,:)- k10_3r(6).*y(184,:)- kcat10(6).*y(184,:);
dy(199,:) = k10_3f(7).*y(198,:).*y(10,:)- k10_3r(7).*y(199,:)- kcat10(7).*y(199,:);
dy(214,:) = k10_3f(8).*y(213,:).*y(10,:)- k10_3r(8).*y(214,:)- kcat10(8).*y(214,:);

dy(238,:) = k10_3f(5).*y(237,:).*y(10,:)- k10_3r(5).*y(238,:)- kcat10_un(5).*y(238,:);
dy(253,:) = k10_3f(6).*y(252,:).*y(10,:)- k10_3r(6).*y(253,:)- kcat10_un(6).*y(253,:);
dy(268,:) = k10_3f(7).*y(267,:).*y(10,:)- k10_3r(7).*y(268,:)- kcat10_un(7).*y(268,:);
dy(283,:) = k10_3f(8).*y(282,:).*y(10,:)- k10_3r(8).*y(283,:)- kcat10_un(8).*y(283,:);


dy(302,:) = k10_1f(4).*e10.*y(32,:)+ k10_2r(4).*y(303,:).*y(4,:)- k10_1r(4).*y(302,:)- k10_2f(4).*y(302,:);
dy(303,:) = k10_2f(4).*y(302,:)+ k10_3r(4).*y(304,:)- k10_2r(4).*y(303,:).*y(4,:)- k10_3f(4).*y(303,:).*y(10,:);
dy(304,:) = k10_3f(4).*y(303,:).*y(10,:)- k10_3r(4).*y(304,:)- kcat10_un(4).*y(304,:);


dy(305,:) = kcat5_un(4).*y(144,:) - k5_2r_un(4).*y(305,:) - k5_3f_un(4).*y(305,:) + k5_3r_un(4).*e5.*y(32,:);

dy(306,:) = k8_1f(4).*e8.*y(32,:)+ k8_2r(4).*y(307,:).*y(4,:)- k8_1r(4).*y(306,:)- k8_2f(4).*y(306,:);
dy(307,:) = k8_2f(4).*y(306,:)+ k8_3r(4).*y(308,:)- k8_2r(4).*y(307,:).*y(4,:)- k8_3f(4).*y(307,:).*y(10,:);
dy(308,:) = k8_3f(4).*y(307,:).*y(10,:)- k8_3r(4).*y(308,:)- kcat8_un(4).*y(308,:);

dy(95,:)  = k3_4f(1).*e3.*y(15,:) - k3_4r(1).*y(95,:);
dy(110,:)  = k3_4f(2).*e3.*y(20,:) - k3_4r(2).*y(110,:);
dy(125,:)  = k3_4f(3).*e3.*y(25,:) - k3_4r(3).*y(125,:);
dy(140,:) = k3_4f(4).*e3.*y(30,:)- k3_4r(4).*y(140,:);
dy(155,:) = k3_4f(5).*e3.*y(36,:)- k3_4r(5).*y(155,:);
dy(170,:) = k3_4f(6).*e3.*y(46,:)- k3_4r(6).*y(170,:);
dy(185,:) = k3_4f(7).*e3.*y(56,:)- k3_4r(7).*y(185,:);
dy(200,:) = k3_4f(8).*e3.*y(66,:)- k3_4r(8).*y(200,:);
dy(215,:) = k3_4f(9).*e3.*y(76,:)- k3_4r(9).*y(215,:);

dy(224,:) = k3_4f(5).*e3.*y(41,:)- k3_4r(5).*y(224,:);
dy(239,:) = k3_4f(6).*e3.*y(51,:)- k3_4r(6).*y(239,:);
dy(254,:) = k3_4f(7).*e3.*y(61,:)- k3_4r(7).*y(254,:);
dy(269,:) = k3_4f(8).*e3.*y(71,:)- k3_4r(8).*y(269,:);
dy(284,:) = k3_4f(9).*e3.*y(81,:)- k3_4r(9).*y(284,:);

dy(96,:)  = k3_5f(1).*y(91,:).*y(15,:)  - k3_5r(1).*y(96,:);
dy(111,:)  = k3_5f(2).*y(91,:).*y(20,:)  - k3_5r(2).*y(111,:);
dy(126,:)  = k3_5f(3).*y(91,:).*y(25,:)  - k3_5r(3).*y(126,:);
dy(141,:) = k3_5f(4).*y(91,:).*y(30,:) - k3_5r(4).*y(141,:);
dy(156,:) = k3_5f(5).*y(91,:).*y(36,:) - k3_5r(5).*y(156,:);
dy(171,:) = k3_5f(6).*y(91,:).*y(46,:) - k3_5r(6).*y(171,:);
dy(186,:) = k3_5f(7).*y(91,:).*y(56,:) - k3_5r(7).*y(186,:);
dy(201,:) = k3_5f(8).*y(91,:).*y(66,:) - k3_5r(8).*y(201,:);
dy(216,:) = k3_5f(9).*y(91,:).*y(76,:) - k3_5r(9).*y(216,:);

dy(225,:) = k3_5f(5).*y(91,:).*y(41,:) - k3_5r(5).*y(225,:);
dy(240,:) = k3_5f(6).*y(91,:).*y(51,:) - k3_5r(6).*y(240,:);
dy(255,:) = k3_5f(7).*y(91,:).*y(61,:) - k3_5r(7).*y(255,:);
dy(270,:) = k3_5f(8).*y(91,:).*y(71,:) - k3_5r(8).*y(270,:);
dy(285,:) = k3_5f(9).*y(91,:).*y(81,:) - k3_5r(9).*y(285,:);


%binding of ACP to TesA
dy(297,:) = k7_inh_f.*e7.*y(4,:) - k7_inh_r.*y(297,:);


%binding of ACP to FabH FabG FabZ FabI FabF
dy(293,:) = k3_inh_f.*e3.*y(4,:) - k3_inh_r.*y(293,:);
dy(294,:) = k4_inh_f.*e4.*y(4,:) - k4_inh_r.*y(294,:);
dy(295,:) = k5_inh_f.*e5.*y(4,:) - k5_inh_r.*y(295,:);
dy(296,:) = k6_inh_f.*e6.*y(4,:) - k6_inh_r.*y(296,:);
dy(298,:) = k8_inh_f.*e8.*y(4,:) - k8_inh_r.*y(298,:);

%binding of ACP to FabA FabB
dy(299,:) = k9_inh_f.*e9.*y(4,:) - k9_inh_r.*y(299,:);
dy(300,:) = k10_inh_f.*e10.*y(4,:) - k10_inh_r.*y(300,:);


%LipidA Reactions
dy(309,:)=-k15.*y(309,:).*y(322,:)-k16.*y(309,:).*y(313,:)+k18.*y(319,:)+k19.*(FtsH_init-y(309,:)-y(319,:));
dy(310,:)=k1.*promote-k11.*y(310,:)-k14.*y(310,:).*y(319,:);
dy(311,:)=-LpxH.*A_kcat3.*y(311,:)./((Km1+y(311,:)).*(1+y(316,:)./Ki))-LpxB.*A_kcat4.*y(311,:).*y(316,:)./(Kma2.*Kmb2+y(311,:).*Kmb2+y(316,:).*Kma2+y(311,:).*y(316,:))+LpxD.*A_kcat10./(1+y(311,:)./ki2).*y(318,:).*y(44,:)./(Kma4.*Kmb4+y(318,:).*Kmb4+y(44,:).*Kma4+y(318,:).*y(44,:));
dy(312,:)=k12.*promot-k13.*y(312,:)-k17.*y(312,:).*(FtsH_init-y(309,:)-y(319,:));
dy(313,:)=A_kcat8.*LpxM.*y(317,:).*y(46,:)./((Km4+y(317,:).*(km_myr+y(46,:))))-A_kcat9.*MsbA.*y(313,:)./(Km5+y(313,:))-y(312,:).*A_kcat11.*CMP_KDO.*y(313,:)./(Kma5.*Kmb5+CMP_KDO.*Kmb5+y(313,:).*Kma5+CMP_KDO.*y(313,:));
dy(314,:)=LpxA.*A_kcat.*UDP_Glc.*y(44,:)./(Kma.*Kmb+UDP_Glc.*Kmb+y(44,:).*Kma+UDP_Glc.*y(44,:))-LpxA.*A_kcat1.*y(314,:).*y(4,:)./(Kma1.*Kmb1+y(314,:).*Kmb1+y(4,:).*Kma1+y(314,:).*y(4,:))-A_kcat2.*y(310,:).*y(314,:)./(Km+y(314,:));
dy(315,:)=A_kcat5.*LpxK.*y(322,:)./(Km2+y(322,:))-y(312,:).*A_kcat6./(1+y(313,:)./ki1).*CMP_KDO.*y(315,:)./(Kma3.*Kmb3+CMP_KDO.*Kmb3+y(315,:).*Kma3+CMP_KDO.*y(315,:));
dy(316,:)=LpxH.*A_kcat3.*y(311,:)./((Km1+y(311,:)).*(1+y(316,:)./Ki))-LpxB.*A_kcat4.*y(311,:).*y(316,:)./(Kma2.*Kmb2+y(311,:).*Kmb2+y(316,:).*Kma2+y(311,:).*y(316,:));
dy(317,:)=A_kcat7.*LpxL.*y(320,:).*y(36,:)./((Km3+y(320,:).*(km_laur+y(36,:))))-A_kcat8.*LpxM.*y(317,:).*y(46,:)./((Km4+y(317,:).*(km_myr+y(46,:))));
dy(318,:)=A_kcat2.*y(310,:).*y(314,:)./(Km+y(314,:))-LpxD.*A_kcat10./(1+y(311,:)./ki2).*y(318,:).*y(44,:)./(Kma4.*Kmb4+y(318,:).*Kmb4+y(44,:).*Kma4+y(318,:).*y(44,:));
dy(319,:)=k15.*y(309,:).*y(322,:)-k18.*y(319,:);
dy(320,:)=y(312,:).*A_kcat6./(1+y(313,:)./ki1).*CMP_KDO.*y(315,:)./(Kma3.*Kmb3+CMP_KDO.*Kmb3+y(315,:).*Kma3+CMP_KDO.*y(315,:))-A_kcat7.*LpxL.*y(320,:).*y(36,:)./((Km3+y(320,:).*(km_laur+y(36,:))));
dy(321,:)=A_kcat9.*MsbA.*y(313,:)./(Km5+y(313,:));
dy(322,:)=LpxB.*A_kcat4.*y(311,:).*y(316,:)./(Kma2.*Kmb2+y(311,:).*Kmb2+y(316,:).*Kma2+y(311,:).*y(316,:))-A_kcat5.*LpxK.*y(322,:)./(Km2+y(322,:));
dy(323,:)=y(312,:).*A_kcat11.*CMP_KDO.*y(313,:)./(Kma5.*Kmb5+CMP_KDO.*Kmb5+y(313,:).*Kma5+CMP_KDO.*y(313,:));


%Phospholipid Reactions
dy(324,:) = kcatP1.*PlsB.*G3P.*y(56,:)./((KmP1+G3P).*(KmP2+y(56,:))) + kcatP1.*PlsB.*G3P.*y(71,:)./((KmP1+G3P).*(KmP3+y(71,:)))...
            -kcatP2*PlsC.*y(324,:).*y(61,:)./((KmP4 + y(324,:)).*(KmP5 + y(61,:))) - kcatP2*PlsC.*y(324,:).*y(71,:)./((KmP4 + y(324,:)).*(KmP6 + y(71,:)));
dy(325,:) = kcatP2*PlsC.*y(324,:).*y(61,:)./((KmP4 + y(324,:)).*(KmP5 + y(61,:))) + kcatP2*PlsC.*y(324,:).*y(71,:)./((KmP4 + y(324,:)).*(KmP6 + y(71,:)));

%FabF and FabB New Initiation reactions.
dy(326,:) = -k8_8f.*e8.*y(326,:) + k8_8r.*y(328,:) - k10_8f.*e10.*y(326,:) + k10_8r.*y(333,:);
dy(327,:) = k8_7f.*e8.*y(10,:) - k8_7r.*y(327,:) - kcat8_CO2.*y(327,:);
dy(328,:) = k8_8f.*e8.*y(326,:) - k8_8r.*y(328,:) - k8_9f.*y(328,:) + k8_9r.*y(329,:).*y(4,:) + kcat8_CO2.*y(327,:);
dy(329,:) = k8_9f.*y(328,:) - k8_9r.*y(329,:).*y(4,:) + k8_5f.*y(330,:) - k8_5r.*y(329,:).*y(9,:) - k8_6f.*y(329,:).*y(10,:) + k8_6r.*y(331,:);
dy(330,:) = k8_4f.*e8.*y(3,:) - k8_4r.*y(330,:) - k8_5f.*y(330,:) + k8_5r.*y(329,:).*y(9,:);
dy(331,:) = k8_6f.*y(329,:).*y(10,:) - k8_6r.*y(331,:) - kcat8_H.*y(331,:);
dy(332,:) = k10_7f.*e10.*y(10,:) - k10_7r.*y(332,:) - kcat10_CO2.*y(332,:);
dy(333,:) = k10_8f.*e10.*y(326,:) - k10_8r.*y(333,:) - k10_9f.*y(333,:) + k10_9r.*y(334,:).*y(4) + kcat10_CO2.*y(332,:);
dy(334,:) = k10_9f.*y(333,:) - k10_9r.*y(334,:).*y(4,:) + k10_5f.*y(335,:) - k10_5r.*y(334,:).*y(9,:) - k10_6f.*y(334,:).*y(10,:) + k10_6r.*y(336,:);
dy(335,:) = k10_4f.*e10.*y(3,:) - k10_4r.*y(335,:) - k10_5f.*y(335,:) + k10_5r.*y(334,:).*y(9,:);
dy(336,:) = k10_6f.*y(334,:).*y(10,:) - k10_6r.*y(336,:) - kcat10_H.*y(336,:);


%Fatty Alcohol Production Pathway (FadD and ACR2)
dy(337,:) = kcatD1(1).*FadD.*y(1,:).*y(21,:)./((Km_ATP+y(1,:)).*(Km_FA(1)+y(21,:))) - kcatD2.*FadD.*y(9,:).*y(337,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(1)+y(337,:)));%6
dy(338,:) = kcatD1(2).*FadD.*y(1,:).*y(26,:)./((Km_ATP+y(1,:)).*(Km_FA(2)+y(26,:))) - kcatD2.*FadD.*y(9,:).*y(338,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(2)+y(338,:)));%8
dy(339,:) = kcatD1(3).*FadD.*y(1,:).*y(31,:)./((Km_ATP+y(1,:)).*(Km_FA(3)+y(31,:))) - kcatD2*FadD.*y(9,:).*y(339,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(3)+y(339,:)));%10
dy(340,:) = kcatD1(4).*FadD.*y(1,:).*y(37,:)./((Km_ATP+y(1,:)).*(Km_FA(4)+y(37,:))) - kcatD2.*FadD.*y(9,:).*y(340,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(4)+y(340,:)));%12
dy(341,:) = kcatD1(5).*FadD.*y(1,:).*y(47,:)./((Km_ATP+y(1,:)).*(Km_FA(5)+y(47,:))) - kcatD2.*FadD.*y(9,:).*y(341,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(5)+y(341,:)));%14
dy(342,:) = kcatD1(6).*FadD.*y(1,:).*y(57,:)./((Km_ATP+y(1,:)).*(Km_FA(6)+y(57,:))) - kcatD2.*FadD.*y(9,:).*y(342,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(6)+y(342,:)));%16
dy(343,:) = kcatD1(7).*FadD.*y(1,:).*y(67,:)./((Km_ATP+y(1,:)).*(Km_FA(7)+y(67,:))) - kcatD2.*FadD.*y(9,:).*y(343,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(7)+y(343,:)));%18
dy(344,:) = kcatD1(8).*FadD.*y(1,:).*y(42,:)./((Km_ATP+y(1,:)).*(Km_FA(8)+y(42,:))) - kcatD2.*FadD.*y(9,:).*y(344,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(8)+y(344,:)));%12_unsat
dy(345,:) = kcatD1(9).*FadD.*y(1,:).*y(52,:)./((Km_ATP+y(1,:)).*(Km_FA(9)+y(52,:))) - kcatD2.*FadD.*y(9,:).*y(345,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(9)+y(345,:)));%14_unsat
dy(346,:) = kcatD1(10).*FadD.*y(1,:).*y(62,:)./((Km_ATP+y(1,:)).*(Km_FA(10)+y(62,:))) -kcatD2.*FadD.*y(9,:).*y(346,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(10)+y(346,:)));%16_unsat
dy(347,:) = kcatD1(11).*FadD.*y(1,:).*y(72,:)./((Km_ATP+y(1,:)).*(Km_FA(11)+y(72,:))) -kcatD2.*FadD.*y(9,:).*y(347,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(11)+y(347,:)));%18_unsat

%FadD Rxn 2 (rate of formation of FA-CoA)
dy(348,:) = kcatD2.*FadD.*y(9,:).*y(337,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(1)+y(337,:))) - kcatACR1(1).*ACR1.*y(348,:)./(Km_ACR1 + y(348,:));%6
dy(349,:) = kcatD2.*FadD.*y(9,:).*y(338,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(2)+y(338,:))) - kcatACR1(2).*ACR1.*y(349,:)./(Km_ACR1 + y(349,:));%8
dy(350,:) = kcatD2.*FadD.*y(9,:).*y(339,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(3)+y(339,:))) - kcatACR1(3).*ACR1.*y(350,:)./(Km_ACR1 + y(350,:));%10
dy(351,:) = kcatD2.*FadD.*y(9,:).*y(340,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(4)+y(340,:))) - kcatACR1(4).*ACR1.*y(351,:)./(Km_ACR1 + y(351,:));%12
dy(352,:) = kcatD2.*FadD.*y(9,:).*y(341,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(5)+y(341,:))) - kcatACR1(5).*ACR1.*y(352,:)./(Km_ACR1 + y(352,:));%14
dy(353,:) = kcatD2.*FadD.*y(9,:).*y(342,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(6)+y(342,:))) - kcatACR1(6).*ACR1.*y(353,:)./(Km_ACR1 + y(353,:));%16
dy(354,:) = kcatD2.*FadD.*y(9,:).*y(343,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(7)+y(343,:))) - kcatACR1(7).*ACR1.*y(354,:)./(Km_ACR1 + y(354,:));%18
dy(355,:) = kcatD2.*FadD.*y(9,:).*y(344,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(8)+y(344,:))) - kcatACR1(8).*ACR1.*y(355,:)./(Km_ACR1 + y(355,:));%12_unsat
dy(356,:) = kcatD2.*FadD.*y(9,:).*y(345,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(9)+y(345,:))) - kcatACR1(9).*ACR1.*y(356,:)./(Km_ACR1 + y(356,:));%14_unsat
dy(357,:) = kcatD2.*FadD.*y(9,:).*y(346,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(10)+y(346,:)))- kcatACR1(10).*ACR1.*y(357,:)./(Km_ACR1 + y(357,:));%16_unsat
dy(358,:) = kcatD2.*FadD.*y(9,:).*y(347,:)./((Km_CoA+y(9,:)).*(Km_FA_AMP(11)+y(347,:)))- kcatACR1(11).*ACR1.*y(358,:)./(Km_ACR1 + y(358,:));%18_unsat

%Fatty Aldehyde Production
dy(359,:) = kcatACR1(1).*ACR1.*y(348,:)./(Km_ACR1 + y(348,:)) - kcatAhr(2).*Ahr.*y(6,:).*y(359,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(2)+y(359,:)));%6
dy(360,:) = kcatACR1(2).*ACR1.*y(349,:)./(Km_ACR1 + y(349,:)) - kcatAhr(3).*Ahr.*y(6,:).*y(360,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(3)+y(360,:)));%8
dy(361,:) = kcatACR1(3).*ACR1.*y(350,:)./(Km_ACR1 + y(350,:)) - kcatAhr(4).*Ahr.*y(6,:).*y(361,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(4)+y(361,:)));%10
dy(362,:) = kcatACR1(4).*ACR1.*y(351,:)./(Km_ACR1 + y(351,:)) - kcatAhr(5).*Ahr.*y(6,:).*y(362,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(5)+y(362,:)));%12
dy(363,:) = kcatACR1(5).*ACR1.*y(352,:)./(Km_ACR1 + y(352,:)) - kcatAhr(6).*Ahr.*y(6,:).*y(363,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(6)+y(363,:)));%14
dy(364,:) = kcatACR1(6).*ACR1.*y(353,:)./(Km_ACR1 + y(353,:)) - kcatAhr(7).*Ahr.*y(6,:).*y(364,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(7)+y(364,:)));%16
dy(365,:) = kcatACR1(7).*ACR1.*y(354,:)./(Km_ACR1 + y(354,:)) - kcatAhr(8).*Ahr.*y(6,:).*y(365,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(8)+y(365,:)));%18
dy(366,:) = kcatACR1(8).*ACR1.*y(355,:)./(Km_ACR1 + y(355,:)) - kcatAhr(9).*Ahr.*y(6,:).*y(366,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(9)+y(366,:)));%12_unsat
dy(367,:) = kcatACR1(9).*ACR1.*y(356,:)./(Km_ACR1 + y(356,:)) - kcatAhr(10).*Ahr.*y(6,:).*y(367,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(10)+y(367,:)));%14_unsat
dy(368,:) = kcatACR1(10).*ACR1.*y(357,:)./(Km_ACR1 + y(357,:)) - kcatAhr(11).*Ahr.*y(6,:).*y(368,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(11)+y(368,:)));%16_unsat
dy(369,:) = kcatACR1(11).*ACR1.*y(358,:)./(Km_ACR1 + y(358,:)) - kcatAhr(12).*Ahr.*y(6,:).*y(369,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(12)+y(369,:)));%18_unsat

%Fatty Alcohol Production
dy(370,:) = kcatAhr(2).*Ahr.*y(6,:).*y(359,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(2)+y(359,:)));%6
dy(371,:) = kcatAhr(3).*Ahr.*y(6,:).*y(360,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(3)+y(360,:)));%8
dy(372,:) = kcatAhr(4).*Ahr.*y(6,:).*y(361,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(4)+y(361,:)));%10
dy(373,:) = kcatAhr(5).*Ahr.*y(6,:).*y(362,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(5)+y(362,:)));%12
dy(374,:) = kcatAhr(6).*Ahr.*y(6,:).*y(363,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(6)+y(363,:)));%14
dy(375,:) = kcatAhr(7).*Ahr.*y(6,:).*y(364,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(7)+y(364,:)));%16
dy(376,:) = kcatAhr(8).*Ahr.*y(6,:).*y(365,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(8)+y(365,:)));%18
dy(377,:) = kcatAhr(9).*Ahr.*y(6,:).*y(366,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(9)+y(366,:)));%12_unsat
dy(378,:) = kcatAhr(10).*Ahr.*y(6,:).*y(367,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(10)+y(367,:)));%14_unsat
dy(379,:) = kcatAhr(11).*Ahr.*y(6,:).*y(368,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(11)+y(368,:)));%16_unsat
dy(380,:) = kcatAhr(12).*Ahr.*y(6,:).*y(369,:)./((Km_Ahr_NADPH+y(6,:)).*(Km_FA_Ald(12)+y(369,:)));%18_unsat

