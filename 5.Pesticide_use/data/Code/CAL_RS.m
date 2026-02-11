%%%%% This is a matlab script that calculates the pesticide risk score of
%%%%% one grid cell following the description in Tang F.H.M., Lenzen M., McBratney A., & Maggi F. (2021) Risk of pesticide pollution at the global scale, Nature Geoscience.
%%%%% Output: Pesticide risk score
%%%%% Authors: Fiona H.M. Tang, Manfred Lenzen, Alexander McBratney, and Federico Maggi
%%%%% Last updated: 17 February 2020
%%%%% Contact name: Fiona H.M. Tang / Federico Maggi
%%%%% Contact email: fiona.tang@sydney.edu.au / federico.maggi@sydney.edu.au

%%%% ---- Load active ingredient properties
load AI_Properties.mat
%%%% ---- END Load active ingredient properties

%%%% ---- Load geo-specific data
load DATA_px1136_py668
%%%% ---- END Load geo-specific data

%%%% --- Define parameters
KH = Henry * 4.03395e-4; % Henry's constant, convert from Pa m3/mol to diamensionless
koc = KOC * 0.001; % adsorption constant to organic carbon, convert mL/g to m3/kg
Ks = Solubility_H2O ./ MM ; % solubility,[mol/m3] % 1 mg/L = 1 g/m3
VP = VP/1000; % vapour pressure, convert mPa to Pa;
Vf = 324.55 ; % [m/h] dilution velocity
d = 0.005 ; %[m] boundary layer thickness
R = 8.314 ; %[J/(mol.K)] gas constant
Da = 0.036 * (76./MM) .^ 0.5; %[m2/h]
Zw = Ks ./ VP ; % [mol/(m3 Pa)]
DEPTH = 0.02 ; % [m] mixing depth is assumed to be 2 cm
fint = 0; %% crop intercept (assumed worst case scenario)
fdrift = 0; % assume no drift loss
if SLP >=20 %% slope factor calculation (Trevisan et al, 2009)
    f1 = 1;
else
    f1 = 0.02153 * SLPtmp + 0.001423 .* SLPtmp .^2;
end
f2 = 1; % assume all fields are adjacent to water bodies
drink_st = 0.1; % EU drinking water standard [micro g/L - EU]
EXCLUDE = [10 15 69]; % compound to exclude from calculation
%%%% ---- END Define parameters

%%%% --- Calculate PECs and RQs
for ncrop = 1:length(Crop_Name)
    
    for ncomp = 1:length(Comp_Name)
        
        jcomp = Comp_PPDB_ID(ncrop,ncomp);
        RATE = ( APR_H(ncrop,ncomp) + APR_L(ncrop,ncomp)) /2; % median application rate [kg/ha]
        
        if ismember(jcomp,EXCLUDE) == 0
            
            %%% Groundwater PECgw and RQgw following Trevisan et al (2009)
            RATEtmp = RATE * 0.1 ; %kg/ha to g/m2
            RF = 1 + ( (BLD * SOC * koc(jcomp) ) +  (SAC .* KH(jcomp)) ) ./ SFC;
            TR = WTD * RF * SFC / GWR;
            AF = exp(-1*TR * log(2)/DT50(jcomp) );
            PECgw = 2.739 * AF * RATEtmp * (1-fint) * (1-fdrift) / (PHI * GWT);  %[micro-g/L]
            RQgw(ncrop,ncomp)   = PECgw / drink_st ;
            
            %%% Soil PECsoil and RQsoil following Trevisan et al (2009)
            RATEtmp = RATE * 100; % convert kg/ha to mg/m2
            PECsoil = RATEtmp .* (1-fint) .* (1-fdrift) ./ (DEPTH .* BLD); %[mg/kg-soil]
            RQsoil(ncrop,ncomp) =  PECsoil./ (LC_EWorm(jcomp)/1000); % placing an assessment factor of 1000
            
            %%% Surface water PECsw and RQsw following Trevisan et al (2009)
            RATEtmp = RATE .* 100; % convert kg/ha to mg/m2
            fw = exp(-3*log(2)/DT50(jcomp))/ (1 + (koc(jcomp) * SOC * BLD) );
            PECsw = ( RATEtmp * (1-fint) * (1-fdrift) * f1 * f2 * fw ) ./ ( max_Rain/1000 );  %[mg/m3]
            PECsw =  PECsw * 0.001; %convert mg/m3 to mg/L
            RQsw(ncrop,ncomp)    =  PECsw./ (LC_Fish(jcomp)/1000); % placing an assessment factor of 1000
            
            %%% Atmosphere PECair and RQair following Trevisan et al (2009)
            RATEtmp = PECsoil; %[mg/kg-soil];
            Za = 1./ (R * avg_Temp) ; %[mol/(Pa m3)]
            Zs = (koc(jcomp) * SOC * BLD * Zw(jcomp)) ./ (1 - PHI); %[mol/(Pa m3)]
            Vs = 1 - PHI;
            Vw = SFC;
            Va = SAC;
            
            Pa = (Za * Va) / (Za * Va + Zw(jcomp) * Vw + Zs * Vs);
            Csa = (RATEtmp * BLD * Pa) / SAC; % [mg/m3]
            J0 = Da(jcomp) * Csa /d; %[mg/(m2 h)]
            PECair = J0 / Vf * 0.001 ; %[mg/L]
            RQair(ncrop,ncomp) =  PECair / (LC_Rat_IH(jcomp)/1000); % placing an assessment factor of 1000
            
        else
            
            RQgw(ncrop,ncomp)   = 0;
            RQsoil(ncrop,ncomp)   = 0;
            RQsw(ncrop,ncomp)   = 0;
            RQair(ncrop,ncomp)   = 0;
        end
        

    end %ncomp
    
end %ncrop


%%%%% calculate the weighted HQ as a function of a.i (HQ_AI_all.mat)
for ncomp = 1:max(max(Comp_PPDB_ID))
    
    [jcrop, jcomp] =  find(Comp_PPDB_ID == ncomp);
    
    tmpgw = zeros(1,length(jcrop));
    tmpsoil = zeros(1,length(jcrop));
    tmpsw = zeros(1,length(jcrop));
    tmpair = zeros(1,length(jcrop));
    for n = 1:length(jcrop)
        tmpgw(n)   = RQgw(jcrop(n),jcomp(n))* CROSA(jcrop(n)) ;
        tmpsoil(n) = RQsoil(jcrop(n),jcomp(n))* CROSA(jcrop(n)) ;
        tmpsw(n)   = RQsw(jcrop(n),jcomp(n))* CROSA(jcrop(n)) ;
        tmpair(n)  = RQair(jcrop(n),jcomp(n))* CROSA(jcrop(n)) ;
    end
    RQ(1,ncomp)   = sum(tmpgw)/sum(CROSA);
    RQ(2,ncomp) = sum(tmpsoil)/sum(CROSA);
    RQ(3,ncomp)   = sum(tmpsw)/sum(CROSA);
    RQ(4,ncomp)  = sum(tmpair)/sum(CROSA);
    
end

%%%% Calculate Risk score (RS)
RP     = log10(sum(RQ,2));
RS     = max(RP); % Risk score of that grid cell
