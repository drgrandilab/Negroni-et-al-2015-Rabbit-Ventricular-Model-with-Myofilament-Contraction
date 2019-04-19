%% rabbit_myofilament_masterCompute
% This file loads the initial conditions, calls the ode solver, and plots
% the results.
%
% Reference: Jorge A Negroni, Stefano Morotti, Elena C Lascano, Aldrin V 
% Gomes, Eleonora Grandi, José L Puglisi, Donald M Bers. ?-adrenergic 
% effects on cardiac myofilaments and contraction in an integrated rabbit 
% ventricular myocyte model. J Mol Cell Cardiol. 2015.
% 
% Please cite the above paper when using this model.

%%

close all;
clear all; 
clc;
%% Parameters for external modules

% ECC and CaM modules
freq = 1.0;                 % [Hz] CHANGE DEPENDING ON FREQUENCY
cycleLength = 1e3/freq;     % [ms]
CaMtotDyad = 418;           % [uM]
BtotDyad = 1.54/8.293e-4;   % [uM]
CaMKIItotDyad = 120;        % [uM] 
CaNtotDyad = 3e-3/8.293e-4; % [uM] 
PP1totDyad = 96.5;          % [uM]
CaMtotSL = 5.65;            % [uM]
BtotSL = 24.2;              % [uM]
CaMKIItotSL = 120*8.293e-4; % [uM]
CaNtotSL = 3e-3;            % [uM]
PP1totSL = 0.57;            % [uM]
CaMtotCyt = 5.65;           % [uM]
BtotCyt = 24.2;             % [uM]
CaMKIItotCyt = 120*8.293e-4;% [uM]
CaNtotCyt = 3e-3;           % [uM] 
PP1totCyt = 0.57;           % [uM]

% ADJUST CAMKII ACTIVITY LEVELS (expression = 'WT', 'OE', or 'KO')
expression = 'WT';
CKIIOE = 0; % Should be zero during 'WT' and 'KO' runs

if strcmp(expression,'OE')
    CKIIOE = 1; % Flag for CKII OE in ECC file (0=WT, 1=OE) - for Ito and INa
    CaMKIItotDyad = 120*6;          % [uM] 
    CaMKIItotSL = 120*8.293e-4*6;   % [uM]
    CaMKIItotCyt = 120*8.293e-4*6;  % [uM]
elseif strcmp(expression,'KO')
    CaMKIItotDyad = 0;              % [uM] 
    CaMKIItotSL = 0;                % [uM]
    CaMKIItotCyt = 0;               % [uM]
end

% Parameters for CaMKII module
LCCtotDyad = 31.4*.9;       % [uM] - Total Dyadic [LCC] - (umol/l dyad)
LCCtotSL = 0.0846;          % [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
RyRtot = 382.6;             % [uM] - Total RyR (in Dyad)
PP1_dyad = 95.7;            % [uM] - Total dyadic [PP1]
PP1_SL = 0.57;              % [uM] - Total Subsarcolemmal [PP1]
PP2A_dyad = 95.76;          % [uM] - Total dyadic PP2A
OA = 0;                     % [uM] - PP1/PP2A inhibitor Okadaic Acid
PLBtot = 38;                % [uM] - Total [PLB] in cytosolic units

% Parameters for BAR module
Ligtot = 0%; 0.1 or 0.02    % [uM] - SET LIGAND CONCENTRATION HERE
LCCtotBA = 0.025;           % [uM] - [umol/L cytosol]
RyRtotBA = 0.135;           % [uM] - [umol/L cytosol]
PLBtotBA = PLBtot;          % [uM] - [umol/L cytosol]
TnItotBA = 70;              % [uM] - [umol/L cytosol]
IKstotBA = 0.025;           % [uM] - [umol/L cytosol]
ICFTRtotBA = 0.025;         % [uM] - [umol/L cytosol]
PP1_PLBtot = 0.89;          % [uM] - [umol/L cytosol]
PLMtotBA = 48;              % [uM] - [umol/L cytosol]
MyototBA = 70;              % [uM] - [umol/L cytosol]
IKrtotBA = 0.025;           % [uM] - [umol/L cytosol]
IClCatotBA = 0.025;         % [uM] - [umol/L cytosol]

% For Recovery from inactivation of LCC
recoveryTime = 10;  % initialize to smallest value
%% Collect all parameters and define mass matrix for BAR module

p = [cycleLength,recoveryTime,CaMtotDyad,BtotDyad,CaMKIItotDyad,CaNtotDyad,PP1totDyad,...
    CaMtotSL,BtotSL,CaMKIItotSL,CaNtotSL,PP1totSL,...
    CaMtotCyt,BtotCyt,CaMKIItotCyt,CaNtotCyt,PP1totCyt...
    LCCtotDyad,RyRtot,PP1_dyad,PP2A_dyad,OA,PLBtot,LCCtotSL,PP1_SL...
    Ligtot,LCCtotBA,RyRtotBA,PLBtotBA,TnItotBA,IKstotBA,ICFTRtotBA,PP1_PLBtot,...
    PLMtotBA,MyototBA,IKrtotBA,IClCatotBA,CKIIOE];

% Need to define Mass matrix for BAR portion
m0 = ones(1,83+6+45+6);   % All state variables in other modules, not BAR
m1 =[0,     0,      0,      1,      1,      1,       1,       1,       1];
m2 =[0,     0,      0,      1,      1,      1];
m3 =[0,     0,      0];
m4 =[1,     1,      0,      0,      1,      1];
m5 =[1,     1,      0,      0,      1,      1,       1,       1      0,      0,      1,      1];

M = diag([m0,m1,m2,m3,m4,m5]); 
%% Establish and define globals

global tStep tArray I_Ca_store I_to_store I_Na_store I_K1_store ibar_store 
global gates_store Jserca_store IKs_store Jleak_store ICFTR_store Incx_store
global I_NaK_store I_kr_store I_Nabk_store
global Lmyo_store Fmyo_store
global IKp_store I_clca_store I_clbk_store Ipca_store Icabk_store Cai_store

tStep = 1; tArray = zeros(1,1e6); I_Ca_store=zeros(1,1e6); I_to_store=zeros(3,1e6);
I_Na_store = zeros(1,1e6); I_K1_store = zeros(1,1e6); ibar_store=zeros(1,1e6);
gates_store = zeros(2,1e6); Jserca_store = zeros(1,1e6); IKs_store = zeros(1,1e6);
Jleak_store = zeros(1e6,2); ICFTR_store = zeros(1,1e6); Incx_store = zeros(1,1e6);
I_NaK_store = zeros(1,1e6); I_kr_store = zeros(1,1e6); I_Nabk_store = zeros(1,1e6);

Fmyo_store = zeros(1,1e6); Lmyo_store = zeros(1,1e6);

IKp_store = zeros(1,1e6); I_clca_store=zeros(1,1e6); I_clbk_store=zeros(1,1e6);
Ipca_store=zeros(1,1e6); Icabk_store=zeros(1,1e6); Cai_store=zeros(1,1e6);
%% Assign initial conditions

% Isometric contraction (mechFlag = 0 in ECC module)
load yfinal_myo_WT_1Hz_isometric                   % control (steady-state, 1-Hz)
%load yfinal_myo_WT_1Hz_isometric_240iso            % 100 nM [ISO] (240-s, 1-Hz) 
%load yfinal_myo_WT_1Hz_isometric_240iso_XBCa       % ISO-XBCa
%load yfinal_myo_WT_1Hz_isometric_240iso_XBcy       % ISO-XBcy
%load yfinal_myo_WT_1Hz_isometric_240iso_XBCa_XBcy  % ISO-XBCa-XBcy
%load yfinal_myo_WT_1Hz_isometric_240iso_titin      % ISO-titin
%load yfinal_myo_WT_1Hz_isometric_240iso_cytofl     % ISO-cytofl

% Isotonic contraction (mechFlag = 1 in ECC module)
%load yfinal_myo_WT_1Hz_isotonic                    % control (steady-state, 1-Hz)
%load yfinal_myo_WT_1Hz_isotonic_240iso             % 100 nM [ISO] (240-s, 1-Hz) 
%load yfinal_myo_WT_1Hz_isotonic_240iso_XBCa        % ISO-XBCa
%load yfinal_myo_WT_1Hz_isotonic_240iso_XBcy        % ISO-XBcy
%load yfinal_myo_WT_1Hz_isotonic_240iso_XBCa_XBcy   % ISO-XBCa-XBcy
%load yfinal_myo_WT_1Hz_isotonic_240iso_titin       % ISO-titin
%load yfinal_myo_WT_1Hz_isotonic_240iso_cytofl      % ISO-cytofl

y0n=yfinal;
%% Run single simulation

tic
tspan = [0 3*1e3]; % [ms]
options = odeset('Mass',M,'RelTol',1e-5,'MaxStep',2); 
[t,y] = ode15s(@rabbit_myofilament_masterODEfile,tspan,y0n,options,p);
yfinal = y(end,:)';
toc

%save yfinal_myo_WT_1Hz_xxx yfinal
%% Rename outputs

tArray = tArray(1:tStep); I_Ca_store = I_Ca_store(1:tStep);
Ito = I_to_store(1,1:tStep); Itof = I_to_store(2,1:tStep);
Itos = I_to_store(3,1:tStep); INa = I_Na_store(1:tStep);
IK1 = I_K1_store(1:tStep); s1 = gates_store(1,1:tStep);
k1 = gates_store(2,1:tStep); Jserca = Jserca_store(1:tStep);
Iks = IKs_store(1:tStep); Jleak = Jleak_store(1:tStep,:);
ICFTR = ICFTR_store(1:tStep); Incx = Incx_store(1:tStep);
INaK = I_NaK_store(1:tStep); Ikr = I_kr_store(1:tStep);
INabk = I_Nabk_store(1:tStep);

Lm = Lmyo_store(1:tStep); Fm = Fmyo_store(1:tStep);

Ikp = IKp_store(1:tStep); Iclca = I_clca_store(1:tStep);
Iclbk = I_clbk_store(1:tStep); Ipca = Ipca_store(1:tStep);
Icabk = Icabk_store(1:tStep); Cai = Cai_store(1:tStep);
%% Plot

color = 'k';

figure, set(gcf,'color','w')
subplot(2,2,1); plot(t,y(:,39),color); hold on;
set(gca,'box','off','tickdir','out','fontsize',12); ylabel('Em (mV)')
subplot(2,2,2); plot(tArray,Fm,color); hold on;
set(gca,'box','off','tickdir','out','fontsize',12); ylabel('Force (mN/mm2)');
subplot(2,2,3); plot(t,y(:,38),color); hold on;
set(gca,'box','off','tickdir','out','fontsize',12); ylabel('[Ca]i (mM)')
xlabel('Time (ms)');
subplot(2,2,4); plot(tArray,Lm,color); hold on;
set(gca,'box','off','tickdir','out','fontsize',12); ylabel('Length (um)');
xlabel('Time (ms)');