% find the Floque mulitipliers of the 8-dimension airfoil model 
% INPUT:
% x_p: the state variable on the chosen poincare section
% T: the evolution time for the orbit
% fs: the sampling frequency of the numerical integration
% F_ls: the tangent vector of the free flight 
% RMap : reset map when hitting the discontinuity surface
% EventFun: event function to detect the impact on the surface
% 
%
% Get the Lyapunov exponent of the Composed map 
function LYE_airfoil_Model
addpath ../public_scripts
addpath ../Postprocessing Scripts
% 初始化
clc
close all
clear

% 导入模型数据 得到基本矩阵
global MHH res gap
Stiffness_ratio=1;
% start computation
gap=0.01;
heave0=0;
alpha0=0.06;
beta0=1e-4/Stiffness_ratio;
damp0=1;
%
Stiff_k1=0;
Stiff_k2=300*sqrt(Stiffness_ratio);

method='Pade';
FLUTEER_MATRIXBASE=Initialization_backlash(gap,heave0,alpha0,beta0,damp0,Stiff_k1,Stiff_k2);
% 定义间隙大小

gap=FLUTEER_MATRIXBASE.TimeInt.gap;
MHH=FLUTEER_MATRIXBASE.MHH;
fs=FLUTEER_MATRIXBASE.TimeInt.fs;
V_care=0.64833*30;
[A_1,A_2,M_new,K_new,C_new,vect_preload]=Get_int_Matrix(FLUTEER_MATRIXBASE,method,V_care);
res=0.72;

F_ls =@(t,y) A_1*y+vect_preload(:,1);
%%


% the start pint on the poincare section
x_p = [-0.00465412709916764,0.0606250812329585,-0.00962039197420514,-0.00102284858298862,-0.000663458443753684,0.00466523614885600,0.0624456248892721,-5.96091134697107e-06]';

%define the period T (found previously by numerical)
T=0.132671235175032;
PlotInd = [3 3 6 1 2];
[MainV,LPE]=ComposedMap_FloqueMPs(x_p,T,fs,F_ls,@RMap,@DetectingContact_events,PlotInd);
keyboard


function y0=RMap(y0)
global MHH res
delta_v=MHH(1:2,1:2)\[MHH(1,3);MHH(2,3)]*(1+res)*y0(6);
y0(6)=-res.*y0(6);
y0(4)=y0(4)+delta_v(1);
y0(5)=y0(5)+delta_v(2);        

%
function [value,isterminal,direction] = DetectingContact_events(t,y)
% Locate the time when flap hit through boundary in both directions
% and stop integration when isterminal ==1
% Events: switch boundary; zero heave velocity; zero pitch velocity; zero
% flap velocity
global gap
value = real([gap^2-y(3)^2;abs(y(6))]);                %  detect switching point
isterminal = [1;0];                              %  stop the integration
direction  = [-1;0]; 