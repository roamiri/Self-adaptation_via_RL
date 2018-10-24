%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Drawing reward function
%
% function FBS_out = draw_R( FBS_in, MBS, mue, Npower, fbsCount, femtocellPermutation, NumRealization, saveNum, CL)

%% Initialization
% clear all;
clc;
% format short
% format compact
%% Generate the UEs
mue = UE(25, 365);
MBS = BaseStation(0 , 0 , 33);
%Generate fbsCount=16 FBSs, FemtoStation is the agent of RL algorithm
FBS_Max = cell(1,10);
for i=1:5
    FBS_Max{i} = FemtoStation_dual_strip((i-1)*10+5,350, MBS, mue, 5, 0);
end

for i=1:5
    FBS_Max{i+5} = FemtoStation_dual_strip((i-1)*10+5,380, MBS, mue, 5, 1);
end
%% Parameters
Pmin = 5;                                                                                                                                                                                                                                                                                                                                                                           %dBm
Pmax = 15; %dBm
%% Minimum Rate Requirements for N MUE users
N = 3;
q_mue = 4.0; q_fue=0.50;
fbsCount = 10;
FBS = FBS_Max(1:fbsCount);
%% Initialize the new Agent (FBS)
%% Calc channel coefficients
fbsNum = size(FBS,2);
G = zeros(fbsNum+1, fbsNum+1); % Matrix Containing small scale fading coefficients
L = zeros(fbsNum+1, fbsNum+1); % Matrix Containing large scale fading coefficients
[G, L] = measure_channel_3GPP(FBS,MBS,mue,100);
%% Main Loop
% calc FUEs and MUEs capacity
% MUE minimum SINR
fbsNum = size(FBS,2);
for j=1:fbsNum
    fbs = FBS{j};
    fbs.P = Pmax;
    FBS{j} = fbs;
end
mue.SINR_min = SINR_MUE_4(G, L, FBS, MBS, mue, -174);
% MUE maximum SINR
for j=1:fbsNum
    fbs = FBS{j};
    fbs.P = Pmin;
    FBS{j} = fbs;
end
mue.SINR_max = SINR_MUE_4(G, L, FBS, MBS, mue, -174);
% Last FUE minimum SINR
for j=1:fbsNum-1
    fbs = FBS{j};
    fbs.P = Pmax;
    FBS{j} = fbs;
end
FBS{fbsNum}.P = Pmin;
SINR_FUE_Vec = SINR_FUE_2(G, L, FBS, MBS, -174);
FBS{fbsNum}.SINR_min = log2(1+SINR_FUE_Vec(fbsNum));
% Last FUE maximum SINR
for j=1:fbsNum-1
    fbs = FBS{j};
    fbs.P = Pmin;
    FBS{j} = fbs;
end
FBS{fbsNum}.P = Pmax;
SINR_FUE_Vec = SINR_FUE_2(G, L, FBS, MBS, -174);
FBS{fbsNum}.SINR_max = log2(1+SINR_FUE_Vec(fbsNum));

x = linspace(mue.SINR_min, mue.SINR_max, 100);
y = linspace(FBS{fbsNum}.SINR_min, FBS{fbsNum}.SINR_max, 100);

% r0 = log2(1+x);
% ri = log2(1+y);

r0 = linspace(0,10,100);
ri = linspace(0,10,100);

[xx,yy] = meshgrid(r0,ri);
beta = FBS{fbsNum}.dMUE/18;

% zz = beta.*yy.*xx.^2 -(yy-q_fue).^2 - (1/beta)*(xx-q_mue).^2;
% zz = -(yy-q_fue).^2 - (xx-q_mue).^2;
% zz = exp(-(xx-q_mue).^2)-exp(-yy);
% zz = (yy-q_fue).^3 + (xx-q_mue).^3;
zz = (beta).^sign(xx-q_mue).*(yy-q_fue).^3 + (beta).^sign(xx-q_mue).*(xx-q_mue).^3;
% zz = yy - (1./beta).*(xx-q_mue).^2;
% zz = beta.*yy - (1/beta)*(xx-q_mue).^2;
zz = zz / max(max(abs(zz)));

figure;
% mesh(zz);
surface(xx,yy,zz);
xlabel('r_{0}','FontSize',12);%, 'FontWeight','bold');
ylabel('r_{k}','FontSize',12);%, 'FontWeight','bold');
box on;
grid on;
% xlim([mue.SINR_min, mue.SINR_max]);
colorbar;
% [px,py] = gradient(zz,.2,.2);
% quiver(x,y,px,py);

% end
