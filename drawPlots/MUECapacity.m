clear;
clc;

%%
sum_FUE_IL = [];
C_MUE_IL = [];
C_FUE_Mat_IL = cell(1,16);
P_FUE_Mat_IL = cell(1,16);
for i=1:16
    fprintf('FBS num = %d\t', i);
    sumfue = 0.;
    c_mue = 0.0;
    c_fue_vec = zeros(1,i);
    p_fue_vec = zeros(1,i);
    Cnt = 0;
    
    for j=1:100
        s = sprintf('Jun13/R_1/pro_%d_%d.mat',i,j);
        filename = strcat(s);
        if exist(s)
            load(filename);
            c_mue = c_mue + QFinal.mue.C;
            sumfue = sumfue + QFinal.sum_CFUE;
            c_fue_vec = c_fue_vec + QFinal.C_FUE;
            pp = zeros(1,i);
            for kk = 1:i
                pp(1,kk) = QFinal.FBS{kk}.P;
            end
            p_fue_vec = p_fue_vec + pp;
            Cnt = Cnt+1;
        end
    end
    fprintf('Total Cnt = %d\n',Cnt);
    C_MUE_IL(i) = c_mue./Cnt;
    sum_FUE_IL = [sum_FUE_IL sumfue/Cnt];
    C_FUE_Mat_IL{i} = c_fue_vec./Cnt;
    P_FUE_Mat_IL{i} = p_fue_vec./Cnt;
end

%%
dirName = 'nopunish';
listing=dir(dirName);


mine = [];    
for i=1:16
    s = sprintf('R_nopunish_mix_L:3,%d,Real:1000',i);
    filename = strcat(dirName , '/', s);
    load(filename);
    C = QFinal.mue.C_profile;
    cc = sum(C(40000:size(C,2)))/(-40000+size(C,2));
    mine = [mine cc];
end

%%
figure;
hold on;
grid on;
box on;
plot( ones(1,14)*1.0, '--b', 'LineWidth',1);
for i=1:14
    vec = C_FUE_Mat_IL{i};
    for j=1:size(vec,2)
        plot(i,vec(j), 'sr', 'LineWidth',1.5,'MarkerSize',10, 'MarkerFaceColor','r', 'MarkerEdgeColor','b');
    end
end
% plot(min_FUE, '--r', 'LineWidth',1,'MarkerSize',10);
% plot(min_FUE_ref, '--b', 'LineWidth',1,'MarkerSize',10);
% title('Capacity of all members of a cluster','FontSize',14, 'FontWeight','bold');
xlabel('Cluster size','FontSize',14, 'FontWeight','bold');
ylabel('Capacity(b/s/Hz)','FontSize',14, 'FontWeight','bold');
% ylabel('Power(dBm)','FontSize',14, 'FontWeight','bold');
% xlim([2 14]);
% ylim([0 20]);
% legend({'\textbf{required QoS ($log_2(q_k)$)}','\textbf{CDP-Q/IL}'},'FontSize',14, 'Interpreter','latex');
%%
figure;
hold on;
grid on;
box on;
plot( ones(1,14)*1.0, '--b', 'LineWidth',1);
plot(C_MUE_IL, '--r', 'LineWidth',1,'MarkerSize',10);
% plot(min_FUE_ref, '--b', 'LineWidth',1,'MarkerSize',10);
% title('Capacity of all members of a cluster','FontSize',14, 'FontWeight','bold');
xlabel('Cluster size','FontSize',14, 'FontWeight','bold');
ylabel('Capacity(b/s/Hz)','FontSize',14, 'FontWeight','bold');
% ylabel('Power(dBm)','FontSize',14, 'FontWeight','bold');
% xlim([2 14]);
% ylim([0 20]);
% legend({'\textbf{required QoS ($log_2(q_k)$)}','\textbf{CDP-Q/IL}'},'FontSize',14, 'Interpreter','latex');
%%
% dirName = 'nopunish_3';
% listing=dir(dirName);
% 
% 
% mine2 = [];    
% for i=1:16
%     s = sprintf('R_nopunish_2:%d,Real:5000',i);
%     filename = strcat(dirName , '/', s);
%     load(filename);
%     C = QFinal.mue.C_profile;
%     cc = sum(C(40000:size(C,2)))/(-40000+size(C,2));
%     mine2 = [mine2 cc];
% end
%%
figure;
hold on;
grid on;
box on;
plot( ones(1,16)*1.0, '--k', 'LineWidth',1 );
% plot(mine, '--*r', 'LineWidth',1,'MarkerSize',10);
plot(R1, '--sb','LineWidth',1,'MarkerSize',10);
% plot(mine2, '--*g');
title('MUE capacity in low interference','FontSize',14, 'FontWeight','bold');
xlabel('FBS Numbers','FontSize',14, 'FontWeight','bold');
ylabel('Capacity(b/s/HZ)','FontSize',14, 'FontWeight','bold');
legend({'threshold','RF1','RF2'},'FontSize',14, 'FontWeight','bold');
% axis([0 16 0 5])
% saveas(gcf,sprintf('FUE_Number_%d.jpg'))
