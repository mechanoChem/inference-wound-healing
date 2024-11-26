clear; clc; close all;
set_plot_defaults_export

basisNamesWF = {'$-\nabla W \cdot \nabla C$','$-C*\nabla W \cdot \nabla C$', '$-C^{2}*\nabla W \cdot \nabla C$', ...
    '$C*\nabla W \cdot v$', '$C \nabla W \cdot v * C$', '$C* \nabla W \cdot v*C^{2}$', ...
    '$C*W$', '$C^{2}*W$'};
% basisNames = {'$-\nabla^2 C$','$-C*\nabla^2 C$', '$-C^{2}*\nabla^2 C$', ...
%     '$-v \cdot \nabla C$', '$-C* v \cdot \nabla C$', '$-C^2* v \cdot \nabla C$', ...
%     '$C$', '$C^{2}$', 'tracking'};
basisNames =  {'$-\nabla^2 C$','$-C*\nabla^2 C$', '$-C^{2}*\nabla^2 C$', ...
    '$-v \cdot \nabla C$', '$-C* v \cdot \nabla C$', '$-C^2* v \cdot \nabla C$', ...
    '$C$', '$C^{2}$'};
figLoc = 'figures/';
concVec = 10000:2000:20000;
concVecLbl = {'initCond10000', 'initCond12000', 'initCond14000', 'initCond16000', 'initCond18000', 'initCond20000'};
timeLabels = {'0 hours Expt.', '0 hours Model', '12 hours', '24 hours', '36 hours', '48 hours'};
timeLabels2 = {'0 hours', '12 hours', '24 hours', '36 hours', '48 hours'};
nTimes = 5;
%% Plot comparison of VSI and Adjoint
% Want to visualize how the model changes at each step
for ii = 1:length(concVec)
    initConc = concVec(ii);
    % adjFwdDir = ['\\wsl$\Ubuntu-20.04\home\pkinn\vsiTestWSL\Cell_Migration_density_rerun10_17_22\results\forward_solution\Adjoint_1D_Time_Independant\initCond', num2str(initConc)];
    % vsiFwdDir = ['\\wsl$\Ubuntu-20.04\home\pkinn\vsiTestWSL\Cell_Migration_density_rerun10_17_22\results\forward_solution\VSI_1D_Time_Independant\initCond', num2str(initConc)];
    % jinFwdDir = ['\\wsl$\Ubuntu-20.04\home\pkinn\vsiTestWSL\Cell_Migration_density_rerun10_17_22\results\forward_solution\UserDefined\initCond', num2str(initConc), '\step0\'];
    adjFwdDir = ['../results/forward_solution/Adjoint_1D_Time_Independant/initCond', num2str(initConc)];
    vsiFwdDir = ['../results/forward_solution/VSI_1D_Time_Independant/initCond', num2str(initConc)];
    jinFwdDir = ['../results/forward_solution/UserDefined/initCond', num2str(initConc), '/step0/'];    
    toImport = 6;

    % adjName = sprintf('%s\\step%i\\%s',adjFwdDir, toImport, 'density.h5');
    adjName = sprintf('%s/step%i/%s',adjFwdDir, toImport, 'density.h5');
    [adjDense, adjMesh] = importFenicsModelDensity1D(adjName, {'0', '1', '2','3'});
    adjDenseAll(:,:,ii) = adjDense;
%     vsiName = sprintf('%s\\step%i\\%s',vsiFwdDir, toImport, 'density.h5');
%     [vsiDense, vsiMesh] = importFenicsModelDensity1D(vsiName, {'0', '1', '2','3'});
%     vsiDenseAll(:,:,ii) = vsiDense;
    jinName = sprintf('%s%s', jinFwdDir, 'density.h5');
    [jinDense, jinMesh] = importFenicsModelDensity1D(jinName, {'0', '1', '2', '3'});
    jinDenseAll(:,:,ii) = jinDense

    %exptLoc = strcat('../results/PreProcess/density', num2str(initConc), '/density_1D_3_3_rolling_win3_refine4.h5');
    exptLoc = strcat('../results/PreProcess/density', num2str(initConc), '/density_1D_3_3_rolling_win1_refine4.h5');
    [exptDensity(:,:,ii), exptMesh] = importFenicsModelDensity1D(exptLoc, {'0', '1','2', '3','4'});
end

%% Calculate MSE
%each density is nTimes x nPoints x nConcentrations. Want to calculate an
%MSE for each concentration.
[nTimeCompare, nSpaceCompare, nConc] = size(adjDenseAll);
for ii = 1:nConc
    jinDif = jinDenseAll(:,:,ii) - exptDensity(2:end, :, ii);
    jinMSE(ii) = 1/(nTimeCompare*nSpaceCompare)*sum(jinDif(:).^2);
    jinRMSE(ii) = sqrt(jinMSE(ii));
    
    adjDif = adjDenseAll(:,:,ii) - exptDensity(2:end, :, ii);
    adjMSE(ii) = 1/(nTimeCompare*nSpaceCompare)*sum(adjDif(:).^2);
    adjRMSE(ii) = sqrt(adjMSE(ii));
end

adjRMSE


%Plot MSEs
figure()
barLbl = {'Seed: 10000', 'Seed: 12000', 'Seed: 14000', 'Seed: 16000', 'Seed: 18000', 'Seed: 20000'};
x = categorical(barLbl);
x = reordercats(x, barLbl)
yBar = [jinRMSE; adjRMSE]';
bar(x, yBar)
legend('Jin et. al.', 'VSI+Optim', 'Location', 'northwest')
ylabel('RMSE (cells/\mum^2)')
exportgraphics(gcf, [figLoc 'jinAdjRMSEcompare.pdf'], 'ContentType','vector');
exportgraphics(gcf, [figLoc 'jinAdjRMSEcompare.jpg'], 'Resolution', 1500);





