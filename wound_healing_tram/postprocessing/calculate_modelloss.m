clear; clc; close all;
set_plot_defaults
%{
1. Locate data
2. Import and format H5/XDMF files
3. Visualize and compare experimental data and forward solution
4. Quantitatively compare expt and fwd solution (calculate loss externally)
%}

well_num = 1:6;
% well_num = 1
well_let = {'A', 'B'};
well_let = {'A'};
% well_let = {'B'};
num_lbl = {'10% R1', '10% R2', '10% serum R3','3% R1', '3% R2', '3% serum R3'};
let_lbl = {'No MMC', '100 ng/ml MMC'};
let_lbl = {'No MMC'};
% let_lbl = {'100 ng/ml MMC'};
c=1;
for ii = 1:length(well_let)
    for jj = 1:length(well_num)
        wells{c} = sprintf('%s%02d', well_let{ii}, well_num(jj));
        L{c} = sprintf('%s - %s', let_lbl{ii}, num_lbl{jj});
        c = c+1;
    end
end
wells{1} = wells{2};
L{1} = L{2};

stepVsi = 5;
stepAdj = stepVsi;
%Replace these after rerunning VSI fwd
importStepsVSI = 5:7;%10:12;
importStepsAdj = importStepsVSI;
nTimes = 72;
for ii = 1:length(wells)
    well = wells{ii};
    for jj = 1:length(importStepsVSI)
        stepVsi = importStepsVSI(jj);
        stepAdj = importStepsAdj(jj);
        vsiFwdPath = 'C:\Users\pkinn\Dropbox (University of Michigan)\vsiglobus\ht_migration\results\forward_solution\VSI_2D_Time_Independent';
        vsiFwdLoc = sprintf('%s/well%s/step%i/',vsiFwdPath, well, stepVsi);
        [vsiDenseGrid(:,:,:,ii,jj), dim1Grid, dim2Grid] = importH5VSI('density.h5', vsiFwdLoc, nTimes);

        exptPath = 'C:\Users\pkinn\Dropbox (University of Michigan)\vsiglobus\ht_migration\results\PreProcess';
        exptLoc = sprintf('%s/%s/', exptPath, well);
        [exptDenseGrid(:,:,:,ii,jj), ~, ~] =  importH5VSI('density_3_3_rolling_win5.h5', exptLoc, nTimes);

        adjFwdPath = 'C:\Users\pkinn\Dropbox (University of Michigan)\vsiglobus\ht_migration\results\forward_solution\Adjoint_2D_Time_Independent';
        adjFwdLoc = sprintf('%s/well%s/step%i/',adjFwdPath, well, stepAdj);
        [adjDenseGrid(:,:,:,ii,jj), ~, ~] = importH5VSI('density.h5', adjFwdLoc, nTimes);
    end
end

%% Calculate time loss
nPts = length(dim1Grid(:));
for ii = 1:length(wells)
    well = wells{ii};
    for jj = 1:length(importStepsVSI)
        vsiTimeLoss(:,ii,jj) = (1/nPts) * sum((vsiDenseGrid(:,:,:,ii,jj) - exptDenseGrid(:,:,:,ii,jj)).^2, [1 2]);
        adjTimeLoss(:,ii,jj) = (1/nPts) * sum((adjDenseGrid(:,:,:,ii,jj) - exptDenseGrid(:,:,:,ii,jj)).^2, [1 2]);
        vsiTimeLoss(:,ii,jj) = sum((vsiDenseGrid(:,:,:,ii,jj) - exptDenseGrid(:,:,:,ii,jj)).^2, [1 2])./sum((exptDenseGrid(:,:,1,ii,jj)).^2, [1 2]);
        adjTimeLoss(:,ii,jj) = sum((adjDenseGrid(:,:,:,ii,jj) - exptDenseGrid(:,:,:,ii,jj)).^2, [1 2])./sum((exptDenseGrid(:,:,1,ii,jj)).^2, [1 2]);

    end
end
maxLoss = max([max(vsiTimeLoss,[], 1) max(adjTimeLoss, [], 1)]);

figure()
tl = tiledlayout(length(wells), length(importStepsVSI));

for ii = 1:length(wells)
    well = wells{ii};
    for jj = 1:length(importStepsVSI)
        step = importStepsAdj(jj);
        nexttile()
        plot(vsiTimeLoss(:,ii,jj), 'r')
        hold on
        plot(adjTimeLoss(:,ii,jj), 'b')
        title(sprintf('Well: %s step: %i', well, step))
        ylim([0 maxLoss(jj)])
    end
end
xlabel(tl, 'Time steps', 'FontSize', 20)
ylabel(tl, 'Mean Square Error', 'FontSize', 20)

colors = {'r', 'g', 'b'};
legText = {};
figure()
tl = tiledlayout(2, 3);
for ii = 1:length(wells)
    well = wells{ii};
    nexttile()
    for jj = 1:length(importStepsAdj)

        step = importStepsAdj(jj);
        plot(vsiTimeLoss(:,ii,jj), 'Color',colors{jj}, 'LineStyle','-')
        hold on
        plot(adjTimeLoss(:,ii,jj), 'Color',colors{jj}, 'LineStyle',':')
        if ii == 1
            legText = [legText sprintf('VSI step %i', step) sprintf('Adj step %i', step)];
        end
    end
    title(sprintf('Well: %s', well))
end
legend(legText)
xlabel(tl, 'Time steps', 'FontSize', 20)
ylabel(tl, 'Mean Square Error', 'FontSize', 20)


figure()
ii = 1
well = wells{ii};
nexttile()
for jj = 1:length(importStepsAdj)

    step = importStepsAdj(jj);
    plot(vsiTimeLoss(:,ii,jj), 'Color',colors{jj}, 'LineStyle','-')
    hold on
    plot(adjTimeLoss(:,ii,jj), 'Color',colors{jj}, 'LineStyle',':')
    if ii == 1
        legText = [legText sprintf('VSI step %i', step) sprintf('Adj step %i', step)];
    end
end
title(sprintf('Well: %s', well))

legend(legText)
xlabel('Time steps', 'FontSize', 20)
ylabel('Mean Square Error', 'FontSize', 20)


