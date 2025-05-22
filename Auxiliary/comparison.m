%% Plots for comparison

function [d,epsilon,k_hat,labels,time,misc] = comparison(X, DATAopts,labels,labelsGT)

viz.markersize = 15;
figure; 
label_ids = unique(labelsGT);
if (DATAopts.intrdim == 1 && strcmp(DATAopts.shape, "three curves")) || DATAopts.intrdim == 2
    subplot(1,2,1)
    hold on;
    for label = label_ids(:).'
        mask = (labelsGT == label);
        plot3(X(mask,1), X(mask,2), X(mask,3), '.', viz);
    end
    title('Oracle labels')
    hold off

    subplot(1,2,2)
    hold on;
    for label = label_ids(:).'
        mask = (labels == label);
        plot3(X(mask,1), X(mask,2), X(mask,3),'.', viz);
    end
    title('LAPD labels')
    hold off

elseif DATAopts.intrdim == 1
    subplot(1,2,1)
    hold on;
    for label = label_ids(:).'
        mask = (labelsGT == label);
        plot(X(mask,1), X(mask,2), '.', viz);
    end
    title('Oracle labels')
    hold off

    subplot(1,2,2)
    hold on;
    for label = label_ids(:).'
        mask = (labels == label);
        plot(X(mask,1), X(mask,2), '.', viz);
    end
    title('LAPD labels')
    hold off
end
end

