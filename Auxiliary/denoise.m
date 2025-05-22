function [denoisedsimplices,rowsTokeep,cutoff] = denoise(n,d,sortedCCmatrix,th,denoisingmethod,SKNN)

% This function removes simplices whose LAPD to its SKNN-th nearest neighbors are large (beyond cutoff). 
% Inputs:  d: intrinsic dimension of manifolds. 
%          n: dataset size. 
%          I: row indices of nonzero entries in the adjacency matrix. 
%          J: column indices of nonzero entries in the adjacency matrix. 
%          W: weights of nonzero entries in the adjacency matrix. 
%          sortedCCmatrix: the original CCmatrix. 
%          th: original scales coming with the sortedCCmatrix. 
%          denoisingmethod: "examinefigure", "automatic_elbow", or "automatic_connectedness". 
%          parallel: parallel computing or not. 
% Outputs: denoisedCCmatrix: the denoised CCmatrix. 
%          newTh: new scales coming with the denoisedCCmatrix. 
%          SKNN: simplices' nearest neighbor to look at for denoising. 
%          cutoff: the threshold to remove simplices with large knn distances. 

[nn,mm] = size(sortedCCmatrix); m = mm-(d+1);

knndistances = zeros(nn, 1); 
    for s=1:nn
        iup = s; %will tell us the index of a nearest neighbor when it is found
        idown = s; %will tell us the index of a nearest neighbor when it is found
        counter = 2; %keeps track of how many nearest neighbors you have found

        for j=1:m %defines the approximate path distance
            % Move up until you find a point not in the same CC then move
            % left and repeat
            while (iup>1 && sortedCCmatrix(iup,j)==sortedCCmatrix(iup-1,j) && counter<SKNN+1)
                iup = iup-1;
                if counter==SKNN
                    knndistances(s)=th(j+1);
                end
                counter=counter+1;
            end
            % Move down until you find a point not in the same CC then move
            % left and repeat
            while (idown<nn && sortedCCmatrix(idown,j)==sortedCCmatrix(idown+1,j) && counter<SKNN+1)
                idown = idown+1;
                if counter==SKNN
                    knndistances(s)=th(j+1);
                end
                counter=counter+1;
            end
        end
        % If we have not found k2 NN b/c we have disconnected components at
        % the largest scale, simply move up/down to add additional
        % neighbors at infinite distance
        while (iup>1 && counter<SKNN+1)
            iup = iup-1;
            if counter==SKNN
                knndistances(s)=Inf;
            end
            counter=counter+1;
        end
        while (idown<nn && counter<SKNN+1)
            idown = idown+1;
            if counter==SKNN
                knndistances(s)=Inf;
            end
            counter=counter+1;
        end
    end
    
    [~, idx] = unique(sortedCCmatrix(:,m+1:end), "rows");
    uniqueDDs = unique(knndistances); Ths = 1.0001*uniqueDDs; 
    sortedknndistances = sort(knndistances); 
    
    if strcmp(denoisingmethod, 'examinefigure')
        close all; 
        figure
        scatter(1:nn, sort(knndistances), 'filled')
        xlim([0 nn]);
        title('LAPD to NN', 'FontSize',16)
        xlabel('Point Index', 'FontSize', 14)
        prompt='At what Knn-LAPD value to threshold for denoising? \n';
        
        cutoff = input(prompt); rowsTokeep = knndistances < cutoff;
        
        uniquenodes = unique(sortedCCmatrix(rowsTokeep,m+1:end)); 
        fprintf('Percent of nodes surviving denoising %d \n', length(uniquenodes) / n);
        fprintf('Percent of simplices surviving denoising %d \n', sum(rowsTokeep) / length(knndistances));
       
    elseif strcmp(denoisingmethod, 'automatic')

        c = randsample(sortedknndistances(floor(nn*0.25):end),floor(0.10*nn));  % downsample 10% of the knndistances.  
        sc = sort(c);
        idxx = knee_pt(sc);
        cutoff = sc(idxx)*1.01;
        
        rowsTokeep = knndistances < cutoff;
        uniquenodes = unique(sortedCCmatrix(rowsTokeep,m+1:end)); 
        valid_num_knns = sum(knndistances < 5);
        if valid_num_knns/length(knndistances) < 0.85
            warning("%.2f%% of simplices has infinity distance to their %d-th neighbor; Graph connection is bad. \n Please relax the filter and re-run.", valid_num_knns/length(knndistances)*100, SKNN);
        end

        while (length(uniquenodes) / n) < 0.95 || sum(rowsTokeep) / valid_num_knns  < 0.85 % valid_num_knns
            
            cutoff = Ths(find(Ths>cutoff,1,'first')); 
            rowsTokeep = knndistances < cutoff; 
            uniquenodes = unique(sortedCCmatrix(rowsTokeep,m+1:end)); 
        end
    else
        error("Unknown denoising method.")
    end

    sortedknndistances = sortedknndistances(sortedknndistances < 1.8);
    close all;      
    figure
    scatter(1:length(sortedknndistances), sortedknndistances, 'filled')
    xlim([0 length(sortedknndistances)]);
    title('LAPD to NN', 'FontSize',16)
    xlabel('Point Index', 'FontSize', 14)
    yline(cutoff,':','Color','red','LineWidth',2);
    text(2000, cutoff+0.05, ['cutoff = ' num2str(cutoff)], 'Color', 'red', 'FontSize', 15);
    
    clear knndistances sortedknndistances
    denoisedsimplices = unique(sortedCCmatrix(rowsTokeep,(m+1):end), 'rows');
    clear sortedCCmatrix
    rowsTokeep = rowsTokeep(idx);
end