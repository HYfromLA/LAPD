function [denoisedCCmatrix,newTh,cutoff] = denoise(n,d,I,J,W,sortedCCmatrix,th,denoisingmethod,SKNN)

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
        length(uniquenodes)
        sum(rowsTokeep) / length(knndistances)
       

    elseif strcmp(denoisingmethod, 'automatic_elbow')   
        c = randsample(sortedknndistances,floor(0.10*nn));  % downsample 10% of the knndistances.  
        sc = sort(c);
        idxx = knee_pt(sc);
        cutoff = sc(idxx)*1.0001; 

        rowsTokeep = knndistances < cutoff;
        uniquenodes = unique(sortedCCmatrix(rowsTokeep,m+1:end)); 

        while length(uniquenodes) < 0.95*n || sum(rowsTokeep) / length(knndistances)  < 0.85 %|| connection < 0.9
            
            cutoff = Ths(find(Ths>cutoff,1,'first')); 
            rowsTokeep = knndistances < cutoff; 
            uniquenodes = unique(sortedCCmatrix(rowsTokeep,m+1:end)); 

            %rowsTokeep1 = rowsTokeep(idx); denoisedGknn = Gknn(rowsTokeep1, rowsTokeep1);         
            %[~,binsizes] = conncomp(graph(denoisedGknn));
            %connection = max(binsizes)/sum(binsizes); 
        end

        close all;      
        figure
        scatter(1:nn, sortedknndistances, 'filled')
        xlim([0 nn]);
        title('LAPD to NN', 'FontSize',16)
        xlabel('Point Index', 'FontSize', 14)
        yline(cutoff,':','Color','red','LineWidth',2);
        text(2000, cutoff+0.05, ['cutoff = ' num2str(cutoff)], 'Color', 'red', 'FontSize', 15);
    
    elseif strcmp(denoisingmethod, 'automatic_connectedness') 
        cutoff = Ths(5);
        rowsTokeep = knndistances < cutoff;  denoisedCCmatrix = sortedCCmatrix(rowsTokeep, :);
        uniquenodes = unique(denoisedCCmatrix(:,m+1:end));
        connection = 0; 
       
         while  length(uniquenodes) < 0.9*n || connection < 0.9 || sum(rowsTokeep) / length(knndistances)  < 0.55
            cutoff = Ths(find(Ths>cutoff,1,'first'));  
            rowsTokeep = knndistances < cutoff;
            denoisedCCmatrix = sortedCCmatrix(rowsTokeep, :);
            uniquenodes = unique(denoisedCCmatrix(:,m+1:end)); 

            rowsTokeep1 = rowsTokeep(idx); denoisedGknn = Gknn(rowsTokeep1, rowsTokeep1);         
            [~,binsizes] = conncomp(graph(denoisedGknn));
            connection = max(binsizes)/sum(binsizes); 
         end

         rowsTokeep1 = rowsTokeep(idx); denoisedGknn = Gknn(rowsTokeep1, rowsTokeep1);       
         denoisedCCmatrix = sortedCCmatrix(rowsTokeep, :);

        close all; 
        figure
        scatter(1:nn, sortedknndistances, 'filled')
        xlim([0 nn]);
        title('LAPD to NN', 'FontSize',16)
        xlabel('Point Index', 'FontSize', 14)
        yline(cutoff,':','Color','red','LineWidth',2);
        text(2000, cutoff+0.05, ['cutoff = ' num2str(cutoff)], 'Color', 'red', 'FontSize', 15);

    end

    clear knndistances sortedknndistances
    Gknn = sparse(I, J, W, nn, nn); Gknn = max(Gknn, Gknn');
    rowsTokeep1 = rowsTokeep(idx); denoisedGknn = Gknn(rowsTokeep1, rowsTokeep1);   
    clear Gknn rowsTokeep1
    [I,J,W] = find(denoisedGknn); 
    clear denoisedGknn 
    denoisedsimplices = unique(sortedCCmatrix(rowsTokeep,(m+1):end), 'rows');
    clear sortedCCmatrix rowsTokeep
    
    [denoisedCCmatrix, newTh] = connectedcomponents(I,J,W,denoisedsimplices,m);
end