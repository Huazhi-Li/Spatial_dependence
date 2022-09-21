%%% Evaluate Kmeans clustering solutions
%%% Calculating the optimal number of clusters based on the ratio of mean correlation and mean standard deviation
%%% Author: Huazhi Li, huazhi.li@vu.nl
%%% Date: 20-07-2021

clear; close all; clc 

% general folder (where the present script is located)
gen_folder= [pwd '\'];
addpath(genpath([pwd 'functions']))

% input data
data_folder= 'ESL_time_frames_and_zones\Updated\';

% output folder
foldsave0 = [gen_folder 'KMeans_evaluation\Updated\'];

%% number of clusters and repetitions

NC= 5:30;
Replicates = 50;

%% Coastal stations zoning

zones = {'IO', 'NA', 'NEA', 'NWA', 'NWP', 'Oceania', 'SA', 'SEA', 'SWA', 'SWP'}; % Global zones

for zz = 3   
    
    zone = zones{zz};    
    NC_opt = zeros(1,3);
    
    dir_data= dir([gen_folder data_folder zone '.mat']);

    for ii= 1: size(dir_data,1)
        
        % load data
        load([gen_folder data_folder dir_data(ii).name])

        name_zone= dir_data(ii).name;
        name_zone= name_zone(1:strfind(name_zone,'.')-1);

        waterlevel = eval([name_zone '.waterlevel']);
        lonlat = eval([name_zone '.lonlat']);
        time  = eval([name_zone '.time']);

        %% Extreme surge levels (esl)

        % Choose quantile
        q= .99;

        % thresholds
        ths= quantile(waterlevel,q);

        % nan to those events below the threshold
        esl= waterlevel; % rows= time, columns= location

        for j= 1: size(lonlat,1)

            esl(esl(:,j)< ths(j),j)= nan;

        end

        % declustering (3-day window)

        for j= 1: size(lonlat,1)

            dec_esl= nan(size(esl(:,j)));

            while sum(isnan(esl(:,j))) ~= length(esl(:,j))
                [vmax,fmax]= max(esl(:,j));
                dec_esl(fmax)= vmax;

                % to don't pick that event again
                if fmax== 1
                    esl(1:fmax+1,j)= nan;
                elseif fmax== length(esl(:,j))
                    esl(fmax-1:fmax,j)= nan;
                else
                    esl(fmax-1:fmax+1,j)= nan;
                end

            end

            esl(:,j)= dec_esl;

        end


        % Kmeans doesn't allow nans
        esl = esl';
        esl(isnan(esl))= 0;
        
%         %% Evaluate the clustering using CalinskiHarabasz and DaviesBouldin Criterion
%         
%         rng('default');
% 
%         % to be able to use 'distance' option in evalclusters function
%         myfunc = @(X,K)(kmeans(X, K,'Distance','correlation','Replicate',50));
% 
%         eva = evalclusters(esl,myfunc,'CalinskiHarabasz','klist',NC);    
%         NC_opt(1) = eva.OptimalK;
%         
%         eva = evalclusters(esl,myfunc,'DaviesBouldin','klist',NC);
%         NC_opt(2) = eva.OptimalK;


        %% Evaluate the clustering using HandMade Criterion
        foldsave= [foldsave0 zone];
        mkdir(foldsave)

        rng('default');
        opts = statset('Display','off');
        
        % correlations depending on the number of clusters
        correlations_clusters= cell(1,length(NC));
        c= 0;

        for num_clusters= NC
            [idx,C,sumd,D]= kmeans(esl,num_clusters,'Distance','correlation',...
                'Replicates',50,'Options',opts); % ,'OnLinePhase','on');
            
            % Make cluster orders
            clusters= unique(idx);

            corrs_centroides= nan(num_clusters,1);
            
            for i= 1: num_clusters

                % To calculate the correlations inside of each cluster
                centroide= esl(idx== clusters(i),:);
                centroide_cell= cell(size(centroide,1),1);

                for ii= 1: size(centroide,1)
                    ci= centroide(ii,:)';
                    ci(ci==0)= nan;
                    centroide_cell{ii}= ci;
                end

                centroide_cell2= repmat(centroide_cell,1,length(centroide_cell));
                R_centroide    = cellfun(@(x,y)corr(x,y,'rows','complete'),...
                    centroide_cell2,centroide_cell2');

                %     Diagonal matrix to calculate the mean correlations without repetitions...
                % .........................................................................
                X= nan(size(R_centroide));

                for ix= 1: length(X)

                    ix2= ix+1;
                    X(ix,ix2:end)= ones(1,length(X(ix,ix2:end)));
                end
                % .........................................................................

                corrs_centroides(i)= nanmean(R_centroide(X==1));
            end

            c= c+1;
            correlations_clusters{c}= corrs_centroides;

            disp([num2str(num_clusters) ' done!']);

        end

        %% Stats mean corr / mean std
        stats_centroids = nan(size(correlations_clusters,2),3);
        
        for irep= 1: size(correlations_clusters,2)

            num_centroid_i= correlations_clusters{irep};

            mean_centroid= mean(num_centroid_i(:));

            std_centroid = std(num_centroid_i(:));
            
            ratio_centroid = mean_centroid/std_centroid;

            stats_centroids(irep,:)= [mean_centroid, std_centroid ratio_centroid];
        end

        NC_opt(3) = NC(stats_centroids(:,3)==max(stats_centroids(:,3)));
        
        save([foldsave '\Stats.mat'],'NC_opt','stats_centroids','-mat');
        
        
        hh= figure;
        set(hh,'units','centimeters','Position',[1 1.8 30 15],'InvertHardCopy','off',...
            'resize','off','PaperPositionMode','auto','PaperType','A0','visible','on','color','w');

        errorbar(stats_centroids(:,1),stats_centroids(:,2),'.','MarkerSize',20,'LineWidth',1.5)
        set(gca,'XTick',1:length(stats_centroids),'XTicklabel',(NC))
        xlabel('Number of clusters')
        ylabel('Mean correlation')
        title(['Ratio of mean correlation to standard deviation, ' zone])
        
        
        ax= gca;
        ax.YAxis.FontSize= 16;
        ax.XAxis.FontSize= 16;
        
        
        print([foldsave '\HandMade criteria STATS'], '-dpng','-r0')
        

        close all
            
    end
end

