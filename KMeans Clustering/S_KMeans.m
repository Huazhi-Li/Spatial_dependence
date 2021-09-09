%%% KMeans clustering
%%% Using Kmeans clustering method to group global extreme sea levels into subgroups based on their statistical correlation 
%%% Over 10 pre-defined coastal zones
%%% Author: Huazhi Li, huazhi.li@vu.nl
%%% Date: 20-07-202

clear; close all; warning off; clc 

% general folder (where the present script is located)
gen_folder= [pwd '\'];

addpath(genpath([pwd 'functions']))

% input data
data_folder= 'ESL_time_frames_and_zones\Updated\';

% folder to save
fsave0 = [gen_folder 'ESL_KMeans\Updated\'];


%% number of clusters and repetitions

NC= [5:30];
Replicates  = 50;

%% Coastal stations zoning

zones = {'IO', 'NA', 'NEA', 'NWA', 'NWP', 'Oceania', 'SA', 'SEA', 'SWA', 'SWP'}; % Global zones

for zz = 7:10
    
    zone = zones{zz}; 
    
    fsave = [fsave0 zone '\'];
    mkdir(fsave)

    dir_data= dir([gen_folder data_folder zone '.mat']);

    for ii= 1: size(dir_data,1)

        load([gen_folder data_folder dir_data(ii).name])

        name_zone= dir_data(ii).name;
        name_zone= name_zone(1:strfind(name_zone,'.')-1);

        waterlevel = eval([name_zone '.waterlevel']);
        lonlat = eval([name_zone '.lonlat']);
        time  = eval([name_zone '.time']);

        %% Extreme sea levels (esl)

        % Choose quantile
        q= .95;

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


        %% Kmeans doesn't allow nans
        esl = esl';
        esl(isnan(esl))= 0;


        %% number of clusters and repetitions
        for iii= 1 : length(NC)

            num_clusters= NC(iii);

            %% Run KMeans
            [KMeans_res]= kmeans_fun(esl,num_clusters,Replicates,lonlat);

            lonlat= KMeans_res{2,1}(:,1:2);
            idx   = KMeans_res{2,1}(:,3);

            ref_series= KMeans_res{2,5};

            %% Plot
            
            % limits of the area
            LonMin= min(lonlat(:,1)); LonMax= max(lonlat(:,1));
            LatMin= min(lonlat(:,2)); LatMax= max(lonlat(:,2));

            % To save the idx's combination
            rgb= pmkmp(num_clusters,'jet');

            hh= figure;
            set(hh,'units','centimeters','Position',[-27.6860    5.8843   20.0025   12.0015],...
                'InvertHardCopy','off',...
                'resize','off','PaperPositionMode','auto','PaperType','A0','visible','on',...
                'color','w');

            m_proj('miller','long',[LonMin-1 LonMax+1],'lat',[LatMin-1 LatMax+1]);
            m_coast('patch',[0.8594    0.8594    0.8594]);
            m_grid('tickdir','in','linest','none','FontName','Times','FontSize',12)

            for i= 1: num_clusters

                hold all;
                m_line(lonlat(idx== ref_series(i,3),1),...
                    lonlat(idx== ref_series(i,3),2),'Color',rgb(i,:)...
                    ,'Marker','.','Linest','none','markersize',12);

                % centroid
                hold on;
                m_line(ref_series(i,1),ref_series(i,2),'Color',rgb(i,:),...
                    'Marker','o','Linest','none','LineWidth',2,'MarkerSize',20)

            end

            ht= title([num2str(num_clusters) ' clusters ' name_zone]);

            set(gca,'FontName','Times','Fontweight','Bold','FontSize',14)
            
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);

            %% Save

            print([fsave ht.String], '-dpng','-r0')

            save([fsave ht.String '.mat'],'time','KMeans_res','lonlat','waterlevel','-mat');
            
            

        end
    end
    
    close all
    
end


