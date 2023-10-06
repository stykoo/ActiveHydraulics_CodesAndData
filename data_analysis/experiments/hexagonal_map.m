% Analyze experimental data on honeycomb lattice
% Author: Camille Jorge <camille.jorge@ens-lyon.fr>


for thisExpe = 7
    %% Rajoute un offset en x et y et rescale le r?seau pour que les centres des canaux des PIV concordent avec les centres des ar?tes
    
    Initialize_26_08;
    
    main_path = ExperimentPath{1,thisExpe};
    Images_path = fullfile(main_path,'Images');
    Analysis_path = fullfile(main_path,'Analysis');
    Geometry_path = fullfile(Analysis_path,'1_Geometrie');
    
    sizebox= 8;
    step = 4;
    dframe = 1;
    % stepTime = 10;
    % stepFrames = stepTime;
    prct= 0.5;
    
    filename_data = fullfile(Geometry_path,['CentroidsHexa_withborders','.mat']);
    load(filename_data,'centroids_tries');
    
    maxX = max(centroids_tries(:,1));
    maxY = max(centroids_tries(:,2));
    minX = min(centroids_tries(:,1));
    minY = min(centroids_tries(:,2));
    
    angle = 0*2*pi/360;% angle estimé du réseau filmé pas droit
    
    coin = centroids_tries(centroids_tries(:,2) == minY,:); %rep?re le point en bas ? gauche du r?seau
    y_tries = sort(centroids_tries(:,2));
    point = centroids_tries(centroids_tries(:,2) == y_tries(5),:);
    
    
    vect = point - coin;
    
    correction = 0*a/20;
    offset_hexa = [(minX)  (minY+a*sqrt(3)/2)];
    
    
    %% Cr?e un r?seau hexagonal de dimension "29*29" de la taille du mask
    
    center_x = [];
    center_y = [];
    
    
    x = 0;
    y=0;
    
    
%     index = 1:37; %r0.6-0.65 1:20
%     n_col = 32; %r0.6-0.65 1:20
    
    
%     index = 1:31; % r0.9 41:50
%     n_col = 30; % r0.9 41:50
%     
%     index = 1:31; %r0.8 21:30 r0.85 31:40
%     n_col = 29; %r0.8 21:30 r0.85 31:40
    
%     index = 1:29; %r0.95 51:60 r1.0 61:70
%     n_col = 28;  %r0.95 51:60 r1.0 61:70
    
    index = 1:27; %r1.05 71:80 r1.1 81:87
    n_col = 27;  %r1.05 71:80 r1.1 81:87
    
%     index = 1:33; % r0.75 96:102 r0.7 88:92
%     n_col = 31; % r0.75 96:102 r0.7 88: 92

%     index = 1:37; % r0.55 104:113
%     n_col = 30; % r0.55 104:113

%     index = 1:35; % r0.75 104:113
%     n_col = 31; % r0.75 104:113
     
    for colonne = index
        if mod(colonne,2)==0
            nbr = n_col-1;
            y = 0;
        else
            nbr = n_col;
            y = -sqrt(3)*a/2;
        end
        
        for cnt = [0 : a :  nbr*a ]
            center_x(end+1) = x;
            center_y(end+1) = y;
            
            y = y +sqrt(3)*a;
        end
        
        x = x +a*1.5;
        
    end
    
    %% Modifie un peu les bords du réseau pour s'affranchir des aberrations
    
    center_x_tries = sort(unique(center_x));
    center_y_tries = sort(unique(center_y));
    
    % correction sur X
    
%     correction = 0;
    
    %minX
    
    % plus petits y
    center_x(center_x == center_x_tries(1) & center_y < center_y_tries(10)) = center_x(center_x == center_x_tries(1) & center_y < center_y_tries(10)) + correction;
    center_y(center_x == center_x_tries(1) & center_y < center_y_tries(10)) = center_y(center_x == center_x_tries(1) & center_y < center_y_tries(10)) + correction;
    
    
    % plus grands y
    center_x(center_x == center_x_tries(1) & center_y > center_y_tries(end-10)) = center_x(center_x == center_x_tries(1) & center_y > center_y_tries(end-10)) + correction;
    center_y(center_x == center_x_tries(1) & center_y > center_y_tries(end-10)) = center_y(center_x == center_x_tries(1) & center_y > center_y_tries(end-10)) - correction;
    
    
    
    % min2X
    
    % plus petits y
    center_x(center_x == center_x_tries(2) & center_y < center_y_tries(10)) = center_x(center_x == center_x_tries(2) & center_y < center_y_tries(10)) + correction;
    center_y(center_x == center_x_tries(2) & center_y < center_y_tries(10)) = center_y(center_x == center_x_tries(2) & center_y < center_y_tries(10)) + correction;
    
    
    
    % plus grands y
    center_x(center_x == center_x_tries(2) & center_y > center_y_tries(end-10)) = center_x(center_x == center_x_tries(2) & center_y > center_y_tries(end-10)) + correction;
    center_y(center_x == center_x_tries(2) & center_y > center_y_tries(end-10)) = center_y(center_x == center_x_tries(2) & center_y > center_y_tries(end-10)) - correction;
    

    %maxX
    
    % plus petits y
    center_x(center_x == center_x_tries(end) & center_y < center_y_tries(10)) = center_x(center_x == center_x_tries(end) & center_y < center_y_tries(10)) - correction;
    center_y(center_x == center_x_tries(end) & center_y < center_y_tries(10)) = center_y(center_x == center_x_tries(end) & center_y < center_y_tries(10)) + correction;
    
    % plus grands y
    center_x(center_x == center_x_tries(end) & center_y > center_y_tries(end-10)) = center_x(center_x == center_x_tries(end) & center_y > center_y_tries(end-10)) - correction;
    center_y(center_x == center_x_tries(end) & center_y > center_y_tries(end-10)) = center_y(center_x == center_x_tries(end) & center_y > center_y_tries(end-10)) - correction;
    
    
    %max2X
    
    % plus petits y
    
    center_x(center_x == center_x_tries(end-1) & center_y < center_y_tries(10)) = center_x(center_x == center_x_tries(end-1) & center_y < center_y_tries(10)) - correction;
    center_y(center_x == center_x_tries(end-1) & center_y < center_y_tries(10)) = center_y(center_x == center_x_tries(end-1) & center_y < center_y_tries(10)) + correction;
    
    % plus grands y
    
    center_x(center_x == center_x_tries(end-1) & center_y> center_y_tries(end-10)) = center_x(center_x == center_x_tries(end-1) & center_y > center_y_tries(end-10)) - correction;
    center_y(center_x == center_x_tries(end-1) & center_y > center_y_tries(end-10)) = center_y(center_x == center_x_tries(end-1) & center_y > center_y_tries(end-10)) - correction;
    
    
    % correction sur Y
%         correction = 2*a/20;
    % minY
    
    % plus petits x
    
    center_x((center_y == center_y_tries(1) |center_y == center_y_tries(1)+correction) & center_x < center_x_tries(8)) = center_x((center_y == center_y_tries(1) |center_y == center_y_tries(1)+correction) & center_x < center_x_tries(8)) + correction;
    center_y((center_y == center_y_tries(1) |center_y == center_y_tries(1)+correction) & center_x < center_x_tries(8)) = center_y((center_y == center_y_tries(1) |center_y == center_y_tries(1)+correction) & center_x < center_x_tries(8)) + correction;
    
    
    % plus grands x
    
    center_x((center_y == center_y_tries(1) |center_y == center_y_tries(1)+correction) & center_x > center_x_tries(end-8)) = center_x((center_y == center_y_tries(1) |center_y == center_y_tries(1)+correction) & center_x > center_x_tries(end-8)) - correction;
    center_y((center_y == center_y_tries(1) |center_y == center_y_tries(1)+correction) & center_x > center_x_tries(end-8)) = center_y((center_y == center_y_tries(1) |center_y == center_y_tries(1)+correction) & center_x > center_x_tries(end-8)) + correction;
    
    %min2y
    
    % plus petits x
    
    center_x((center_y == center_y_tries(2) |center_y == center_y_tries(2)+correction) & center_x < center_x_tries(8)) = center_x((center_y == center_y_tries(2) |center_y == center_y_tries(2)+correction) & center_x < center_x_tries(8)) + correction;
    center_y((center_y == center_y_tries(2) |center_y == center_y_tries(2)+correction) & center_x < center_x_tries(8)) = center_y((center_y == center_y_tries(2) |center_y == center_y_tries(2)+correction) & center_x < center_x_tries(8)) + correction;
    
    
    % plus grands x
    
    center_x((center_y == center_y_tries(2) |center_y == center_y_tries(2)+correction) & center_x > center_x_tries(end-7)) = center_x((center_y == center_y_tries(2) |center_y == center_y_tries(2)+correction) & center_x > center_x_tries(end-7)) - correction;
    center_y((center_y == center_y_tries(2) |center_y == center_y_tries(2)+correction) & center_x > center_x_tries(end-7)) = center_y((center_y == center_y_tries(2) |center_y == center_y_tries(2)+correction) & center_x > center_x_tries(end-7)) + correction;
    
    correction = correction/2; %r1.0
    % maxY
    
    % plus petits x
    
    center_x((center_y == center_y_tries(end) |center_y == center_y_tries(end)-correction) & center_x < center_x_tries(8)) = center_x((center_y == center_y_tries(end) |center_y == center_y_tries(end)-correction) & center_x < center_x_tries(8)) + correction;
    center_y((center_y == center_y_tries(end) |center_y == center_y_tries(end)-correction) & center_x < center_x_tries(8)) = center_y((center_y == center_y_tries(end) |center_y == center_y_tries(end)-correction) & center_x < center_x_tries(8)) - correction;
    
    
    % plus grands x
    
    center_x((center_y == center_y_tries(end) |center_y == center_y_tries(end)+correction) & center_x > center_x_tries(end-7)) = center_x((center_y == center_y_tries(end) |center_y == center_y_tries(end)+correction) & center_x > center_x_tries(end-7)) - correction;
    center_y((center_y == center_y_tries(end) |center_y == center_y_tries(end)+correction) & center_x > center_x_tries(end-7)) = center_y((center_y == center_y_tries(end) |center_y == center_y_tries(end)+correction) & center_x > center_x_tries(end-7)) - correction;
    
    %max2Y
    
    % plus petits x
    
    center_x((center_y == center_y_tries(end-1) |center_y == center_y_tries(end-1)-correction) & center_x < center_x_tries(8)) = center_x((center_y == center_y_tries(end-1) |center_y == center_y_tries(end-1)-correction) & center_x < center_x_tries(8)) + correction;
    center_y((center_y == center_y_tries(end-1) |center_y == center_y_tries(end-1)-correction) & center_x < center_x_tries(8)) = center_y((center_y == center_y_tries(end-1) |center_y == center_y_tries(end-1)-correction) & center_x < center_x_tries(8)) - correction;
    
    
    % plus grands x
    
    center_x((center_y == center_y_tries(end-1) |center_y == center_y_tries(end-1)+correction) & center_x > center_x_tries(end-7)) = center_x((center_y == center_y_tries(end-1) |center_y == center_y_tries(end-1)+correction) & center_x > center_x_tries(end-7)) - correction;
    center_y((center_y == center_y_tries(end-1) |center_y == center_y_tries(end-1)+correction) & center_x > center_x_tries(end-7)) = center_y((center_y == center_y_tries(end-1) |center_y == center_y_tries(end-1)+correction) & center_x > center_x_tries(end-7)) - correction;
    
    
    %% Donne les centres des arr?tes des hexagones
    
    centers = [center_x+offset_hexa(1) ;center_y+offset_hexa(2)];
    % centers = [center_x ;center_y];
    edges_centers = [;];
    
    bord_gauche = centers(1,:)==min(centers(1,:));
    bord_bas = centers(2,:)==min(centers(2,:));
    bord_droit = centers(1,:)==max(centers(1,:));
    bord_haut = centers(2,:)==max(centers(2,:));
    
    
    % edge_center_1 = a*[0; sqrt(3)/2];
    % edge_center_2 = a*[-3/4; sqrt(3)/4];
    % edge_center_3 = a*[-3/4; -sqrt(3)/4];
    % edge_center_4 = a*[0; -sqrt(3)/2];
    % edge_center_5 = a*[3/4; -sqrt(3)/4];
    % edge_center_6 = a*[3/4; sqrt(3)/4];
    
    edge_center_1 = a*[0; sqrt(3)/2];
    edge_center_2 = a*[-3/4; sqrt(3)/4];
    edge_center_3 = a*[-3/4; -sqrt(3)/4];
    edge_center_4 = a*[0; -sqrt(3)/2];
    edge_center_5 = a*[3/4; -sqrt(3)/4];
    edge_center_6 = a*[3/4; sqrt(3)/4];
    
    for hexa = 1:length(centers)
        
        %     if bord_gauche(hexa) == 1
        %
        %         edge_center_1 = a*[0; sqrt(3)/2];
        %         edge_center_2 = 0.8*a*[-3/4; sqrt(3)/4];
        %         edge_center_3 = 0.8*a*[-3/4; -sqrt(3)/4];
        %         edge_center_4 = a*[0; -sqrt(3)/2];
        %         edge_center_5 = a*[3/4; -sqrt(3)/4];
        %         edge_center_6 = a*[3/4; sqrt(3)/4];
        %
        %     elseif bord_droit(hexa) == 1
        %
        %         edge_center_1 = a*[0; sqrt(3)/2];
        %         edge_center_2 = a*[-3/4; sqrt(3)/4];
        %         edge_center_3 = a*[-3/4; -sqrt(3)/4];
        %         edge_center_4 = a*[0; -sqrt(3)/2];
        %         edge_center_5 = 0.8*a*[3/4; -sqrt(3)/4];
        %         edge_center_6 = 0.8*a*[3/4; sqrt(3)/4];
        %
        %     elseif bord_bas(hexa) == 1
        %
        %         edge_center_1 = a*[0; sqrt(3)/2];
        %         edge_center_2 = a*[-3/4; sqrt(3)/4];
        %         edge_center_3 = 0.8*a*[-3/4; -sqrt(3)/4];
        %         edge_center_4 = 0.8*a*[0; -sqrt(3)/2];
        %         edge_center_5 = 0.8*a*[3/4; -sqrt(3)/4];
        %         edge_center_6 = a*[3/4; sqrt(3)/4];
        %
        %     elseif bord_haut(hexa) == 1
        %
        %         edge_center_1 = 0.8*a*[0; sqrt(3)/2];
        %         edge_center_2 = 0.8*a*[-3/4; sqrt(3)/4];
        %         edge_center_3 = a*[-3/4; -sqrt(3)/4];
        %         edge_center_4 = a*[0; -sqrt(3)/2];
        %         edge_center_5 = a*[3/4; -sqrt(3)/4];
        %         edge_center_6 = 0.8*a*[3/4; sqrt(3)/4];
        %
        %     else
        %
        %         edge_center_1 = a*[0; sqrt(3)/2];
        %         edge_center_2 = a*[-3/4; sqrt(3)/4];
        %         edge_center_3 = a*[-3/4; -sqrt(3)/4];
        %         edge_center_4 = a*[0; -sqrt(3)/2];
        %         edge_center_5 = a*[3/4; -sqrt(3)/4];
        %         edge_center_6 = a*[3/4; sqrt(3)/4];
        %
        %     end
        
        
        if hexa ==1
            
            edges_centers(:,end+1) = centers(:,hexa) + edge_center_1;
            edges_centers(:,end+1) = centers(:,hexa) + edge_center_2;
            edges_centers(:,end+1) = centers(:,hexa) + edge_center_3;
            edges_centers(:,end+1) = centers(:,hexa) + edge_center_4;
            edges_centers(:,end+1) = centers(:,hexa) + edge_center_5;
            edges_centers(:,end+1) = centers(:,hexa) + edge_center_6;
            
        else
            
            d1 = edges_centers-(centers(:,hexa) + edge_center_1);
            d2 = edges_centers-(centers(:,hexa) + edge_center_2);
            d3 = edges_centers-(centers(:,hexa) + edge_center_3);
            d4 = edges_centers-(centers(:,hexa) + edge_center_4);
            d5 = edges_centers-(centers(:,hexa) + edge_center_5);
            d6 = edges_centers-(centers(:,hexa) + edge_center_6);
            
            if sum(sqrt(d1(1,:).^2+d1(2,:).^2)<30)==0
                edges_centers(:,end+1) = centers(:,hexa) + edge_center_1;
            end
            if sum(sqrt(d2(1,:).^2+d2(2,:).^2)<30)==0
                edges_centers(:,end+1) = centers(:,hexa) + edge_center_2;
            end
            if sum(sqrt(d3(1,:).^2+d3(2,:).^2)<30)==0
                edges_centers(:,end+1) = centers(:,hexa) + edge_center_3;
            end
            if sum(sqrt(d4(1,:).^2+d4(2,:).^2)<30)==0
                edges_centers(:,end+1) = centers(:,hexa) + edge_center_4;
            end
            if sum(sqrt(d5(1,:).^2+d5(2,:).^2)<30)==0
                edges_centers(:,end+1) = centers(:,hexa) + edge_center_5;
            end
            if sum(sqrt(d6(1,:).^2+d6(2,:).^2)<30)==0
                edges_centers(:,end+1) = centers(:,hexa) + edge_center_6;
            end            
            
        end
        
        
    end

edges_centers_unique = transpose(unique(transpose(round(edges_centers,8)),'Rows'));
    

figure;
scatter(centers(1,:),centers(2,:))
axis equal
hold on 
scatter(edges_centers_unique(1,:),edges_centers_unique(2,:),'*')


%% Tourne le r?seau de l'angle du r?seau pour coller ? la PIV

% angle = acos(dot(vect,[1 0])/norm(vect)); 
% angle = -0.01; 
% %%%% !!!!! ATTENTION !!!!!
% %%%% Il faut tourner le réseau par rapport à son point en bas à gauche et
% %%%% pas par rapport à l'origine de la base !
% % Vire l'offset et place le centre en 0 
% centers_angle(1,:) = centers(1,:) - offset_hexa(1) - mean(centers(1,:) - offset_hexa(1));
% centers_angle(2,:) = centers(2,:) - offset_hexa(2) - mean(centers(2,:) - offset_hexa(2));
% edges_centers_unique_angle(1,:) = edges_centers_unique(1,:) - offset_hexa(1) - mean(centers(1,:) - offset_hexa(1));
% edges_centers_unique_angle(2,:) = edges_centers_unique(2,:) - offset_hexa(2) - mean(centers(2,:) - offset_hexa(2));
% % % Tourne 
% centers_angle = centers_angle.*[cos(angle) ;-sin(angle)] + centers_angle.* [sin(angle) ;cos(angle)];
% edges_centers_unique_angle = edges_centers_unique_angle .* [cos(angle) ;-sin(angle)] +   edges_centers_unique_angle .* [sin(angle) ;cos(angle)];
% % % remet l'offset
% centers_angle(1,:) = centers_angle(1,:) + offset_hexa(1) + mean(centers(1,:) - offset_hexa(1));
% centers_angle(2,:) = centers_angle(2,:) + offset_hexa(2) + mean(centers(2,:) - offset_hexa(2));
% edges_centers_unique_angle(1,:) = edges_centers_unique_angle(1,:) + offset_hexa(1) + mean(centers(1,:) - offset_hexa(1));
% edges_centers_unique_angle(2,:) = edges_centers_unique_angle(2,:) + offset_hexa(2) + mean(centers(2,:) - offset_hexa(2));

% Quand on veut pas d'angle
centers_angle = centers;
edges_centers_unique_angle = edges_centers_unique;

Mask_path = fullfile(Images_path,'MaskWhite.tif');
Mask = imread(Mask_path);

figure;
imagesc(Mask)
set(gca, 'YDir', 'normal')
hold on 
scatter(centers_angle(1,:),centers_angle(2,:))
axis equal
hold on 
scatter(edges_centers_unique_angle(1,:),edges_centers_unique_angle(2,:),'*')


%% Cr?e une image de la taille de la PIV

PIVLab_data_path=fullfile(Analysis_path,'PIVLab',['PIVLab_postprocess_',num2str(sizebox),'_dframe_',num2str(1),'_step_',num2str(step),'.mat']);
load(PIVLab_data_path,'U_filt','V_filt')
filename_data4 = fullfile(Geometry_path,['EcoulementBinaireCirculation_double_MeanTimeAll_VecteursCanaux_prct_',num2str(prct),'_box_',num2str(sizebox),'_step_',num2str(step),'_withborders','.mat']);
load(filename_data4,'MaskFinal')


facteur1 = length(V_filt(:,1))/size(Mask,1);
facteur2 = size(V_filt(1,:,1),2)/size(Mask,2);

centers_angle_PIV = centers_angle .* [facteur1; facteur2];
edges_centers_unique_angle_PIV = edges_centers_unique_angle .* [facteur1; facteur2];

U = nanmean(U_filt,3);
V = nanmean(V_filt,3);



mean_x = mean(centers_angle_PIV(1,:));
mean_y = mean(centers_angle_PIV(2,:));
centers_angle_PIV(1,:) = centers_angle_PIV(1,:) -mean_x;
centers_angle_PIV(2,:) = centers_angle_PIV(2,:) -mean_y;
edges_centers_unique_angle_PIV(1,:) = edges_centers_unique_angle_PIV(1,:) -mean_x;
edges_centers_unique_angle_PIV(2,:) = edges_centers_unique_angle_PIV(2,:) -mean_y;

figure;
scatter(centers_angle_PIV(1,:)*cos(angle)+centers_angle_PIV(2,:)*sin(angle)+mean_x,centers_angle_PIV(2,:)*cos(angle)-centers_angle_PIV(1,:)*sin(angle)+mean_y)
axis equal
hold on 
scatter(edges_centers_unique_angle_PIV(1,:)*cos(angle)+edges_centers_unique_angle_PIV(2,:)*sin(angle)+mean_x,edges_centers_unique_angle_PIV(2,:)*cos(angle)-edges_centers_unique_angle_PIV(1,:)*sin(angle)+mean_y,'*')
hold on 
quiver(U.*MaskFinal,V.*MaskFinal,'k')
% quiver(U,V,'k')


%% Calcule la vitesse moyenne autour de chaque centre de canal 

flux_X = [];
flux_Y = [];
U_flow = U.*MaskFinal;
V_flow =V.*MaskFinal;
U_flow(isnan(U_flow)==1)=0;
V_flow(isnan(V_flow)==1)=0;


for channel = 1:length(edges_centers_unique_angle_PIV)
%     y = round(edges_centers_unique_angle_PIV(1,channel)); % wtf les x et les y sont invers?s ? je pige r
%     x = round(edges_centers_unique_angle_PIV(2,channel));

    y = round(edges_centers_unique_angle_PIV(1,channel)*cos(angle)-edges_centers_unique_angle_PIV(2,channel)*sin(angle)+mean_x); % wtf les x et les y sont invers?s ? je pige r
    x = round(edges_centers_unique_angle_PIV(2,channel)*cos(angle)+edges_centers_unique_angle_PIV(1,channel)*sin(angle)+mean_y);
    
    flux_X(end+1) = (U_flow(x,y) + U_flow(x+1,y) + U_flow(x+1,y+1) + U_flow(x,y+1) + U_flow(x-1,y+1) + U_flow(x-1,y) + U_flow(x-1,y-1) + U_flow(x,y-1) + U_flow(x+1,y-1))/9;% + U_flow(x+2,y)+ U_flow(x+2,y+1)+ U_flow(x+2,y+2)+ U_flow(x+1,y+2)+ U_flow(x,y+2)+ U_flow(x-1,y+2)+ U_flow(x-2,y+2)+ U_flow(x-2,y+1)+ U_flow(x-2,y)+ U_flow(x-2,y-1)+ U_flow(x-2,y-2)+ U_flow(x-1,y-2)+ U_flow(x,y-2)+ U_flow(x+1,y-2)+ U_flow(x+2,y-2)+ U_flow(x+2,y-1))/25;
    flux_Y(end+1) = (V_flow(x,y) + V_flow(x+1,y) + V_flow(x+1,y+1) + V_flow(x,y+1) + V_flow(x-1,y+1) + V_flow(x-1,y) + V_flow(x-1,y-1) + V_flow(x,y-1) + V_flow(x+1,y-1))/9;%  + V_flow(x+2,y)+ V_flow(x+2,y+1)+ V_flow(x+2,y+2)+ V_flow(x+1,y+2)+ V_flow(x,y+2)+ V_flow(x-1,y+2)+ V_flow(x-2,y+2)+ V_flow(x-2,y+1)+ V_flow(x-2,y)+ V_flow(x-2,y-1)+ V_flow(x-2,y-2)+ V_flow(x-1,y-2)+ V_flow(x,y-2)+ V_flow(x+1,y-2)+ V_flow(x+2,y-2)+ V_flow(x+2,y-1))/25;
    
%     flux_X(end+1) = (U_flow(x,y) + U_flow(x+1,y) + U_flow(x+1,y+1) + U_flow(x,y+1) + U_flow(x-1,y+1) + U_flow(x-1,y) + U_flow(x-1,y-1) + U_flow(x,y-1) + U_flow(x+1,y-1)+ U_flow(x+2,y)+ U_flow(x+2,y+1)+ U_flow(x+2,y+2)+ U_flow(x+1,y+2)+ U_flow(x,y+2)+ U_flow(x-1,y+2)+ U_flow(x-2,y+2)+ U_flow(x-2,y+1)+ U_flow(x-2,y)+ U_flow(x-2,y-1)+ U_flow(x-2,y-2)+ U_flow(x-1,y-2)+ U_flow(x,y-2)+ U_flow(x+1,y-2)+ U_flow(x+2,y-2)+ U_flow(x+2,y-1))/25;
%     flux_Y(end+1) = (V_flow(x,y) + V_flow(x+1,y) + V_flow(x+1,y+1) + V_flow(x,y+1) + V_flow(x-1,y+1) + V_flow(x-1,y) + V_flow(x-1,y-1) + V_flow(x,y-1) + V_flow(x+1,y-1)+ V_flow(x+2,y)+ V_flow(x+2,y+1)+ V_flow(x+2,y+2)+ V_flow(x+1,y+2)+ V_flow(x,y+2)+ V_flow(x-1,y+2)+ V_flow(x-2,y+2)+ V_flow(x-2,y+1)+ V_flow(x-2,y)+ V_flow(x-2,y-1)+ V_flow(x-2,y-2)+ V_flow(x-1,y-2)+ V_flow(x,y-2)+ V_flow(x+1,y-2)+ V_flow(x+2,y-2)+ V_flow(x+2,y-1))/25;
    
%     flux_X(end+1) = (U_flow(x,y) + U_flow(x+1,y) + U_flow(x+1,y+1) + U_flow(x,y+1) + U_flow(x-1,y+1) + U_flow(x-1,y) + U_flow(x-1,y-1) + U_flow(x,y-1) + U_flow(x+1,y-1)+ U_flow(x+2,y)+ U_flow(x+2,y+1)+ U_flow(x+1,y+2)+ U_flow(x,y+2)+ U_flow(x-1,y+2)+ U_flow(x-2,y+1)+ U_flow(x-2,y)+ U_flow(x-2,y-1)+ U_flow(x-1,y-2)+ U_flow(x,y-2)+ U_flow(x+1,y-2)+ U_flow(x+2,y-1))/13;
%     flux_Y(end+1) = (V_flow(x,y) + V_flow(x+1,y) + V_flow(x+1,y+1) + V_flow(x,y+1) + V_flow(x-1,y+1) + V_flow(x-1,y) + V_flow(x-1,y-1) + V_flow(x,y-1) + V_flow(x+1,y-1)+ V_flow(x+2,y)+ V_flow(x+2,y+1)+ V_flow(x+1,y+2)+ V_flow(x,y+2)+ V_flow(x-1,y+2)+ V_flow(x-2,y+1)+ V_flow(x-2,y)+ V_flow(x-2,y-1)+ V_flow(x-1,y-2)+ V_flow(x,y-2)+ V_flow(x+1,y-2)+ V_flow(x+2,y-1))/13;

end


%% Attribue ? chaque vecteur flux un des 6 vecteurs possibles autour de l'hexagone (de norme 1)

flux_X_clean = [];
flux_Y_clean = [];
for channel = [1:length(flux_X)]
    
    if flux_X(channel)~=0 || flux_Y(channel)~=0
        liste_dot_product = [(flux_X(channel)) (flux_X(channel)/2 +flux_Y(channel)*sqrt(3)/2) (-flux_X(channel)/2 + flux_Y(channel)*sqrt(3)/2)];
        direction = find(abs(liste_dot_product)==max(abs(liste_dot_product)));
        sens = sign(liste_dot_product(direction));
    else
        direction=0;
        sens=0;
    end
    
    if direction==1
        flux_X_clean(end+1) = 1 * sens;
        flux_Y_clean(end+1) = 0;
    elseif direction == 2
        flux_X_clean(end+1) = 1/2 * sens;
        flux_Y_clean(end+1) = sqrt(3)/2 * sens;
    elseif direction == 3
        flux_X_clean(end+1) = -1/2 * sens;
        flux_Y_clean(end+1) = sqrt(3)/2 * sens;
    elseif direction == 0
        flux_X_clean(end+1) = 0;
        flux_Y_clean(end+1) = 0;
    end
end

%% Nouveau r?seau de pas 1 pour que toutes les figures soient de la m?me taille

edges_centers_unique = edges_centers_unique./a;

X = edges_centers_unique(1,:);
Y = edges_centers_unique(2,:);

figure;
quiver(X-flux_X_clean./2,Y-flux_Y_clean./2,flux_X_clean,flux_Y_clean,'linewidth',2,'color','r')
axis equal
hold on 
scatter(centers(1,:)./a,centers(2,:)./a,'*','k')

flux = [X; Y; flux_X_clean; flux_Y_clean];

filename_data = fullfile(Geometry_path,['hexagonal_map','.mat']);
save(filename_data,'flux'); 

end

% close all
%clear all
