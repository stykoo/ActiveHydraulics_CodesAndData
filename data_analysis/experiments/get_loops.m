% Obtain loops from experimental data
% Author: Camille Jorge <camille.jorge@ens-lyon.fr>

Initialize_20_01;

for thisExpe = 71%[1:87 96:102]

%% Charge la carte des flux

    main_path = ExperimentPath{1,thisExpe};
    Analysis_path = fullfile(main_path,'Analysis');
    Geometry_path = fullfile(Analysis_path,'1_Geometrie');
    
    filename_data = fullfile(Geometry_path,['hexagonal_map','.mat']);
    load(filename_data,'flux');


    X = flux(1,:);
    Y = flux(2,:);


    Y_tries = unique(sort(Y));
    

    all_loops = {'positions and flux vectors'};
    memory = [0;0];

    for i = 1 : length(X) %on explore le graphe colonne par colonne 
        j=i;
        if ~any(memory(1,:)==X(j) & memory(2,:)==Y(j)) && (flux(3,j)~=0 || flux(4,j)~=0) % choix d'un point jamais explor�


            memory(:,end+1) = [X(j) Y(j)];


            %% Construction de la boucle

            current_loop = [X(j); Y(j); flux(3,j); flux(4,j)];

            % D'abord on cherche les 4 (ou moins) plus proches voisins qui ne
            % sont pas sur le bord
            %neighbors = [X(sqrt((X-X(j)).^2+(Y-Y(j)).^2)<=0.87 & sqrt((X-X(j)).^2+(Y-Y(j)).^2)>=0.01 & ~any(memory(1,:)==X(j) & memory(2,:)==Y(j))); Y(sqrt((X-X(j)).^2+(Y-Y(j)).^2)<=0.87 & sqrt((X-X(j)).^2+(Y-Y(j)).^2)>=0.01 & ~any(memory(1,:)==X(j) & memory(2,:)==Y(j)))];
            neighbors = [X(sqrt((X-X(j)).^2+(Y-Y(j)).^2)<=1.2 & sqrt((X-X(j)).^2+(Y-Y(j)).^2)>=0.01); Y(sqrt((X-X(j)).^2+(Y-Y(j)).^2)<=1.2 & sqrt((X-X(j)).^2+(Y-Y(j)).^2)>=0.01)];
            admitted_neighbors =0;
            % Ensuite on regarde quels sont les points suivant dans la boucle

            for n = 1:size(neighbors,2)

                index = find(X == neighbors(1,n) & Y == neighbors(2,n));

                if ~any(memory(1,:)==X(index) & memory(2,:)==Y(index))
                    memory(:,end+1) = [X(index) Y(index)];%Si le point n'était pas en mémoire on l'ajoute dans la mémoire
                

                    if flux(3,index)~=0 || flux(4,index)~=0

                        current_loop(:,end+1) = [X(index) Y(index) flux(3,index) flux(4,index)];
                        admitted_neighbors= admitted_neighbors+1;

                    end


                end

            

            end


            % Dans la plupart des cas, on ajoute 2 points à la boucle. On
            % choisit le dernier pour explorer la boucle dans un sens

            
            j = find(X ==  current_loop(1,end) & Y == current_loop(2,end)); % On choisit comme nouveau point de départ le dernier point ajouté à la current_loop

            while admitted_neighbors > 0
                % Tant qu'on n'atteint pas un bord ou qu'on ne boucle pas on
                % continue l'exploration !

                
                neighbors = [X(sqrt((X-X(j)).^2+(Y-Y(j)).^2)<=1.2 & sqrt((X-X(j)).^2+(Y-Y(j)).^2)>=0.01); Y(sqrt((X-X(j)).^2+(Y-Y(j)).^2)<=1.2 & sqrt((X-X(j)).^2+(Y-Y(j)).^2)>=0.01)];
                admitted_neighbors =0;

                for n = 1:size(neighbors,2)

                    index = find(X == neighbors(1,n) & Y == neighbors(2,n));
                    
                    if ~any(memory(1,:)==X(index) & memory(2,:)==Y(index))

                        memory(:,end+1) = [X(index) Y(index)];%Si le point n'était pas en mémoire on l'ajoute dans la mémoire


                        if flux(3,index)~=0 || flux(4,index)~=0
                            % On compte le nombre de chemins attribués.
                            % si c'est 1 :ok
                            % si c'est 2 ou + : pb conservation de la masse
                            % si c'est 0 : stop le while
                            current_loop(:,end+1) = [X(index) Y(index) flux(3,index) flux(4,index)];
                            admitted_neighbors =admitted_neighbors+1;

                        end
                    end

                end

                if admitted_neighbors > 1
                    disp("ERROR MASS CONSERVATION OR WTF") %normalement on ne peut trouver qu'un seul ou zéro voisin qui convient
                end

                j = find(X ==  current_loop(1,end) & Y == current_loop(2,end)); %On passe au point suivant 


            end

            %% A ce moment, soit on a bouclé une boucle, soit on a parcouru une partie d'un polymère traversant
            % On va donc vérifier si (X(j),Y(j)) n'est pas plus proche voisin du premier point de current_loop

            %             if sqrt((X(j)-current_loop(1,1))^2 + (Y(j)-current_loop(2,1))^2) >= 0.87
            if size(current_loop,2)>1
                if sqrt((X(j)-current_loop(1,2))^2 + (Y(j)-current_loop(2,2))^2) >= 0.87
                    % si cette condition edt v�rifi�ee, les derniers et premiers
                    % points de current_loop ne sont pas plus  proches voisins, et
                    % nous avons affaire à un polymère traversant
                    
                    j = find(X ==  current_loop(1,2) & Y == current_loop(2,2));
                    % On choisit comme nouveau point de départ le deuxième point
                    % ajouté à current_loop, qui correspond à la direction pas
                    % encore explorée du polymère traversant.
                    % On recommence donc toutes les opérations effectuées pour la
                    % première partie du polymère traversant.
                    
                    neighbors = [X(sqrt((X-X(j)).^2+(Y-Y(j)).^2)<=1.2 & sqrt((X-X(j)).^2+(Y-Y(j)).^2)>=0.01 & ~any(memory(1,:)==X(j) & memory(2,:)==Y(j))); Y(~any(memory(1,:)==X(j) & memory(2,:)==Y(j)) & sqrt((X-X(j)).^2+(Y-Y(j)).^2)<=1.2& sqrt((X-X(j)).^2+(Y-Y(j)).^2)>=0.01)];
                    admitted_neighbors = 1;
                    
                    while admitted_neighbors > 0
                        % Tant qu'on n'atteint pas un bord ou qu'on ne boucle pas on
                        % continue l'exploration !
                        neighbors = [X(sqrt((X-X(j)).^2+(Y-Y(j)).^2)<=1.2 & sqrt((X-X(j)).^2+(Y-Y(j)).^2)>=0.01); Y(sqrt((X-X(j)).^2+(Y-Y(j)).^2)<=1.2 & sqrt((X-X(j)).^2+(Y-Y(j)).^2)>=0.01)];
                        admitted_neighbors =0;
                        
                        for n = 1:size(neighbors,2)
                            
                            index = find(X == neighbors(1,n) & Y == neighbors(2,n));
                            
                            if ~any(memory(1,:)==X(index) & memory(2,:)==Y(index))
                                memory(:,end+1) = [X(index) Y(index)];%Si le point n'était pas en mémoire on l'ajoute dans la mémoire
                                
                                if flux(3,index)~=0 || flux(4,index)~=0
                                    % On compte le nombre de chemins attribués.
                                    % si c'est 1 :ok
                                    % si c'est 2 ou + : pb conservation de la masse
                                    % si c'est 0 : stop le while
                                    current_loop(:,end+1) = [X(index) Y(index) flux(3,index) flux(4,index)];
                                    admitted_neighbors = admitted_neighbors+1;
                                    
                                end
                            end
                            
                            
                        end
                        
                        if admitted_neighbors > 1
                            disp("ERROR MASS CONSERVATION OR WTF")
                            % normalement on ne peut trouver qu'un seul ou zéro
                            % voisin qui convient
                        end
                        
                        
                        j = find(X ==  current_loop(1,end) & Y == current_loop(2,end));
                        
                        
                        
                    end
                end
            end

            % A priori à ce moment là current_loop contient la boucle/le
            % polymère traversant complet et on peut l'ajouter à all_loops.

            all_loops{end+1}= current_loop; %#ok<*SAGROW>


        elseif ~any(memory(1,:)==X(j) & memory(2,:)==Y(j)) && (flux(3,j)==0 && flux(4,j)==0)

            memory(:,end+1) = [X(j) Y(j)];

        end

    end


    filename_data = fullfile(Geometry_path,['loops','.mat']);
    save(filename_data,'all_loops');
    figure;
    title(num2str(thisExpe))
    for i = 2:length(all_loops)
        flux_test = all_loops{1,i};
        quiver(flux_test(1,:)-flux_test(3,:)./2,flux_test(2,:)-flux_test(4,:)./2,flux_test(3,:),flux_test(4,:),'linewidth',2,'AutoScale','off')
        hold on
        axis equal
    end
    
    hold off
end
