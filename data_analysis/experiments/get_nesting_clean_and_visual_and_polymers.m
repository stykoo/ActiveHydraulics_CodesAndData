% Obtain nesting level from experimental data
% Author: Camille Jorge <camille.jorge@ens-lyon.fr>

Initialize_20_01;

r_liste = [];
replica_liste = [];
var_liste = [];
mean_liste =[];
max_liste = [];
min_liste = [];
for z = 4%1:13%[1:3 5:11]
    if z == 1
        liste_index = 2;%3:10;%1:10;
        r = 0.6;
        coeff = 33*37/(33*37);
        l=37;
        b=33;
        tit = 'r0.6';
    elseif z == 2
        liste_index = 11:20;
        r = 0.65;
        coeff = 33*37/(33*37);
        l=37;
        b=33;
        tit = 'r0.65';
    elseif z == 3
        liste_index = 88:90;
        r = 0.7;
        coeff = (33*37)/(33*32);
        l=33;
        b=32;
        tit = 'r0.7';
    elseif z ==4
        liste_index = 96:102;
        tit = 'r0.75';
        r = 0.75;
        coeff = (33*37)/(33*32);
        l=33;
        b=32;
    elseif z ==5
        liste_index = 21:30;
        r = 0.8;
        coeff = (33*37)/(30*31);
        l=31;
        b=30;
        tit = 'r0.8';
    elseif z ==6
        liste_index = 31:40;
        r = 0.85;
        coeff = (33*37)/(30*31);
        l=31;
        b=30;
        tit = 'r0.85';
    elseif z ==7
        liste_index = 41:50;
        r = 0.9;
        l=31;
        b=31;
        coeff = (33*37)/(31*31);
        tit = 'r0.9';
    elseif z ==8
        liste_index = 51:60;
        r = 0.95;
        coeff = (33*37)/(29*29);
        l=29;
        b=29;
        tit = 'r0.95';
    elseif z ==9
        liste_index = 61:70;
        r = 1.0;
        coeff = (33*37)/(29*29);
        l=29;
        b=29;
        tit = 'r1.0';
    elseif z ==10
        liste_index = 71:80;
        r = 1.05;
        coeff = (33*37)/(28*27);
        l=27;
        b=28;
        tit = 'r1.05';
    elseif z ==11
        liste_index = 81:87;
        r = 1.1;
        coeff = (33*37)/(28*27);
        l=27;
        b=28;
        tit = 'r1.1';
    elseif z ==12
        tit = 'r0.55';
        liste_index = 104:113;
        coeff = (33*37)/(37*30);
        l=47;
        b=30;
        tit = 'r1.1';
        r=0.55;
    elseif z ==13
        tit = 'r0.625';
        liste_index = 114:123;
        coeff = (33*37)/(35*31);
        l=45;
        b=31;
        tit = 'r1.1';
        r=0.625;
    end
    
    
    
    for thisExpe = liste_index%[21:30]%[11:20 66:75 84:93 94:103 104:113 114:121]%[1:10 41:50 76:83 ]
        
        replica_liste(end+1) = thisExpe-liste_index(1)+1;
        
        main_path = ExperimentPath{1,thisExpe};
        Analysis_path = fullfile(main_path,'Analysis');
        Geometry_path = fullfile(Analysis_path,'1_Geometrie');
        
        filename_data = fullfile(Geometry_path,['hexagonal_map','.mat']);
        load(filename_data,'flux');
        
        
        %% Crée la matrice C correspondant au nesting
        
        C=zeros(l,b);
        
        i = 0;
        
        % figure;
        
        liste_x = unique(flux(1,:));
        liste_y = unique(flux(2,:));
        
        disp(length(2:2:length(unique(flux(1,:)))))
        disp(length(unique(flux(1,:))))
        
        for x = 2:2:length(unique(flux(1,:)))
            
            
            i = i+1;
            j = 0;
            
            for y = 1:2:length(unique(flux(2,:)))-1
%             for y = 1:length(unique(flux(2,:)))-1
                
                
                
                if ~isempty(find(flux(1,:)==liste_x(x) & flux(2,:)==liste_y(y), 1))
                    
                    index = find(flux(1,:)==liste_x(x) & flux(2,:)==liste_y(y));
                    
                    j = j+1;
                    
                    if flux(3,index) > 0
                        
                        C(i,1:b >=j) = C(i,1:b >=j) +1 ;
                        
                        
                    elseif flux(3,index) < 0
                        
                        C(i,1:b >=j) = C(i,1:b >=j) -1 ;
                        
                    end
                    
                end
                
                
            end
            
        end
        
        
        %% Crée un patchwork hexagonal de la taille du réseau
        
        % La couleur du patchwork est donnée par une matrice C
        
        xhex=[1 0.5 -0.5 -1 -0.5 0.5]*1.2; % x-coordinates of the vertices 1.35
        yhex=[0 sqrt(3)/2 sqrt(3)/2 0 -sqrt(3)/2 -sqrt(3)/2]*1.2; % y-coordinates of the vertices 1.17
        
        
        
        figure;
        for k=1:l
            m=k-1;
            for i=mod(m,2)+1:b
                j=i-1;
                
                                if mod(k,2)==1
                    patch(xhex+sqrt(3)*m,(yhex+mod(k,2))+2*j,C(k,i)*coeff,'LineStyle','none') % test pour ré-aligner
                else
                    patch(xhex+sqrt(3)*m,(yhex+mod(k,2))+2*j,C(k,i-1)*coeff,'LineStyle','none')
                end
%                     patch(xhex+sqrt(3)*m,(yhex+mod(k,2))+2*j,C(k,i-(1-mod(k,2)))*coeff,'LineStyle','none')
                hold on
            end
        end
        axis equal
        colorbar()
        caxis([-12 12])
        colormap(jet)
        
%         hold off 

        flux = flux *2/sqrt(3);
        
        quiver(flux(1,:)-flux(3,:)/2,flux(2,:)-flux(4,:)/2 +2,flux(3,:),flux(4,:),'k','linewidth',2,'ShowArrowHead','off')
        set_plot()
        
%         pathname = 'A:\Camille\Résultats boucles 01_03\nesting\cartes';
%         filename = ['cartes niveau + boucles ',num2str(tit),' ',num2str(thisExpe - liste_index(1)+1),'.png'];
%         saveas(gcf,fullfile(pathname,filename));
%         filename = ['cartes niveau + boucles ',num2str(tit),' ',num2str(thisExpe - liste_index(1)+1),'.fig'];
%         saveas(gcf,fullfile(pathname,filename));

%         pathname = 'A:\Camille\films_plots_propre\contour_lines';
%         filename = ['contour_lines ',num2str(tit),' ',num2str(thisExpe-min(liste_index)+1),'.pdf'];
%         saveas(gcf,fullfile(pathname,filename));

        
%         var_liste(end+1) = var(C,1,"all");
%         mean_liste(end+1)= mean((C*coeff).^2,'all');
        r_liste(end+1) = r;
        max_liste(end+1) = max(max(C*coeff));
        min_liste(end+1) = min(min(C*coeff));
        
%         figure;
%         histogram(C,'Normalization','Probability')
        
        
    end
    
    
    
    
end

r_liste_mean = [0.6 0.65 0.7 0.8 0.85 0.9 0.95 1.0 1.05 1.1 0.55 0.625];
max_liste_mean = [];
diff_liste_mean = [];
var_liste_mean = [];
mean_liste_mean = [];
% 
for i = r_liste_mean
%     var_liste_mean(end+1) = mean(var_liste(r_liste==i));
%     mean_liste_mean(end+1) =  mean(mean_liste(r_liste==i));
%     max_liste_mean(end+1) =  mean(max_liste(r_liste==i));
    diff_liste_mean(end+1) =  mean(max_liste(r_liste==i)-min_liste(r_liste==i));
%     
end

r = [0.66 0.7 0.74 0.77 0.84 0.99 1.07 1.14 1.22 1.3 1.44];
figure;
plot(r,diff_liste_mean([11 1 12 2:8 10]),'o--')
xlabel('$\epsilon$', 'interpreter', 'latex')
ylabel('$\rm Degrees ~of ~nesting$', 'interpreter', 'latex')
% legend('repliques','moyenne')
ylim([1 12])
xlim([0.6 1.5])
set_plot()

% % 
% figure;
% plot(r_liste,var_liste,'o',r_liste_mean,var_liste_mean,'o--')
% xlabel('$\rm aspect ~ratio$', 'interpreter', 'latex')
% ylabel('$<\delta h^2>$', 'interpreter', 'latex')
% legend('repliques','moyenne')
% set_plot()

% figure;
% plot(r_liste,mean_liste,'o',r_liste_mean,mean_liste_mean,'o--')
% xlabel('$\rm aspect ~ratio$', 'interpreter', 'latex')
% ylabel('$<h^2>$', 'interpreter', 'latex')
% legend('repliques','moyenne')
% ylim([0 10])
% set_plot()

% figure;
% for index = 1:10
%     
%     plot(r_liste(replica_liste == index),max_liste(replica_liste == index),'o','Color',[index/10, index/10, index/10, 0.2])
%     hold on
%     
%     
% end
% 
% xlabel('$\rm aspect ~ratio$', 'interpreter', 'latex')
% ylabel('\rm max(h)', 'interpreter', 'latex')
% 
% set_plot()

% figure;
% plot(r_liste,max_liste,'o',r_liste_mean,max_liste_mean,'o--')
% xlabel('$\rm aspect ~ratio$', 'interpreter', 'latex')
% ylabel('\rm max(h)', 'interpreter', 'latex')
% legend('repliques','moyenne')
% ylim([0 10])
% set_plot()
