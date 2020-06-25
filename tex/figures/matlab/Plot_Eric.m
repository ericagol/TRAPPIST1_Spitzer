close all;
% clear all;
clc;
% MR = load('/Users/cdorn/Documents/Output_publications/Papers/co-authored/Grimm_trappist/PDFs/MR.dat');
MR = load('/Users/cdorn/Documents/Output_publications/Papers/co-authored/Agol_Eric_paper2020/T1_mass_radius_posterior_1M.txt');
addpath('/Users/cdorn/Documents/Postdoc/DATABASE_exoplanets')
%M_star [M_sun], R_star [R_sun], M_b, M_c, M_d, M_e, M_f, M_g, M_h, R_b, R_c, R_d, R_e, R_f, R_g, R_h


                MRplot=openfig('/Users/cdorn/Documents/Postdoc/DATABASE_exoplanets/MR.fig');
                MRdata=get(gca,'Children');
                XData=get(MRdata,'XData');  %get the x data
                YData=get(MRdata,'YData');  %get the y data
                close all
                
XMR = [0*MR(:,3) MR(:,3) MR(:,3+7)];
XMR(end+1:end+length(MR),1:3) = [0*MR(:,3)+1 MR(:,4) MR(:,4+7)];
XMR(end+1:end+length(MR),1:3) = [0*MR(:,3)+2 MR(:,5) MR(:,5+7)];
XMR(end+1:end+length(MR),1:3) = [0*MR(:,3)+3 MR(:,6) MR(:,6+7)];
XMR(end+1:end+length(MR),1:3) = [0*MR(:,3)+4 MR(:,7) MR(:,7+7)];
XMR(end+1:end+length(MR),1:3) = [0*MR(:,3)+5 MR(:,8) MR(:,8+7)];
XMR(end+1:end+length(MR),1:3) = [0*MR(:,3)+6 MR(:,9) MR(:,9+7)];
MR  = XMR;
%% preparing figure
clf
figure(1), hold on;




set(gca,'FontSize',16,'XMinorTick','on','YMinorTick','on');

lhx = get(gca,'xlabel');
set(lhx,'FontSize',16);
lhy = get(gca,'ylabel');
set(lhy,'FontSize',16);

box on;


        xlabel('planet mass M [M_{\oplus}]');
        ylabel('planet radius R [R_{\oplus}]');
        
        mrnewmartin = load('mass-radius_relationship_interior_UNTERBORN17_5percent_water_steam.dat');
        plot(mrnewmartin(:,1),mrnewmartin(:,2),'color',[252 148 3]./254,'LineStyle','-','LineWidth',2)
 mrnewmartin2 = load('mass-radius_relationship_interior_UNTERBORN17_0.01percent_water_steam.dat');       
 plot(mrnewmartin2(:,1),mrnewmartin2(:,2),'color',[252 148 3]./254,'LineStyle','--','LineWidth',2)
mrf6 = load('massradius_djoeke_big.ddat');
% M_v6 = mrf6(mrf6(:,3)==0.25095511880408494,1);%M_v6 = sort(M_v6);
% R_v6 = mrf6(mrf6(:,3)==0.25095511880408494,2);%R_v6 = sort(R_v6);
% plot(M_v6,R_v6,'color',[0 0. 0.8],'LineStyle','--','LineWidth',2)
% 
% M_v6 = mrf6(mrf6(:,3)==0.14894884665530439,1);%M_v6 = sort(M_v6);
% R_v6 = mrf6(mrf6(:,3)==0.14894884665530439,2);%R_v6 = sort(R_v6);
% plot(M_v6,R_v6,'color',[0 0. 0.8],'LineStyle','-','LineWidth',2)

M_v6 = mrf6(mrf6(:,3)==5.2137643930878358E-002,1);%M_v6 = sort(M_v6);
R_v6 = mrf6(mrf6(:,3)==5.2137643930878358E-002,2);%R_v6 = sort(R_v6);
plot(M_v6,R_v6,'color',[53 153 181]./254,'LineStyle','-','LineWidth',2)

% mrf4 = load('/Users/cdorn/Documents/Postdoc/InverseProblem/program/ExoplanetMCMC/Xu_profiles/massradius_femg075_water05_ratio_new.ddat');
% M_v4 = mrf4(:,6);M_v4 = sort(M_v4);
% R_v4 = mrf4(:,7);R_v4 = sort(R_v4);
% plot(M_v4,R_v4,'color',[0 0. 0.8],'LineStyle','-','LineWidth',2)

mrf5 = load('/Users/cdorn/Documents/Postdoc/InverseProblem/program/ExoplanetMCMC/Xu_profiles/massradius_femg0_mgsi_1-0_ratio_new.ddat');
M_v5 = mrf5(:,1);M_v5 = sort(M_v5);
R_v5 = mrf5(:,2);R_v5 = sort(R_v5);
plot(M_v5,R_v5,'color',[0.1 0 0],'LineStyle','-.','LineWidth',2)

% mrf = load('/Users/cdorn/Documents/Postdoc/InverseProblem/program/ExoplanetMCMC/Xu_profiles/massradius_femg0-5_mgsi_1-02_new.ddat');
% M_v3 = mrf(:,1);M_v3 = sort(M_v3);
% R_v3 = mrf(:,2);R_v3 = sort(R_v3);
% plot(M_v3,R_v3,'color',[0.1 0 0],'LineStyle',':','LineWidth',2)

mrf2 = load('/Users/cdorn/Documents/Postdoc/InverseProblem/program/ExoplanetMCMC/Xu_profiles/massradius_femg0-75_mgsi1-02_new.ddat');
M_v2 = mrf2(:,1);M_v2 = sort(M_v2);
R_v2 = mrf2(:,2);R_v2 = sort(R_v2);
plot(M_v2,R_v2,'color',[0 0 0.1],'LineStyle','--','LineWidth',2);%0.9 0.3 0.8


M_v6 = mrf6(mrf6(:,3)==0.0,1);%M_v6 = sort(M_v6);
R_v6 = mrf6(mrf6(:,3)==0.0,2);%R_v6 = sort(R_v6);
plot(M_v6,R_v6,'color',[0 0.1 0],'LineStyle','-','LineWidth',2)
% mrc = load('/Users/cdorn/Documents/Postdoc/InverseProblem/program/ExoplanetMCMC/Xu_profiles/massradius_solarratio_new.ddat');
% M_v1 = mrc(:,1);M_v1 = sort(M_v1);
% R_v1 = mrc(:,2);R_v1 = sort(R_v1);
% plot(M_v1,R_v1,'color',[0 0.1 0],'LineStyle','-','LineWidth',2)


                plot(XData{4,:}',YData{4,:}','k.','LineWidth',2);

legend('Fe/Mg=0.75, Mg/Si=1.02, m_{steam}/M=0.05','Fe/Mg=0.75, Mg/Si=1.02, m_{steam}/M=10^{-4}','Fe/Mg=0.75, Mg/Si=1.02, m_{water}/M=0.05','Fe/Mg=0.0 , Mg/Si=1.02','Fe/Mg=0.75, Mg/Si=1.02 (suggested by U17)','Fe/Mg=0.83, Mg/Si=1.02 (solar abundance)','pure Fe','Location','SouthEast')
legend('boxoff')
legend('AutoUpdate','off')



uia = unique(MR(:,1));
%    mrf4 = load('/Users/cdorn/Documents/Postdoc/InverseProblem/program/ExoplanetMCMC/Xu_profiles/massradius_femg075_water05_ratio_new.ddat');
% M_v4 = mrf4(:,6);M_v4 = sort(M_v4);
% R_v4 = mrf4(:,7);R_v4 = sort(R_v4);
% plot(M_v4,R_v4,'color',[0 0. 0.8],'LineStyle','-','LineWidth',2)
%%
for uu = 1: length(uia);
   
%     fprintf('Planet %i: \n',uu)
% hold on
[yM, xM] = ecdf(MR(MR(:,1)==uia(uu),2));
[yR, xR] = ecdf(MR(MR(:,1)==uia(uu),3));
[yg, xg] = ecdf((MR(MR(:,1)==uia(uu),2))./(MR(MR(:,1)==uia(uu),3)).^3);
% [h, statsM] = cdfplot(MR(MR(:,1)==uia(uu),2));
% [h, statsR] = cdfplot(MR(MR(:,1)==uia(uu),3));
% 
medianM(uu) = xM(find(yM>0.5,1,'first'));
stdM_p =  (xM(find(yM>0.5+(0.6827/2),1,'first'))-medianM(uu));
stdM_m =  abs(medianM(uu) - xM(find(yM>0.5-(0.6827/2),1,'first')));
% 
% 
medianR(uu) = xR(find(yR>0.5,1,'first'));
stdR_p =  abs(xR(find(yR>0.5+(0.6827/2),1,'first'))-medianR(uu));
stdR_m =  abs(medianR(uu) - xR(find(yR>0.5-(0.6827/2),1,'first')));
% 
mediang = xg(find(yg>0.5,1,'first'));

stdG_p =  abs(xg(find(yg>0.5+(0.6827/2),1,'first'))-mediang);
stdG_m =  abs(mediang - xg(find(yg>0.5-(0.6827/2),1,'first')));

 
  cdo = corr(MR(MR(:,1)==uia(uu),2),MR(MR(:,1)==uia(uu),3));
%  
  
fprintf('Planet %i: &%.3f & %.3f & %.3f & %.3f & %.3f &%.3f & %.3f & %.3f &%.3f & %.3f \\\\ \n',uu,medianM(uu),stdM_m,stdM_p,medianR(uu),stdR_m,stdR_p,cdo,mediang,stdG_m,stdG_p)


end


%%
load('/Users/cdorn/Documents/Postdoc/DATABASE_exoplanets/data.mat')

icount = 63;
% TRAPPIST-1b
icount = icount+1;
data(icount).name='Trappist-1b';
data(icount).radius=1.121;%1.125; %
data(icount).Rerr=[0.032,0.031];%[0.065,0.06];
data(icount).mass=1.017;%1.0915;  
data(icount).Merr=[0.143,0.154];%[0.1923, 0.1863];
data(icount).temp = 400;
% TRAPPIST-1c
icount = icount+1;
data(icount).name='Trappist-1c';
data(icount).radius=1.095;%1.101;%
data(icount).Rerr=[0.031,0.03];%[0.064,0.059];%
data(icount).mass=1.156;%1.2720 ; 
data(icount).Merr=[0.131,0.142];%[ 0.2289,0.1983];
data(icount).temp = 341;
% TRAPPIST-1d
icount = icount+1;
data(icount).name='Trappist-1d';
data(icount).radius=0.784;%0.788;%
data(icount).Rerr=[0.023,0.023];%[0.046,0.044];%
data(icount).mass=0.297;%0.3178 ;
data(icount).Merr=[0.035 0.039];%[ 0.0719, 0.0700];
data(icount).temp = 288;
% TRAPPIST-1e
icount = icount+1;
data(icount).name='Trappist-1e';
data(icount).radius=0.91;%0.912;%
data(icount).Rerr=[0.027,0.026];%[0.053,0.048];%
data(icount).mass=0.772;%0.6748 ;
data(icount).Merr=[0.075 0.079];%[ 0.1025, 0.1065];
data(icount).temp = 251;
% TRAPPIST-1f
icount = icount+1;
data(icount).name='Trappist-1f';
data(icount).radius=1.046;%1.050;
data(icount).Rerr=[0.03,0.029];%[0.060,0.057];%
data(icount).mass=0.934;%0.8850; 
data(icount).Merr=[0.078,0.08];%;[0.1284, 0.1230];
data(icount).temp = 219;
% TRAPPIST-1g
icount = icount+1;
data(icount).name='Trappist-1g';
data(icount).radius=1.148;%1.153;%
data(icount).Rerr=[0.033,0.032];%[0.066,0.061];%
data(icount).mass=1.148;%1.1328 ; %1.34
data(icount).Merr=[0.095,0.098];%[ 0.1574, 0.1575];
data(icount).temp = 198;
% TRAPPIST-1h
icount = icount+1;
data(icount).name='TRAPPIST-1h';
data(icount).radius=0.773;%0.778;%
data(icount).Rerr=[0.027,0.026];%[0.049,0.041];%
data(icount).mass=0.331; 
data(icount).Merr=[0.049, 0.056];
data(icount).temp = 168; % K




    hold on

    
xloc = 1.4;
yloc = 1;

for uu=[2 3 4 5 6 7 1]
    clear N C

    temp = data(63+uu).temp;
mases = (MR(MR(:,1)==uia(uu),2));
radi  = (MR(MR(:,1)==uia(uu),3));
% edges = {min(mases):0.01:max(mases),min(radi):0.01:max(radi)};
[N,C] = hist3([mases,radi],[50 ,50]);
 evalin('base',sprintf('luluu_%i = imagesc(C{1},C{2},N''*0+temp);',uu));
      
% luluu_ = imagesc(C{1},C{2},N'*0+temp);
 Ttem(uu) = temp;
N = N./max(max(N));
% luluu.AlphaData = N';
 evalin('base',sprintf('luluu_%i.AlphaData = N'';',uu));
cornyplot([mases, radi])
%     idcc = ((CC(:,2)==uia(uu))&(CC(:,1)==1));
%     plot((CC(idcc,3)),CC(idcc,4),'color',[0.2 .2 .2])
    
%     lx = [1.27 1.27 1.37 1.37 1.27];
%     ly = [yloc yloc+0.012 yloc+0.012 yloc yloc]- 0.02*uu +0.02;
%     lhr = fill(lx,ly,temp*ones(5,1),'EdgeColor','none');
%     plot(lx,ly,'color',[0.3 .3 .3])

end 

 
         

%           scatter(medianM,medianR,50,Ttem,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',1.5)
            scatter(medianM,medianR,60,Ttem,'filled')
            scatter(medianM,medianR,20,'k','filled')
  
% xloc = 1.14;
% yloc = 0.734;



sz = 20;
yloc = yloc + 0.005;
t=text(1.48,1.07,'b');
s=t.FontSize;
t.FontSize=sz;
t=text(1.3,1.05,'c');
s=t.FontSize;
t.FontSize=sz;
t=text(0.41,0.75,'d');
s=t.FontSize;
t.FontSize=sz;
t=text(0.7,0.88,'e');
s=t.FontSize;
t.FontSize=sz;
t=text(1.07,0.996,'f');
s=t.FontSize;
t.FontSize=sz;
t=text(1.35,1.2,'g');
s=t.FontSize;
t.FontSize=sz;
t=text(0.35,0.7,'h');
s=t.FontSize;
t.FontSize=sz;
hold on
%%SOLAR SYSTEM

M   = [0.0553;0.815;1;0.107 ;14.536;17.147];
R   = [0.383 ;0.950;1;0.532 ; 3.981; 3.865];
rho = [0.984 ;0.951;1;0.7134; 0.230; 0.297];
temp= [449   ;328  ;279;226 ;64   ;51   ];

sz = 15;
rhoEarth  = 5.514; % ref. density
massEarth = 5.97236473e24; %kg
radEarth  = 6378100;%m
        scatter(M,R,72,'k','filled');
        scatter(M,R,60*ones(length(M),1),temp,'filled');
                        t=text((M(3)+0.01) ,(R(3)-0.01),'{Earth}');
                s=t.FontSize;
                t.FontSize=sz;
                t=text((M(2)+0.01) ,(R(2)-0.01),'{Venus}');
                s=t.FontSize;
                t.FontSize=sz;
                                t=text((M(4)+0.01) ,(R(4)-0.01),'{Mars}');
                s=t.FontSize;
                t.FontSize=sz;
% set(gca, 'XScale', 'log')
ax = gca;
xlim([0.1 1.6])
ylim([0.4 1.3])
% xlim([0.15 1.6])
ax.XTick = ([0.2 0.4 0.6 0.8 1. 1.2 1.4 1.6 2.0]);
% ylim([0.6 1.3])
% set(gcf,'Position',[1     1   786   628]);
 set(gcf, 'PaperPosition', [0 0 30 25]);
cbar =  colorbar;
ylabel(cbar,'equilibrium temperature [K]','FontSize',20)
 caxis([150 400])
 colormap(parula(8))
  print('-dpng','-r300',sprintf('Figure_MR.png'))
% _______________

% plot(MR(:,3),MR(:,3+7),'.')
% plot(MR(:,4),MR(:,4+7),'.')
% plot(MR(:,5),MR(:,5+7),'.')
% plot(MR(:,6),MR(:,6+7),'.')
% plot(MR(:,7),MR(:,7+7),'.')
% plot(MR(:,8),MR(:,8+7),'.')
% plot(MR(:,9),MR(:,9+7),'.')