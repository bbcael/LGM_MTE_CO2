clear all; close all; clc; % clear workspace

%% STEP 1 estimate ef(T) -- Step 1A: phytoplankton

nb = 200; % bootstrap iteration number: enough for significant figures in calculations/paper

% load data from Browning & Moore 2023
opts = delimitedTextImportOptions("NumVariables", 6);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["Var1", "Tp", "Var3", "Var4", "Gp", "Var6"];
opts.SelectedVariableNames = ["Tp", "Gp"];
opts.VariableTypes = ["string", "double", "string", "string", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, ["Var1", "Var3", "Var4"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var3", "Var4"], "EmptyFieldRule", "auto");
tbl = readtable("/Users/ep_cael/Documents/work/1priorities/lgm/Bioassay_Figure4_data.csv", opts);
Tp = tbl.Tp; % T = temperature, p = phytoplankton
Gp = tbl.Gp; % G = growth rate
Gp(Gp<1e-4) = NaN; % remove negative values
i = find(~isnan(Tp)); % remove missing values
Tp = Tp(i); Gp = Gp(i);
i = find(~isnan(Gp));
Tp = Tp(i); Gp = Gp(i);
clear i opts tbl;

maxpct = 100.*round(1-1./(length(Gp)./5),2); % largest safe %ile to estimate via quantile regression: 99th is ok
Tp = 1./(8.62e-5.*(Tp+273)); % convert temperature to relevant units

% quantile regression for 90th to maxpct %iles for robustness
for i = 1:(maxpct-89);
    pct = 89+i;
    [p,stats] = quantreg(Tp,log(Gp),pct./100,1,nb); % perform quantile regression
    Ep_90toMaxPct(i) = -p(1); % temperature dependence
    Ip_90toMaxPct(i) = p(2); % intercept
    Ep_unc(i) = stats.pse(1); % Ep uncertainty
    Ip_unc(i,:) = stats.pboot(2); % Ip uncertainty
end
Ep_90toMaxPct_mean_std = [mean([Ep_90toMaxPct]) std([Ep_90toMaxPct])]
Ep_95_unc = Ep_unc(6)

clear Ep_90toMaxPct_mean_std i j p stats maxpct pct Ep_unc Ep_95_unc Ip_unc;

%% figure 1a

%{
scatter(Tp,Gp,50,[.5 0 .5],'filled'); 
hold on;
plot(linspace(38,43),exp(-Ep_90toMaxPct(6).*linspace(38,43)+Ip_90toMaxPct(6)),'k','linewidth',3);
for i = 1:length(Ep_90toMaxPct);
plot(linspace(38,43),exp(-Ep_90toMaxPct(i).*linspace(38,43)+Ip_90toMaxPct(i)),'--k');
end
lgnd = legend('Data, Browning \& Moore 2023','95th Percentile, $E_p = 0.34\pm0.05$ eV','90-99th Percentiles, $E_p = 0.31\pm0.03$ eV')
set(lgnd,'interpreter','latex','fontsize',15,'location','northeast')
set(gca,'ticklabelinterpreter','latex','fontsize',15)
box on;
xlabel('Temperature 1/kT [eV$^{-1}$]','interpreter','latex')
ylabel('Phytoplankton Growth Rate [d$^{-1}$]','interpreter','latex')
axis([38.1 42.9 0 2.5])
set(gca,'xtick',[39 40 41 42],'ytick',[0.5 1 1.5 2 2.5])
text(38.9,2.4,'24$^\circ$C','interpreter','latex','fontsize',15)
text(39.9,2.4,'17$^\circ$C','interpreter','latex','fontsize',15)
text(40.9,2.4,'10$^\circ$C','interpreter','latex','fontsize',15)
text(41.95,2.4,'3$^\circ$C','interpreter','latex','fontsize',15)
%}

clear i Tp Gp lgnd Ep_90toMaxPct Ip_90toMaxPct;

%% STEP 1 estimate ef(T) -- Step 1B: zooplankton

load huntley_lopez_appendix_A.mat; % load copepod growth rates from Huntley & Lopez 1992, z = zooplankton
Gz = growthrate'; Tz = temperature'; clear growthrate temperature;
x0 = [25 0.65]; % initial fit values -- don't make a difference
Tz = 1./(8.62e-5.*(Tz+273)); % convert to relevant unts
mte = fittype(@(a,b,x) exp(a)*exp(-b*x)); % model to fit to data
for i = 1:nb;  % estimate parameters: bootstrap for uncertainty because fitting nonlinear model
        j = randi(130,1,130); % bootstrap resample
        [mdl,gof] = fit(Tz(j),Gz(j),mte,'StartPoint',x0); % estimate parameters
        coeffs = coeffvalues(mdl); Iz(i) = coeffs(1); Ez(i) = coeffs(2); % parameter estimates
end
Ez_mean_std = [mean(Ez) std(Ez)] % print mean & standard error of temperature dependence

clear Ez_mean_std coeffs gof i j mdl mte x0;

%% figure 1b

%{
for i = 1:length(Ez); % get uncertainty in fit
    Z(i,:) = exp(Iz(i))*exp(-Ez(i).*linspace(38,43));
end
ZU = prctile(Z,84.13); ZL = prctile(Z,15.87); % +/-1sigma error bars
clear i Z;
a = area(linspace(38,43)',[ZL; (ZU-ZL)]'); % plot as area
a(1).FaceColor = [1 1 1];
a(2).FaceColor = [0 .5 .5];
a(2).FaceAlpha = 0.2;
hold on;
b = scatter(Tz,Gz,50,[0 .5 .5],'filled'); 
hold on;
c = plot(linspace(38,43),exp(-linspace(38,43).*mean(Ez)+mean(Iz)),'k','linewidth',3);
lgnd = legend([b c],'Data, Huntley \& Lopez 1992','Mean, $E_z = 0.69\pm0.03$ eV')
set(lgnd,'interpreter','latex','fontsize',15,'location','northeast')
set(gca,'ticklabelinterpreter','latex','fontsize',15)
box on;
xlabel('Temperature 1/kT [eV$^{-1}$]','interpreter','latex')
ylabel('Phytoplankton Growth Rate [d$^{-1}$]','interpreter','latex')
axis([38.1 42.9 0 1.2])
set(gca,'xtick',[39 40 41 42],'ytick',[0.25 0.5 0.75 1 1.25])
text(38.9,1.15,'24$^\circ$C','interpreter','latex','fontsize',15)
text(39.9,1.15,'17$^\circ$C','interpreter','latex','fontsize',15)
text(40.9,1.15,'10$^\circ$C','interpreter','latex','fontsize',15)
text(41.95,1.15,'3$^\circ$C','interpreter','latex','fontsize',15)
xlabel('Temperature 1/kT [eV$^{-1}$]','interpreter','latex')
ylabel('Zooplankton Growth Rate [d$^{-1}$]','interpreter','latex')
%}

clear a b c i Gz Tz Ez Iz lgnd ZL ZU;

%% STEP 1 estimate ef(T) -- Step 1C: export efficiency

load cael_follows_export_efficiency.mat; % load data from Cael & Follows 2016
load maiti_export_efficiency.mat; % load additional data from Maiti et  al. 2013
ef = [EFb; EFh; EFl1; EFl2; EFm]; % combine five data compilations
Tef = [Tb; Th; Tl1; Tl2; Tm];
clear EFb EFh EFl1 EFl2 EFm Tb Th Tl1 Tl2 Tm;

a = 12.35:.001:12.65; % find alpha value that explains ef's T-dependence -- this range ok in practice, n.b. no estimates at edge values
for i = 1:nb; % bootstrap iterations
    for j = 1:length(a);
        b = randi(length(ef),1,length(ef)); % bootstrap resample
        x = Tef(b); y = ef(b)./(1-exp(a(j)).*exp(-.33./(8.62e-5.*(x+273)))); % normalize ef
        [rs(j),ps(j)] = corr(x,y,'type','spearman');
        [ss(j),qs(j)] = corr(x,abs(y-mean(y)),'type','spearman');
    end
    [~,indx] = min(abs(rs)+abs(ss));
    RS(i) = rs(indx); % correlations of normalized ef w/ T
    PS(i) = ps(indx); % p-value thereof
    QS(i) = qs(indx); % correlations of normalized ef spread w/ T
    SS(i) = ss(indx); % p-value thereof
    A(i) = a(indx); % estimates of alpha
    i
    Y(i) = sum(y>1)./length(y);
end

corrs_pvals = [mean(RS) mean(PS) mean(SS) mean(QS)] % to report for figure 3
alpha_mean_std = [mean(A) std(A)] % to use later
beta = mean(ef./(1-exp((A)).*exp(-.33./(8.62e-5.*(Tef+273))))); % to use later -- makes no difference for results though
[mean(beta) std(beta)]

clear RS PS SS QS a i j b x y indx rp pp rs ps sp qp ss qs A beta ans Y;

% n.b. for source of uncertainty #1 – alpha – re-run calculations with alpha +/-1 std. err. 

%% figure 1c

Tef = 1./(8.62e-5.*(Tef+273)); % convert to relevant unts
%{
efplot = ef;
efplot(efplot>1) = 0.999; % plot ef>1 values in [0,1] axis
EF0 = 1-exp(12.48).*exp(-.33.*linspace(38,43)); % mean, upper, & lower 
EFL = 1-exp(12.44).*exp(-.33.*linspace(38,43));
EFU = 1-exp(12.52).*exp(-.33.*linspace(38,43))-EFL;
a = area(linspace(38,43)',[EFL; EFU]'); % plot +/-1sigma as error bar as in figure 2
a(1).FaceColor = [1 1 1];
a(2).FaceColor = [.5 .5 0];
a(2).FaceAlpha = 0.2;
hold on;
b = plot(linspace(38,43),EF0,'k','linewidth',3);
c = scatter(Tef,efplot,50,[.5 .5 0],'filled'); 
lgnd = legend([c b],'Data, Multiple Sources','$\phi(T), \alpha = 12.48\pm0.04$')
set(lgnd,'interpreter','latex','fontsize',15,'location','northwest')
set(gca,'ticklabelinterpreter','latex','fontsize',15)
box on;
axis([38.1 42.9 0 1])
text(38.2,.78,'$\rho (ef/\phi,T) = -0.01$, $p \gg 0.1$','interpreter','latex','fontsize',16)
text(38.2,.71,'$\rho (\vert ef/\phi - \langle ef/\phi \rangle \vert,T) = 0.02$, $p \gg 0.1$','interpreter','latex','fontsize',16)
xlabel('Temperature 1/kT [eV$^{-1}$]','interpreter','latex')
ylabel('Export Efficiency [$\sim$]','interpreter','latex')
set(gca,'xtick',[39 40 41 42],'ytick',[0.2 0.4 0.6 0.8 1])
text(38.9,.96,'24$^\circ$C','interpreter','latex','fontsize',15)
text(39.9,.96,'17$^\circ$C','interpreter','latex','fontsize',15)
text(40.9,.96,'10$^\circ$C','interpreter','latex','fontsize',15)
text(41.95,.96,'3$^\circ$C','interpreter','latex','fontsize',15)
%}

clear a alpha_mean_std b c corrs_pvals EF0 EFL EFU lgnd nb efplot;

%% figure 1d

%{
beta = ef./(1-exp(12.48).*exp(-.33.*Tef)); % normalize data
binedges = 38.3:.25:42.8;
for i = 1:(length(binedges)-1);
    tb(i) = (binedges(i)+binedges(i+1))./2;
    eb(i) = mean(beta(Tef<binedges(i+1) & Tef>binedges(i)));
    nb(i) = length(beta(Tef<binedges(i+1) & Tef>binedges(i)));
end
beta(beta>2.05) = 2.04; % fit one outlier on plot
c = scatter(Tef,beta,50,[.7 .4 .05],'filled'); % plot
hold on;
min(nb)
scatter(tb,eb,200,'k','filled');
axis([38.1 42.9 0 2.05]);
box on;
lgnd = legend('Data, Multiple Sources','$T$-binned Mean')
set(lgnd,'interpreter','latex','fontsize',15,'location','northwest')
set(gca,'ticklabelinterpreter','latex','fontsize',15)
set(gca,'xtick',[39 40 41 42],'ytick',[0 0.5 1 1.5 2])
text(38.9,1.95,'24$^\circ$C','interpreter','latex','fontsize',15)
text(39.9,1.95,'17$^\circ$C','interpreter','latex','fontsize',15)
text(40.9,1.95,'10$^\circ$C','interpreter','latex','fontsize',15)
text(41.95,1.95,'3$^\circ$C','interpreter','latex','fontsize',15)
xlabel('Temperature 1/kT [eV$^{-1}$]','interpreter','latex')
ylabel('$\phi(T)$-normalized Export Efficiency [$\sim$]','interpreter','latex')
%}
clear tb nb eb binedges beta ans c i ef Tef lgnd;

%% STEP 2 estimate ∆ef(∆T) -- Step 2A: load data – 2 NPP, 3 SSTmod, & 4 SSTlgm products, & convert to ef

load clim4cael.mat; % P = net primary production
P1 = cafes; P1 = permute(P1,[2 1 3]); P1 = flip(P1,2); % CAFE NPP climatology; orient to match SSTlgm
p1(181:360,:,:) = P1(1:180,:,:); p1(1:180,:,:) = P1(181:360,:,:); P1 = p1; clear p1; % rotate to match SSTlgm
P2 = cbpm; P2 = permute(P2,[2 1 3]); P2 = flip(P2,2); % CbPMv2 climatology; orient to match SSTlgm
p2(181:360,:,:) = P2(1:180,:,:); p2(1:180,:,:) = P2(181:360,:,:); P2 = p2; clear p2; % rotate to match SSTlgm
clear bbp cafes cbpm fmicro logChl z_eu; % clear other extra data
Tm1 = ncread('sst.mon.ltm.nc','sst'); Tm1 = flip(Tm1,2); % COBE2; Tm = modern SST; orient to match SSTlgm
ym1 = ncread('sst.mon.ltm.nc','lat'); ym1 = repmat(ym1,1,360)'; % get lat grid
xm1 = ncread('sst.mon.ltm.nc','lon'); xm1 = repmat(xm1,1,180); ym1 = fliplr(ym1); % get lon grid
cd ersst_v5  % ERSSTv5
Tm2 = ncread('sst.mon.ltm.nc','sst'); Tm2 = flip(Tm2,2); % same as above
ym2 = -(180./89).*(-44:44)'; ym2 = repmat(ym2,1,180)'; ym2 = fliplr(ym2); % same as above
xm2 = ncread('sst.mon.ltm.nc','lon'); xm2 = repmat(xm2,1,89); xm2 = xm2+1; % same as above
cd ..
Tm2(Tm2<-10) = NaN; Tm1(Tm1<-10) = NaN; % remove bad values
Tm0 = ncread('baseline_LGMRholo_sstice_clim.nc','SST_cpl'); % Tm from Osman et al. 2021 reconstruction
y = ncread('baseline_LGMRholo_sstice_clim.nc','lat'); y = repmat(y,1,144)'; % get lat grid
x = ncread('baseline_LGMRholo_sstice_clim.nc','lon'); x = repmat(x,1,96); % get lon grid

% for each month, interpolate P1, P2, Tm1, & Tm2 onto Tm0 grid
for i = 1:12;
    [i 1]
    prout = P1(:,:,i);
    F = scatteredInterpolant(double(xm1(:)),double(ym1(:)),double(prout(:)));
    prt(:,:,i) = F(double(x),double(y));
end
P1 = prt;
clear i prout F prt;
for i = 1:12;
    [i 2]
    prout = P2(:,:,i);
    F = scatteredInterpolant(double(xm1(:)),double(ym1(:)),double(prout(:)));
    prt(:,:,i) = F(double(x),double(y));
end
P2 = prt;
clear i prout F prt;
for i = 1:12;
    [i 3]
    prout = Tm1(:,:,i);
    F = scatteredInterpolant(double(xm1(:)),double(ym1(:)),double(prout(:)));
    prt(:,:,i) = F(double(x),double(y));
end
Tm1 = prt;
clear i prout F prt;
for i = 1:12;
    [i 4]
    prout = Tm2(:,:,i);
    F = scatteredInterpolant(double(xm2(:)),double(ym2(:)),double(prout(:)));
    prt(:,:,i) = F(double(x),double(y));
end
Tm2 = prt;
clear i prout F prt xm1 xm2 ym1 ym2 ans;

TgO = ncread('LGMRlgm_sstice_clim.nc','SST_cpl'); % g = glacial, O = Osman SSTlgm
TgJ = ncread('annan_sstice_clim.nc','SST_cpl'); % J = J. Annan SSTlgm
TgD = ncread('amrhein_sstice_clim.nc','SST_cpl'); % D = D. Amrhein SSTlgm
TgT = ncread('lgmDA_sstice_clim.nc','SST_cpl'); % T = Tierney SSTlgm

% remove non-shared grid cells -- irrelevant coastal ocean bits related to gridding
Tm0(isnan(Tm1)) = NaN; Tm1(isnan(Tm0)) = NaN; Tm2(isnan(Tm0)) = NaN; Tm2(isnan(Tm1)) = NaN;
Tm1(isnan(Tm2)) = NaN; Tm0(isnan(Tm2)) = NaN; P1(isnan(Tm0)) = NaN; P2(isnan(Tm0)) = NaN;
TgT(isnan(Tm0)) = NaN; TgD(isnan(Tm0)) = NaN; TgJ(isnan(Tm0)) = NaN; TgO(isnan(Tm0)) = NaN;

% clip NPP outliers -- doesn't affect final results but these values are likely erroneous so should do anyways
P1(P1>prctile(P1(:),99)) = prctile(P1(:),99); 
P2(P2>prctile(P2(:),99)) = prctile(P2(:),99);

% convert T to ef
efD = .37.*(1-exp(12.48).*exp(-.33./(8.62e-5.*(273+TgD))));
efJ = .37.*(1-exp(12.48).*exp(-.33./(8.62e-5.*(273+TgJ))));
efO = .37.*(1-exp(12.48).*exp(-.33./(8.62e-5.*(273+TgO))));
efT = .37.*(1-exp(12.48).*exp(-.33./(8.62e-5.*(273+TgT))));
ef1 = .37.*(1-exp(12.48).*exp(-.33./(8.62e-5.*(273+Tm1))));
ef2 = .37.*(1-exp(12.48).*exp(-.33./(8.62e-5.*(273+Tm2))));
ef0 = .37.*(1-exp(12.48).*exp(-.33./(8.62e-5.*(273+Tm0))));

%% STEP 2 estimate ∆ef(∆T) -- Step 2B: mask out areas where nutrients don't support ∆ef change

yNi = ncread('woa23_all_n00_01.nc','lat'); % load world ocean atlas nitrate
xNi = ncread('woa23_all_n00_01.nc','lon');
xNi = mod(xNi,360); % same longitude notation as Tm0
NO3i = ncread('woa23_all_n00_01.nc','n_an');
NO3i = NO3i(:,:,1); % WOA annual mean surface nitrate
yNii = ncread('GLODAPv2.2016b.NO3.nc','lat');
xNii = ncread('GLODAPv2.2016b.NO3.nc','lon');
xNii = mod(xNii,360); % same longitude units as Tm0
NO3ii = ncread('GLODAPv2.2016b.NO3.nc','NO3');
NO3ii = NO3ii(:,:,1); % GLODAP annual mean surface nitrate
NO3i(NO3i>.57347 & NO3i<.57349) = NaN; % remove an anomalous value found at 360 points from distribution
xgrad = 0:.01:max(NO3i(:)); % n.b. insensitive to step size: = 0.70 for 0.002 (i.e. 1/5) or 0.70 for 0.05 (i.e. 5x)
[ypdf,xpdf] = ksdensity(NO3i(:),xgrad); % PDF of NO3
ypdf = interp1(xpdf,ypdf,xgrad,'spline'); % smooth via spline for better gradient calculation
dy=gradient(ypdf,xgrad(2)-xgrad(1)); % calculate first derivative
ddy=gradient(dy,xgrad(2)-xgrad(1)); % calculate second derivative
[~,indx] = max(ddy); % find maximum second derivative
Ni_thresh = xgrad(indx); % define threshold like this -- double or halve for sensitivity
xgrad = 0:.01:max(NO3ii(:)); % n.b. insensitive to step size: = 0.59 for 0.002 (i.e. 1/5) or 0.60 for 0.05 (i.e. 5x)
[ypdf,xpdf] = ksdensity(NO3ii(:),xgrad); % same as above
ypdf = interp1(xpdf,ypdf,xgrad,'spline');
dy=gradient(ypdf,xgrad(2)-xgrad(1));
ddy=gradient(dy,xgrad(2)-xgrad(1));
[~,indx] = max(ddy);
Nii_thresh = xgrad(indx);
clear xgrad ypdf xpdf dy ddy indx;

yNi = repmat(yNi,1,360)'; yNii = repmat(yNii,1,360)'; xNi = repmat(xNi,1,180); xNii = repmat(xNii,1,180); % make grid
F = scatteredInterpolant(double(xNi(:)),double(yNi(:)),NO3i(:)); % interpolate to common grid
NO3i = F(x,y);
F = scatteredInterpolant(double(xNii(:)),double(yNii(:)),NO3ii(:));
NO3ii = F(x,y);
clear F yNi yNii xNi xNii;
Nmask1 = NO3i; q = find(NO3i<Ni_thresh); Nmask1(q) = 0; q = find(NO3i>Ni_thresh); Nmask1(q) = 1; % define masks
Nmask2 = NO3ii; q = find(NO3ii<Nii_thresh); Nmask2(q) = 0; q = find(NO3ii>Nii_thresh); Nmask2(q) = 1; % define masks
clear NO3i NO3ii Ni_thresh Nii_thresh q;

%% figure 2a

%{

efm = 1/3.*(ef1 + ef2 + ef0); % average of three SSTmod ef maps
efg = 1./4.*(efD + efJ + efO + efT); % average of four SSTlgm ef maps
def = efg-efm; % difference
def = nanmean(def,3); % annual mean
x(x>356) = 360; % visual tweak for map continuity
def(Nmask1==0) = NaN; % remove nitrate-masked areas

figure;
axesm robinson
axesm('robinson','MapLatLimit',[-90 90],'MapLonLimit',[20 380])
framem off
axis off
h = contourfm(y,x,double(def),[0:.01:.07]);
hold on;
colormap(turbo)
c = colorbar;
set(c,'fontsize',24,'ytick',[0:.01:.07],'ticklabelinterpreter','latex')
caxis([0 .07])
hold on;
landareas = shaperead('landareas.shp','UseGeoCoords',true);
geoshow(landareas,'FaceColor',[.9 .85 .8]);
title('$\langle ef(\mathrm{LGM}) \rangle - \langle ef(\mathrm{PI})\rangle$','fontsize',24,'interpreter','latex')
%}

clear efm efg def h c landareas; 

%% figure 2b

%{
figure;
ttl = {'Osman et al.','Tierney et al.','Amrhein et al.','Annan et al.'};
efm = 1/3.*(ef1 + ef2 + ef0); % average of three SSTmod ef maps
for i = 1:4; % repeat above but for individual reconstructions
    clear efg;
    if i==1; efg = efO; def = efg-efm; def = nanmean(def,3);
    elseif i == 2; efg = efT; def = efg-efm; def = nanmean(def,3);
    elseif i == 3; efg = efD; def = efg-efm; def = nanmean(def,3);
    elseif i == 4; efg = efJ; def = efg-efm; def = nanmean(def,3);
    end
    x(x>356) = 360; def(Nmask1==0) = NaN;
    subplot(2,2,i);
    axesm robinson
    axesm('robinson','MapLatLimit',[-90 90],'MapLonLimit',[20 380])
    framem off
    axis off
    h = contourfm(y,x,double(def),[0:.01:.07]);
    hold on;
    colormap(turbo)
    caxis([0 .07])
    landareas = shaperead('landareas.shp','UseGeoCoords',true);
    geoshow(landareas,'FaceColor',[.9 .85 .8]);
    title(ttl(i),'fontsize',24,'interpreter','latex')
end
%}

clear efm efg def x h c landareas ttl; 

%% STEP 2 estimate ∆ef(∆T) -- Step 2C: get global ∆ef

y = cosd(y); % latitude weighting
y = repmat(y,1,1,12); % make same size as other fields
Nmask1 = repmat(Nmask1,1,1,12);
Nmask2 = repmat(Nmask2,1,1,12);

dEF = zeros(4,3,2,2); % 4x SSTlgm, 3x SSTmod, 2x NPP, 2x NO3
dT = zeros(4,3,2,2); % associated GMSST changes
for i = 1:4; % just doing it as a nasty for loop for clarity/simplicity
    if i==1;
    ef_lgm = efD;
    SST_lgm = TgD;
    elseif i==2; 
    ef_lgm = efJ;
    SST_lgm = TgJ;
    elseif i==3; 
    ef_lgm = efO;
    SST_lgm = TgO;
    elseif i==4;
    ef_lgm = efT;
    SST_lgm = TgT;
    end
    for j = 1:3;
        if j==1;
        ef_mod = ef0;
        SST_mod = Tm0;
        elseif j==2;
        ef_mod = ef1;
        SST_mod = Tm1;
        elseif j==3;
        ef_mod = ef2;
        SST_mod = Tm2;
        end
        for k = 1:2;
            if k==1;
            npp = P1;
            elseif k==2;
            npp = P2;
            end
            for l = 1:2;
                if l==1;
                no3 = Nmask1;
                elseif l==2;
                no3 == Nmask2;
                end
            
                landmask = ones(size(ef_lgm)); landmask(isnan(ef_lgm)) = NaN;
                ef_lgm(no3==0) = ef_mod(no3==0); % ef doesn't increase if nutrients don't allow it
                EF_lgm = nansum(landmask(:).*y(:).*npp(:).*ef_lgm(:)./nansum(landmask(:).*y(:).*npp(:)));
                EF_mod = nansum(landmask(:).*y(:).*npp(:).*ef_mod(:)./nansum(landmask(:).*y(:).*npp(:)));
                GMSST_lgm = nansum(landmask(:).*y(:).*SST_lgm(:)./nansum(landmask(:).*y(:)));
                GMSST_mod = nansum(landmask(:).*y(:).*SST_mod(:)./nansum(landmask(:).*y(:)));
                dEF(i,j,k,l) = (EF_lgm-EF_mod)./EF_mod; % relative change in global export efficiency
                dT(i,j,k,l) = GMSST_lgm-GMSST_mod; % GMSST change
                EFmod(i,j,k,l) = EF_mod; % global modern export efficiency -- need for step 3

            end
        end
    end
end
EFmod_mean_std = [mean(EFmod(:)) std(EFmod(:))]
clear i j k l landmask EF_lgm EF_mod GMSST_lgm GMSST_mod no3 npp SST_mod SST_lgm ef_lgm ef_mod ef* Nmask* P* Tg* Tm* y x EFmod_mean_std EFmod ans;

%% STEP 3 estimate ∆CO2(∆ef) -- Step 3A: version 1 of parameter estimation

dE = [-.558 .745 .845 .124 .205 .201 .111 -.017 -.011 -.127 .188 .329 -.025 .1 .124 -.330 .454 .653 -.236 .230 0]./2.35; % relative chages in export efficiency from Lauderdale & Cael 2021 GRL
dC = [-63 70 93 12 11 13 -2 -1 -2 -18 10 29 -8 6 9 -40 36 66 -22 25 0]./269; % relative changes in CO2 from same

nb = 200; % bootstrap iteration
for i = 1:nb;
    i
    j = randi(length(dE),1,length(dE));
    dCb = dC(j);
    dEb = dE(j);
    mdl1 = fitlm(dEb,dCb,'Intercept',false); % model II geometric mean w/o intercept
    mdl2 = fitlm(dCb,dEb,'Intercept',false);
    g1 = table2array(mdl1.Coefficients);
    g1 = g1(1);
    g2 = table2array(mdl2.Coefficients);
    g2 = g2(1);
    gamma(i) = sqrt(g1./g2); % distribution of gamma values
end
[mean(gamma) std(gamma)]
 
clear nb i j dCb dEb mdl1 mdl2 g1 g2 gamma ans;

%% figure 4

%{
c = plot(linspace(-.4,.4),-.88.*linspace(-.4,.4)+.019,'k','linewidth',3);
hold on;
b = scatter(dE,-dC,100,[.7 .4 .05],'filled'); 
lgnd = legend([b c],'Model Experiments, Lauderdale \& Cael 2021','$\gamma = 0.88 \pm 0.04$')
set(lgnd,'interpreter','latex','fontsize',15,'location','northeast')
set(gca,'ticklabelinterpreter','latex','fontsize',15)
box on;
xlabel('Relative change in export efficiency','interpreter','latex')
ylabel('Relative change in atmospheric CO$_2$','interpreter','latex')
axis([-.3 .4 -.4 .3])
clear c b lgnd;
%}
clear dC dE;

%% STEP 3 estimate ∆CO2(∆ef) -- Step 3B: version 2 of parameter estimation

% Cael et al. 2017 L&O:L box model, gamma = ef*P*V/(Ib*Psi), use estimates of each quantity
EF = .14; % from above, n.b. similar to or in exact agreement with other estimates
P = 53; % Johnson et al. '21
V = 1.33; % Charette et al. '10
Ib = 3300; % Goodwin & Cael '21
Psi = 108; % Liang et al. '17
V = V.*1e18; Psi = Psi.*1e6; P = P./31557600; % convert to relevant units
gamma_alt = (EF.*P.*V)./(Ib.*Psi)
clear EF P V Ib Psi gamma_alt

%% STEP 4 estimate CO2 change -- Step 4A: compile all sources of uncertainty 

% number of estimates of CO2 change: 4x SSTlgm by 3x SSTmod by 2x NPP by 2x NO3 by 3x alpha best/high/low by 3x beta best/high/low = 432

alpha_unc = 0.11; % relative uncertainty from alpha, estimated by recalculating dEF above with +/-1sigma alpha estimates
gamma = 0.88; gamma_unc = 0.05; % coefficient of proportionality & relative uncertainty from above
dEF = dEF(:); dT = dT(:); % other four sources of uncertainty accounted for in these estimates

dEF = repmat(dEF,3,1); dT = repmat(dT,3,1); % fold in gamma uncertainty
dEF(1:48) = gamma.*dEF(1:48).*278; dEF(49:96) = (gamma+gamma_unc).*dEF(49:96).*278; dEF(97:144) = (gamma-gamma_unc).*dEF(97:144).*278; % gamma uncertainty
dEF = repmat(dEF,3,1); dT = repmat(dT,3,1); % fold in alpha uncertainty
dEF(145:288) = dEF(145:288).*(1+alpha_unc); dEF(289:432) = dEF(289:432).*(1-alpha_unc); % alpha uncertainty

dCO2 = dEF; % CO2 change distribution
dCO2dT = -dEF./dT; % CO2/GMSST change distribution

clear alpha_unc gamma gamma_unc gamma gamma_unc dT dEF;

%% figure 3a

%{
figure;
h = histogram(dCO2dT,[5.54 6 6.5 7 7.5 8 8.5 9.05],'Normalization','PDF');
set(gca,'ticklabelinterpreter','latex','fontsize',16,'ytick',[],'xtick',[4:8]);
ylabel('Probability Density','interpreter','latex');
axis([5.44 9.15 0 .47])
xlabel('CO$_2$ change per GMSST change [ppm/K]','interpreter','latex')
box on;
h.FaceColor = [.7 .6 .1];
h.LineWidth = 1.5;
hold on;
g = plot(linspace(5,11),.9*normpdf(linspace(5,11),mean(dCO2dT),std(dCO2dT)),'--k','linewidth',3)
lgnd = legend('Estimates','7.2$\pm$0.8 ppm/K');
set(lgnd,'fontsize',16,'interpreter','latex','location','northeast')
%}

[mean(dCO2dT) std(dCO2dT)]

clear ans dC dCO2dT dE g h lgnd;

%% figure 3b

%{
figure;
h = histogram(dCO2,linspace(12.7,37.9,9),'Normalization','PDF');
set(gca,'ticklabelinterpreter','latex','fontsize',16,'ytick',[],'xtick',[10:5:40]);
ylabel('Probability density','interpreter','latex');
axis([12.4 38.1 0 .08])
xlabel('CO$_2$ change [ppm]','interpreter','latex')
box on;
h.FaceColor = [.7 .4 .05];
h.LineWidth = 1.5;
hold on;
lgnd = legend('Estimates, 23$\pm$6 ppm');
set(lgnd,'fontsize',16,'interpreter','latex','location','northeast')
%}
clear h g lgnd ans;