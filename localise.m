
xcm_all=xcm_all(1:total_molecules);
ycm_all=ycm_all(1:total_molecules);
xf_all=xf_all(1:total_molecules);
yf_all=yf_all(1:total_molecules);
a0_all=a0_all(1:total_molecules);
r0_all=r0_all(1:total_molecules);
off_all=off_all(1:total_molecules);
framenum_all=framenum_all(1:total_molecules);
xf_err_all=xf_err_all(1:total_molecules);
yf_err_all=yf_err_all(1:total_molecules);
a0_err_all=a0_err_all(1:total_molecules);
r0_err_all=r0_err_all(1:total_molecules);
off_err_all=off_err_all(1:total_molecules);
grab_sum_all=grab_sum_all(1:total_molecules);
npix_all=pi*(r0_all.^2) ;    % area of molecule in square pixels
N=npix_all.*a0_all; % number of photons for each molecule
lp2=((r0_all*q).^2+(q^2)/12)*1./N+8*pi*((r0_all*q).^4)*(bkgn^2)/(q^2)*1./(N.*N);
lp=1.3*sqrt(lp2); 

lp(lp<0.005)=[];

 figure
hist(lp*1000,100); %makes a histogram plot
        xlabel('Loc. Prec. (nm)'); %labels the x axis fo the histogram plot
        ylabel('# molecules'); %labels the y axis of the histogram plot
          xlim([0 100])
%         ylim([1 10]) 
        pbaspect([1 1 1])
        hold on
        set(gca,'FontSize',15)
        
      
        