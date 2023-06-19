figure
tot=1000 % put number of molecules used to render

xf_all=xf_all(1:tot);
yf_all=yf_all(1:tot);
lp=lp(1:tot);
N=N(1:tot);

nfiles=1;%No. of data files to be combined
n_species=1; %No. of data files to be combined
Int_weights=[0.02 .005 0.2];  %|Weight Changed from 0.0002 to 0.02|
S1_RGB=[1 1 1];% RGB color units for species rendering

exf=4; %expansion factor
size_fac=2; %factor to scale size of points


ishift=0;%20; %xshift
jshift=0;%80; %yshift
left1=1;
right1=401;
bottom1=1;
top1=221;
xstart=left1;
xend=right1;
ystart=bottom1;
yend=top1;


xw=xend-xstart;%width of image        
yw=yend-ystart;%length of image       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%scale bar options
include_bar='y';
x0_bar=200; %coordinates of scale bar
y0_bar=100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Other Plots 
plot_hist_lp='n';
plot_hist_ratio='y';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xf_all_cum=zeros(0);
yf_all_cum=zeros(0);
nrat_cum=zeros(0);


for file_index=1:nfiles

if file_index>1
    xf_all_cum=[xf_all_cum; xf_all];
    yf_all_cum=[yf_all_cum; yf_all];
    lp_cum=[lp_cum; lp];
    if n_species>1
       nrat_cum=[nrat_cum; nrat];
    end
    N_cum=[N_cum; N];
else
    xf_all_cum=xf_all;
    yf_all_cum=yf_all;
    lp_cum=lp;
    if n_species>1
       nrat_cum=nrat;
    end
    N_cum=N;
end
%clear xf_all yf_all lp nrat N
end

xf_all=xf_all_cum;
yf_all=yf_all_cum;
lp=lp_cum;
if n_species>1
   nrat=nrat_cum;
end
N=N_cum;

xf=exf*(xf_all-ishift);
yf=exf*(yf_all-jshift);
lppix=lp/q*exf;   % width of each molecule in pixels, based on localization-uncertainty 
lp2pix=lppix.*lppix;     

%|......................................................................|
%|Determining the position of molecule in pixel numbers|
xc=uint16(xf);    %|x-coordinate of the molecule|
xc1=xc;
yc=uint16(yf);    %|y-coordinate of the molecule|
yc1=yc;
%lppix=lppix ;     %|localization precision of the molecule|
lppix1=lppix./(max(lppix)-min(lppix)) ;     %|Normalization (kind-of) for generating COLORMAP|
lppix1=lppix1./max(lppix1);

if n_species==1
    [image_current,n_rendered]=func_1color(xf_all,yf_all,N,lppix,exf,size_fac, Int_weights(1), ishift, jshift, xw, yw);
    impts=zeros(yw*exf,xw*exf,3);
    impts(:,:,1)=S1_RGB(1)*image_current;
    impts(:,:,2)=S1_RGB(2)*image_current;
    impts(:,:,3)=S1_RGB(3)*image_current;   
end

impts=1.2*impts;
thresh=impts*0+1;
impts=thresh.*(impts>thresh)+impts.*(impts<=thresh);

if(include_bar=='y')
        um1=1/q*exf;
        um1_round=round(um1);
        bar_width=12;
        for i=1:um1_round
            for j=1:bar_width
                xi=i+x0_bar;
                yi=j+y0_bar;
                impts(yi,xi,:)=1;
            end
        end

    nm250=0.25/q*exf;
    nm250_round=round(nm250);
    bar_width=4;
    y0_bar=y0_bar+bar_width;
    for i=1:um1_round
        for j=1:bar_width
            ival=round(i/nm250_round);
            ival=mod(ival,2);
            xi=i+x0_bar;
            yi=j+y0_bar;
            impts(yi,xi,:)=double(ival);
        end
    end
end

%%%%%%%%%%%%%%%
%%%%%%%%
% subplot(1,4,4)
AA=flip(impts,1);
Q=AA(:,:,1);
Q=AA(:,:,1);
imshow(Q,[]), axis image,colormap hot ,box on
%title(titletext);
% axis image, colormap gray