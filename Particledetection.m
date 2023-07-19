%%%%%%% PSF dectection%%%%%%%%%%%%%%
%%%%%%%%rewritten by Jigmi%%%%%%%%%%%%
clc %Clears the command window
warning off
global xpix ypix wbox; %sets global variables that can be referenced in other m files

imagefile(:,:,15) = strcat('E:\3Dvolummedata\D2HA\cell4\','15','.tif');
% imagefile(:,:,2) = strcat('G:\3Dvolummedata\D2HA\cell4\','3','.tif');
% imagefile(:,:,3) = strcat('G:\3Dvolummedata\D2HA\cell4\','4','.tif');
% imagefile(:,:,4) = strcat('G:\3Dvolummedata\D2HA\cell4\','5','.tif');
% imagefile(:,:,5) = strcat('G:\3Dvolummedata\D2HA\cell4\','6','.tif');
% imagefile(:,:,6) = strcat('G:\3Dvolummedata\D2HA\cell4\','7','.tif');
% imagefile(:,:,7) = strcat('G:\3Dvolummedata\D2HA\cell4\','8','.tif');
% imagefile(:,:,8) = strcat('G:\3Dvolummedata\D2HA\cell4\','9','.tif');
% imagefile(:,:,5) = strcat('G:\3Dvolummedata\D2HA\cell4\','6','.tif');

data_dir = 'E:\3Dvolummedata\D2HA\cell4\'; %this is apparently the path to find the video file (tif?) that contains the data
base_name='10'; %Defining variable, the use at this point in the code is unknwon
imagefile1 = strcat(data_dir,base_name,'.tif'); %This is the total file path to the tif file including the filename and extention.



rbox=4; %Box radiu
q=0.60; %pixel size in um
wvlnth=571/1000; %convert wavelength from nm to um -- where um is micrometers  
NA=1.3; %NA is the numerical aperture- this stands for the resolving power of the microscope
psf_scale=1.2; %The scale of the point spread function
pix_to_pho =17; %40; %Thirty one point two pixels per photon? This is to let us turn pixel counts into actual photons. Perhaps a photon hits more than one pixel each time?
min_thresh  = 23; %minimum threshold for a "bright object"  

%%%%% used for  corelation approach%%%%%%%%
B1 = imread('temp.tif'); %%%%templat
B1=B1/pix_to_pho;
cutoff=0.32;%%%nominal 0.35
pad=5;



box_overlap_factor = 2; %if center to center closer, don't include either -- if the fitting boxes for two bright pixels overlap, eliminate both in the data because we do not have a model suitable for overlapping fits.
w_mask = round(rbox*box_overlap_factor); %width of masking when a "high" pixel is identified -- The mask is used to check to see if two pixels are close enough to warrant throwing away. 
wbox=2*rbox+1; %Box width is going to be one unit larger than rbox, so 7 wide. The +1 might be to make it so there is a center pixel.
[xpix0,ypix0] = meshgrid(-2*rbox:2*rbox,-2*rbox:2*rbox); %This creates a domain in x and y pixels. The specific domain is from -6 to 6 units in x and y.
[xpix,ypix] = meshgrid(-rbox:rbox,-rbox:rbox);%This creates another region from -3 to 3 in x and y. Or from -rbox to rbox.
psf_w0 = psf_scale*0.55*wvlnth/NA/1.17; %PSF is the point spread function. I believe in this case it will be an airy function. 
psf_std=psf_w0/2; %standard deviation of point spread function
psf_w02=(psf_w0/q)*(psf_w0/q); %square of 1/e^2 radius in pixels-- so this is the point spread function divided by pixel size squared. 
rball=6; %radius of rolling ball -- The rolling ball is used to model the molecules that we are imaging next line defines the ball
se = strel('ball',rball,rball,0); %structural element, i.e. rolling ball  -- creates a matrix representing 3d geometry-- a ball in this case. this is a strel data type
FWHM=1; %FWHM of gaussian smoothing in pixels -- This is fairly straightforward, the full width at half maximum of some gaussian. 
rk=(FWHM)/sqrt(2*log(2)); %1/e^2 smoothing radius in pixels -- not sure what the smoothing radius is for right now.
kw=20; %kernal width of smoothing function -- Don't know what a kernel width is.
[Xgs,Ygs]=meshgrid(-kw/2:kw/2,-kw/2:kw/2); %creating a domain in xgs and ygs that goes from -10 to 10 in x and y. (kw/2)
kd=sqrt(Xgs.*Xgs+Ygs.*Ygs); %Multiplies each element of the xgs matrix by the xgs matrix and each element of the ygs matrix times the ygs matrix. this is e^sum of those
gs=exp(-2*kd.*kd/(rk*rk)); %gs is a small gaussian of FWHM of...
gs=gs/sum(sum(gs)); %smoothing function normalized to have area = 1 -- so GS is being divided by the 'integral' of gs.

%initialize arrays
n_init=500000; %Five hunderd thousand is a lot. Seems like maybe a counting index? Maybe?
xcm_all=zeros(n_init,1); %X array of all zeroes from index 1 (or 0 maybe) to 500000. One row.
ycm_all=zeros(n_init,1); %Set to mimic the x array
xf_all=zeros(n_init,1); %Same stuff. Xf stands for??
yf_all=zeros(n_init,1); %Yf stands for??
a0_all=zeros(n_init,1); %a0 might be some kind of constant. Not sure
r0_all=zeros(n_init,1); %r0 might be radius related
off_all=zeros(n_init,1); %off_all seems like a context sensitive thing. What does it mean?
framenum_all=zeros(n_init,1); %Which frame of the data are we looking at? These are big arrays!
xf_err_all=zeros(n_init,1); %Xf error maybe. Zeroes.
yf_err_all=zeros(n_init,1); %Same for y
a0_err_all=zeros(n_init,1); %error for a0
r0_err_all=zeros(n_init,1); %error for r0
off_err_all=zeros(n_init,1); %error for off all?
grab_sum_all=zeros(n_init,1); %total sum matrix? Needs contet.

total_molecules=0; %there are no molecules!
n_fail_a0=0; % maybe this is the number of failed curve fits.
n_fail_outbox=0; %Related to failure but not sure how.



%auto calc bkgnd pixelation noise noise if bkgn is commented out above
if(exist('bkgn','var')==0) %If you did not comment out the bkgn statement and thus bkgn does not exist
    bkgn=bg_noise_calc01(imagefile1,pix_to_pho); %find the background noise automatically using an outside function.
        
%     answer = questdlg(sprintf(['Noise: ',num2str(bkgn,'%2.2f'),'. Yes to continue, No to restart.']));
    answer = questdlg(sprintf(['Noise: ',num2str(bkgn,'%2.2f'),', Thresh: ',num2str(min_thresh),'. Yes to continue, No to restart.']));

        %Displays a question dialogue and stores the result in the variable
        %answer. Sprintf displays the value of a variable in string
        %format. Asks whether or not the calculated background is okay.
    if ~strcmp(answer,'Yes') %if it IS okay
        return; %Return the calculated background to matlab
    end %end apparently must be called after 'if'
end %end apparently must be called after 'if'



tic %begins time for performance analysis.
for plane=1%1:10  %% z plane #
for fno=15%:1:5%data stack address
current_imagefile=imagefile(:,:,fno);
fn01=fno;
k=num2str(fn01);
kk=strcat(k,'.tif')
info=imfinfo(kk);
jkl=length(info)
for fileloop=plane:1:jkl %Select increment to process same plane in next subsequent cycles 
    drawnow; %Update the waitbar on the screen.
    fprintf ('loading image ... %d\n', fileloop)
    i1=double(imread(current_imagefile,fileloop))/pix_to_pho; % covert into photon count 
    [yw,xw]=size(i1); %How big is it? store that in this matrix thing!
    y_cord=1:yw;
    x_cord=1:xw;
   
    %%%%% bacground subtraction rolling ball algo %%%%%%%
    i1_gs = uint16(conv2(i1,gs,'same')); %smoothed 
    bkg = double(imopen(i1_gs,se)); %morphological oepning 
    iprod=i1-bkg; %subtracts the background from the original image 
    iprod=iprod.*(iprod>0); %set negative values to 0 
    high_pixel_mask = zeros(yw,xw); %initializes the high pixel mask matrix to zeroes
    n_boxes=0; % the number of boxes?
    boxes_xy=zeros(10000,3); %make a matrix for the x and y positions of these boxes
%%%%%%%%%%%%%%%%%%%%% threshold  cutoff approach%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%
% [iy,ix] = find(iprod >= min_thresh); %find all pixels in iprod above the threshold (returns the indecies of such pixels.)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%correlation approach cut off approach %%%%%%%%% 
    normx_corrmap=normxcorr2(B1(:,:,1),iprod);
    [x1, y1]= size(normx_corrmap);
    [xt1, yt1]=size(B1);
    gg=normx_corrmap(floor(xt1/2)+1:x1-floor(xt1/2),floor(yt1/2)+1:y1-floor(yt1/2));
    gg1=gg(pad:size(gg,1)-pad-1,pad:size(gg,2)-pad-1);  
    my_image= padarray(gg1,[pad pad],0);
    %maxptx = max(normx_corrmap(:));
    %surf(my_image),shading flat
    %my_image=gg2;
    image_thresholded = my_image;
    %image_thresholded(my_image>=.4) = my_;
    image_thresholded(my_image<cutoff) = 0;

    [iy,ix] = find(image_thresholded >= cutoff);    
%     figure,surf(gg1),shading flat
%     figure,surf(i1),shading flat
%     figure,surf(image_thresholded),shading flat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% select only  useful high thresholds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    high_pix_inframe = size(iy,1); 
    high_pix_index = 1;
    while high_pix_index <= high_pix_inframe %So long as truth is truth, so long as the sky be blue...
        pix_val = 0; %intitalize/restart this variable!
        
        if (iprod(iy(high_pix_index),ix(high_pix_index))>pix_val && high_pixel_mask(iy(high_pix_index),ix(high_pix_index))==0  && iy(high_pix_index)<yw-rbox-1 && iy(high_pix_index)>rbox+1  && ix(high_pix_index)<xw-rbox-1 && ix(high_pix_index)>rbox+1)
            pix_val = iprod(iy(high_pix_index),ix(high_pix_index));
            high_pixel_y = iy(high_pix_index);  
            high_pixel_x = ix(high_pix_index);
            high_pix_index = high_pix_index + 1;
        else
            high_pix_index = high_pix_index + 1;
            continue
        end

        x0_box=high_pixel_x-rbox;%Make a
        y0_box=high_pixel_y-rbox;%new region
        x1_box=high_pixel_x+rbox;%centered around
        y1_box=high_pixel_y+rbox;%the bright pixel (from x0 to x1 and from y0 to y1

        x0_mask=high_pixel_x-w_mask; %this is the left-bound edge of our bounding box
        if x0_mask < 1 %make sure the mask is
            x0_mask = 1; %at LeAST one pixel wide.
        end %endif!
        x1_mask=high_pixel_x+w_mask; %x0 and x1 are opposite sides of the box.
        if x1_mask > xw
            x1_mask = xw;
        end
        y0_mask=high_pixel_y-w_mask; %same things done for y0 and y1
        if y0_mask < 1
            y0_mask = 1;
        end
        y1_mask=high_pixel_y+w_mask;
        if y1_mask > yw
            y1_mask = yw;
        end%endif endif endif...

        high_pixel_mask(y0_mask:y1_mask,x0_mask:x1_mask)=1; %so here we make a matrix that is all ones that coveres the region of interest for a particular bright pixel.

        grab=iprod(y0_box:y1_box,x0_box:x1_box); %maybe this extracts a portion of the image around th region of interest?
        grab_sum=sum(sum(grab)); %this is perhaps the total value of brightness within this region.

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%calculate C.O.M. for intial guess fit %%%%%%%%%%%%%%%%%
        
        xm_sum=0;
        ym_sum=0;
        m_sum=0;

        for i=x0_box:x1_box %so long as i (the index) is between x0 and x1 box. So, within the x region
            for j=y0_box:y1_box %And so long as our index j is within the y region:
                xind=floor(i);%Start with x index at the lowest i
                yind=floor(j);%start with y index at the lowest j
                intens=iprod(yind,xind); % not sure what this means exactly.
                xm_sum=xm_sum+xind*intens; %these depend on the intensity (intens maybe) and is some sum of values. Perhaps the total brightness in the image?
                ym_sum=ym_sum+yind*intens;
                m_sum=m_sum+intens;
            end%endif
        end%endif

        x_cm=xm_sum/m_sum; 
        y_cm=ym_sum/m_sum;

        xc_box=(x0_box+x1_box)*0.5;
        yc_box=(y0_box+y1_box)*0.5;

        xguess=x_cm-xc_box;
        yguess=y_cm-yc_box;

        for i=1:wbox
            for j=1:wbox
                k=(i-1)*wbox+j;
                xymerge(k)=0;
                zmerge(k)=grab(i,j);
            end
        end 

        beta0=[xguess,yguess,50,psf_w0/q,min(grab(:))]; % x,y,a0,r0,offset -- this is perhaps for the gaussian, 
        [betafit,resid,J,COVB,mse] = nlinfit(xymerge,zmerge,@gaussian_merge3,beta0);%nlinfit is a non linear regression fit 
        ci = nlparci(betafit,resid,'covar',COVB); %calculate error estimates on parameters -- how does it do this again? Comparing data to the fit?
        ci_err=(ci(:,2)-ci(:,1))/2; %Confusing syntax

        yf=betafit(2)+yc_box;
        xf=betafit(1)+xc_box;
        a0=betafit(3);
        r0=abs(betafit(4));
        off=betafit(5);
        
        failed=0; %flag for failed localization
        if(a0 < 0)
            n_fail_a0=n_fail_a0+1;%this goes back to the fail matrix  and says whether or not the fit was successful. Has it been localized?
            failed=1;
        end
        if(xf > x1_box || xf < x0_box || yf > y1_box || yf < y0_box) %this is also checking for failure but where it checked a0 before it checks to see if the
            %gaussian goes outside the region of interest, maybe?
            n_fail_outbox=n_fail_outbox+1;
            failed=1;
        end
        
        

        if(failed==0) %assign if fit criteria is satisfied -- so if it did NOT fail, go ahead and write the fit to a series of variables
            total_molecules=total_molecules+1;
            xcm_all(total_molecules)=x_cm;
            ycm_all(total_molecules)=y_cm;
            xf_all(total_molecules)=xf;
            yf_all(total_molecules)=yf;
            a0_all(total_molecules)=a0;
            r0_all(total_molecules)=r0;
            off_all(total_molecules)=off;
            
            
            framenum_all(total_molecules)=fileloop;
            xf_err_all(total_molecules)=ci_err(1);
            yf_err_all(total_molecules)=ci_err(2);
            a0_err_all(total_molecules)=ci_err(3);
            r0_err_all(total_molecules)=ci_err(4);
            off_err_all(total_molecules)=ci_err(5);
            grab_sum_all(total_molecules)=grab_sum;
                       
            n_boxes=n_boxes+1; %increase the box index
            boxes_xy(n_boxes,1)=high_pixel_x;
            boxes_xy(n_boxes,2)=high_pixel_y;
            boxes_xy(n_boxes,3)=1;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% display detections %%%%%%%%%%%%%%%%%%
    imshow(i1,[]) %get the data ready o draw as an image
    title(['Frame: ' num2str(fileloop) ', ' num2str(total_molecules) ' molecules']); %put some information on the image as to where it came from and the data within
    hold on %freeze that output
    draw_boxes(n_boxes,boxes_xy,rbox); %draw boxes is an outside function. It might show the fitted gaussian regions superimposed on the image.
    hold off %unfreeze the output
    drawnow; %draw what you've got ready to draw!

end %endforloop
end

%%%%%%%%%%%%%%%%%%%%%%%%%save data points%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fnm=sprintf('Workspace_%d.mat',plane);
% save(fnm);
%  clear xcm_all ycm_all xf_all yf_all a0_all r0_all
%  total_molecules=0;
end
toc


