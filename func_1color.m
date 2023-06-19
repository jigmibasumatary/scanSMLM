function [image,n_rendered]=func_1color(xf_current_frame,yf_current_frame,N_current_frame,lppix_current_frame,exf, size_fac, weight, xshift, yshift, xw, yw)
   
    xf=exf*(xf_current_frame-xshift);
    yf=exf*(yf_current_frame-yshift);
    xc=uint16(xf);
    yc=uint16(yf);
    lp2pix_current_frame=lppix_current_frame.^2;
        
    image=zeros(yw*exf,xw*exf);
    n_rendered=0;
    nend=length(xf_current_frame);
    for i=1:nend
        if xc(i)>=1 && yc(i)>=1 && xc(i)<xw*exf && yc(i)<yw*exf && N_current_frame(i)>0
          wide=ceil(size_fac*lppix_current_frame(i)*1.5+1);
%           if wide> 20 
%               wide=20;
% %                  wide=10
%           end
          if xc(i)-wide>=1 && xc(i)+wide<xw*exf && yc(i)-wide>=1 && yc(i)+wide<yw*exf
            n_rendered=n_rendered+1;
            for j=xc(i)-wide:xc(i)+wide
              for k=yc(i)-wide:yc(i)+wide
                dx=double(j)-xf(i);
                dy=double(k)-yf(i);
                int=pi*lp2pix_current_frame(i)*size_fac;

%                 a=exp(-2*(dx*dx+dy*dy)/(size_fac*size_fac*lp2pix_current_frame(i)));
                  a=exp(-2*(dx*dx+dy*dy)/8)*200*weight;  
        
                image(k,j)=image(k,j)+a;
              end
            end
          end
        end    
    end
    image=image/5;
    image=1.*(image>1)+image.*(image<=1);