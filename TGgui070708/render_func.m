impts=zeros(yw*exf,xw*exf);
n_rendered=0;
weight=str2double(get(handles.weight,'String'));
size_fac=str2double(get(handles.size_fac_edit,'String'));

for i=nstart:nend
    if xc(i)>=1 && yc(i)>=1 && xc(i)<xw*exf && yc(i)<yw*exf && N(i)>0
      wide=ceil(size_fac*lppix(i)*1.5+1);
%       if wide>20 
%           wide=20;
%       end
      if xc(i)-wide>=1 && xc(i)+wide<xw*exf && yc(i)-wide>=1 && yc(i)+wide<yw*exf
        n_rendered=n_rendered+1;
        for j=xc(i)-wide:xc(i)+wide
          for k=yc(i)-wide:yc(i)+wide
            dx=double(j)-xf(i);
            dy=double(k)-yf(i);
            int=pi*lp2pix(i)*size_fac;
            a=exp(-2*(dx*dx+dy*dy)/(size_fac*size_fac*lp2pix(i)))*N(i)*weight/int;
            impts(k,j)=impts(k,j)+a;
          end
        end
      end
    end
    waitbarxmod(i/nend,w); %update
end

thresh=impts*0+1;
impts=thresh.*(impts>thresh)+impts.*(impts<=thresh);
clear thresh;