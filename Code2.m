%%
clear;
clc;
close all;
image_name = 'colours.jpg';
I_colour = imread(image_name);
% I_colour = uint8(zeros(300,450,3));
% I_colour(1:150,:,1) = 255;
I = (rgb2gray(I_colour));
I = double(I);
iso = I/255;
I_colour  = double(I_colour);
figure(1), imshow(I_colour/255);

hand = images.roi.AssistedFreehand;
draw(hand);
binaryImage = hand.createMask();
binaryImage = double(binaryImage);

source = binaryImage;

figure(2),imshow(binaryImage);
figure(3),imshow(iso);
%%
[m,n] = size(I);
Inew = zeros(m,n);
Inew = double(Inew);
Inew_colour = double(zeros(m,n,3));
for i=1:m
    for j=1:n
        if binaryImage(i,j) == 0
            Inew(i,j) = I(i,j);
            Inew_colour(i,j,:) = I_colour(i,j,:);
        else
            Inew(i,j) = 255;
            Inew_colour(i,j,:) = 255;
        end
    end
end
figure, imshow(Inew_colour/255);
C = double(~binaryImage);
%%
while(check(binaryImage))
    [xc,yc] = get_contour(binaryImage);
    [slopex,slopey] = get_normal(xc,yc);
    [Cont_C,Cont_D,gmag] = gettarget(Inew,xc,yc,binaryImage,C,slopex,slopey);
    [s,t,cmax] = get_max(Cont_C,Cont_D,xc,yc);
    [p,q] = sample(Inew_colour,source,binaryImage,s,t,gmag);
    [Inew,Inew_colour,binaryImage,C] = update(Inew,Inew_colour,binaryImage,C,p,q,s,t,cmax);
    figure(5),imshow(Inew_colour/255);
    hold on
    plot(t,s,'r+');
    plot(q,p,'g+');
    hold off
end
%%
function [xc,yc] = get_contour(binaryImage)
xc = [];
yc = [];
[m,n] = size(binaryImage);
for i=1:m
    for j = 1:n
        if binaryImage(i,j) == 1
            if i==1 || i== m || j == 1 || j == n
                xc = [xc,i];
                yc = [yc,j];
            else
                if binaryImage(i-1,j)==0 ||  binaryImage(i,j-1)==0 || binaryImage(i,j+1)==0 || binaryImage(i+1,j)==0
                    xc = [xc,i];
                    yc = [yc,j];
                else
                    if binaryImage(i-1,j-1)==0 ||  binaryImage(i+1,j-1)==0 || binaryImage(i-1,j+1)==0 || binaryImage(i+1,j+1)==0
                        xc = [xc,i];
                        yc = [yc,j];
                    end
                end
            end
        end
    end
end
end
%%
function [slopex,slopey] = get_normal(xc,yc)
if size(xc,2)==1
    slopex = 1/sqrt(2);
    slopey = 1/sqrt(2);
else
xCenter = mean(xc);
yCenter = mean(yc);
angles = atan2d((yc-yCenter),(xc-xCenter));
[~, sortIndexes] = sort(angles);
xc = xc(sortIndexes);  % Reorder x and y with the new sort order.
yc = yc(sortIndexes);
t = cumsum([0;diff(xc(:)).^2 + diff(yc(:)).^2]);
splx = spline(t,xc);    
sply = spline(t,yc);
slopex = ppval(fnder(splx),t);
slopey = ppval(fnder(sply),t);
slopex = slopex./((slopex.^2+slopey.^2).^(0.5));
slopey = slopey./((slopex.^2+slopey.^2).^(0.5));
tanslope = slopey./slopex;
%figure(100),quiver(xc',yc',slopex,slopey)
end
end
%%
function [Cont_C,Cont_D,gmag] = gettarget(Inew,xc,yc,binaryImage,C,slopex,slopey)
len = size(xc,2);
Cont_C = double(zeros(1,len));
Cont_D = double(zeros(1,len));
[m,n] = size(binaryImage);
for i = 1:len
    p = xc(1,i);
    q = yc(1,i);
    sum = 0;
    for j = max(1,p-4):min(p+4,m)
        for k = max(1,q-4):min(q+4,n)
            if binaryImage(j,k)==1
                continue
            end
            sum = sum + C(j,k);
        end
    end
    Cont_C(1,i) = sum/81;
end

[Gx, Gy] = imgradientxy(Inew);
[gmag,gdir] = imgradient(Gx,Gy);
se = strel('ball',3,3);
J = imdilate(255*uint8(binaryImage),se);
Gx = Gx.*(~double(J/255));
Gy = Gy.*(~double(J/255));
gmag = gmag.*(~double(J/255));
for i = 1:len
    p = xc(1,i);
    q = yc(1,i);
    a = slopex(i); 
    b = slopey(i);
    G1 = [b, -a]/sqrt(a^2+b^2);
    high = -1*Inf;
    for j = max(1,p-4):min(p+4,m)
        for k = max(1,q-4):min(q+4,n)
            c = Gx(j,k);
            d = Gy(j,k);
            G2 = [d, -c];
            cont_d = abs(dot(G1,G2));
            if high<cont_d
                high = cont_d;
            end
        end
    end
    Cont_D(1,i) = high;
end
end
%%
function res = distance(Inew_colour,source,bImage,x,y,s,t)
    [m,n,~] = size(Inew_colour);
%     [m,n] = size(Inew_colour);
    res = 0;
    flag = 0;
    for a = -4:4
        for b = -4:4
            if s+a<1 || s+a>m || t+b<1 || t+b>n
                continue
            end
            if source(x+a,y+b)== 0
                if bImage(s+a,t+b) == 0
                    flag = 1;
%                     res = res + (double(Inew_colour(x+a,y+b))-double(Inew_colour(s+a,t+b)))^2;
                    res = res + (double(Inew_colour(x+a,y+b,1))-double(Inew_colour(s+a,t+b,1)))^2;
                    res = res + (double(Inew_colour(x+a,y+b,2))-double(Inew_colour(s+a,t+b,2)))^2;
                    res = res + (double(Inew_colour(x+a,y+b,3))-double(Inew_colour(s+a,t+b,3)))^2;
                end
            else
                res = res + Inf;
            end
        end
    end
    if flag == 0 
        res = Inf;
    end
end

function [p,q] = sample(Inew_colour,source,bImage,s,t,gmag)
   	[m, n,~] = size(Inew_colour);
%     [m,n] = size(Inew_colour);
    min_dist = Inf;
    grad = Inf;
    for i = 1:m
        for j = 1:n
            if i-4<1 || j-4 <1 || i+4>m || j+4>n
                continue
            end
            d = distance(Inew_colour,source,bImage,i,j,s,t);
%             max_grad = max(max(gmag(max(i-4,1):min(i+4,m),max(j-4,1):min(j+4,n))));
            max_grad = gmag(i,j);
            if d<min_dist || ((d==min_dist) && (grad>max_grad)) 
                grad = max_grad;
                p = i;
                q = j;
                min_dist = d;
            end
        end
    end  

end
%%
function [Inew,Inew_colour,bImage,C] = update(Inew,Inew_colour,bImage,C,p,q,s,t,cmax)
    [m,n] = size(Inew);
    for i = -4:4
        for j = -4:4
            if s+i<1 || s+i>m || t+j<1 || t+j>n
                continue
            end
            if bImage(s+i,t+j)==0
                continue
            end
            Inew(s+i,t+j) = Inew(p+i,q+j);
            Inew_colour(s+i,t+j,:) = Inew_colour(p+i,q+j,:);
            C(s+i,t+j) = cmax;
            bImage(s+i,t+j) = 0;
        end
    end
end
%%
function [condition] = check(binaryImage)
    condition = 0;
    [m,n] = size(binaryImage);
    for i = 1:m
        for j=1:n
            if binaryImage(i,j)==1
                condition = 1;
            end
        end
    end
end
%%
function [s,t,cmax] = get_max(Cont_C,Cont_D,xc,yc)
	len = size(Cont_C,2);
	P = double(zeros(1,len));
	for i = 1:len
		P(1,i) = Cont_C(1,i)*Cont_D(1,i);
%         P(1,i) = Cont_D(1,i);
    end
	[~,pos] = max(P);
	s = xc(1,pos);
	t = yc(1,pos);
    fprintf('%f %f\n',s,t);
	cmax = Cont_C(pos);
end
%%

