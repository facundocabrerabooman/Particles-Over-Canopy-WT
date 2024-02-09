function Imp = remove_part(Im,PartRadius)
%%
%st = strel('disk',1);
stdil = strel('disk',round(PartRadius/2));

%Imd = imdilate(imopen(Im,st),stdil);
Imd = imdilate(Im,stdil);
Imp = Imd;

%% refine particle image
% Rp = PartRadius;
% 
% Rmin = min(0.1*Rp,2);
% Rmax = max(4*Rp,10);
% 
% [C, R] = imfindcircles(Imd, cast([Rmin Rmax],class(Im)));
% if size(C,1) == 1
%     roi=images.roi.Circle('Center',C,'Radius',R);
%     mask = createMask(roi,size(Imd,1),size(Imd,2));
% else
%     mask = Imd;
% end
% Imp = immultiply(Im,cast(mask,class(Im)));