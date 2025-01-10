function Imp=preprocessing(Im,th,part_radius,gain,bkg)

%Im = imcomplement(Im);


%se=strel('disk',strel_size); %opening to remove big elements
%imo=imopen(Im,se);
Ims=imsubtract(double(Im),double(bkg));
%Ims = Im;
%Im = bpass2(Im,lnoise,part_size);
%Im = imgaussfilt(Im,part_radius);
%Im=imnlmfilt(Im);

% Th = th*mean(mean(Im));
Ims(Ims<th)=0; %thresholding to remove the background
%Im=medfilt2(Im,[3,3]);
%Im = imadjust(Im,stretchlim(Im),[]);

%Im=imsharpen(Im,'Radius',6);
%Im(Im<th)=0; %thresholding to remove the background

%se=strel('disk',1); %opening to remove little elements
%Im=imerode(Im,se);

%Imbw=imregionalmax(Im,6);
%Im = Im/mean(Im(Imbw))*gain;

tt = class(Ims);
%Im = cast(double(Im)/max(max(double(Im)))*gain,tt);
Ims = cast(double(Ims)*gain,tt);


Imp= remove_part(Ims,part_radius);
%[Imt,Imp]= remove_part_byintensity(Im,radii_range,intensity_thr);
 %[Imt,Imp]= remove_part(Im,9);

