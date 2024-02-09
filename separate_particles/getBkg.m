function bkg = getBkg(subfolders,fname_prefix,Istart,Iend,increment,demosaicpattern)


fname = [subfolders filesep fname_prefix filesep 'B' num2str(Istart,'%05d') '.tif' ];

Im = imread(fname);%,'ppm','uint16');
cc = class(Im);
if ~isempty(demosaicpattern)
    Im = rgb2gray(demosaic(Im,demosaicpattern));
    cc = class(Im);
end

Im = double(0*Im);

Nim = 0;
for k = Istart : increment : Iend
    Nim=Nim+1;
    fname = [subfolders filesep fname_prefix filesep 'B' num2str(k,'%05d') '.tif' ];
    if ~isempty(demosaicpattern)
        Im = Im+double(rgb2gray(demosaic(imread(fname),demosaicpattern)));
    else
        Im = Im+double((imread(fname)));
    end
end
% figure;imagesc(Im/Nim)
bkg = cast(Im/Nim,cc);

