function save_three_ortho_jpgs(img, jpeg_name)
%
% generate 3 jpeg montages given an SPM volume and a filename
%
% the full and proper way to do this requires reslicing the file according
% to its voxel-to-world transformation -- the idea here is a quick-n-dirty 
% generation of images for QA eyeballing -- a section montage in each of
% the three V(:,:,:) planes. These *probably* correspond to sagittal,
% coronal, and axial, but without reslicing, we can't guarantee it.

suffix = { '01', '02', '03' };

% if img is a filename, try to load the spm volume

if (ischar(img))
	try
		header = spm_vol(img);
		img = spm_read_vols(header);
	catch
		fprintf('%s: Cannot load image file\n', mfilename);
		return;
	end
end

% if img is 4D, use the first volume

if (length(size(img)) == 4)
	 img_to_display = img(:,:,:,1);
else
	 img_to_display = img;
end

% jpeg_name can be a full path; bust it up for later use

[ path,fname,~ ] = fileparts(jpeg_name);

% dimension 1

% testing indicates 6x6 images generally gives a nice montage

k_max = size(img_to_display,3);
delta_k = floor(k_max/36);
nrow = 6; ncol = 6;	

% tr_3Dto2d goes: img, k_start, k_end, delta_k, nrows, ncols

section_montage = tr_3Dto2D(img_to_display,1,k_max,delta_k,nrow,ncol); 
section_montage = section_montage/max(max(section_montage(:))); % scaling
section_montage = add_survey_markers(section_montage);
outname = fullfile(path, sprintf('%s_%s.jpg', fname, suffix{1}));
imwrite(section_montage, outname);

% tr_3Dto2D  works on the first two dimensions of the input
% so we have to rotate the image by one to get the next 

img_to_display = shiftdim(img_to_display,1);

% dimension 2

k_max = size(img_to_display,3);
delta_k = floor(k_max/36);
section_montage = tr_3Dto2D(img_to_display,1,k_max,delta_k,nrow,ncol);
section_montage = section_montage/max(max(section_montage(:)));
section_montage = add_survey_markers(section_montage);
outname = fullfile(path, sprintf('%s_%s.jpg', fname, suffix{2}));
imwrite(section_montage, outname);

img_to_display = shiftdim(img_to_display,1);

% dimension 3

k_max = size(img_to_display,3);
delta_k = floor(k_max/36);
section_montage = tr_3Dto2D(img_to_display,1,k_max,delta_k,nrow,ncol); 
section_montage = section_montage/max(max(section_montage(:)));
section_montage = add_survey_markers(section_montage);
outname = fullfile(path, sprintf('%s_%s.jpg', fname, suffix{3}));
imwrite(section_montage, outname);
	

end		

function img = add_survey_markers(img)
	dx = floor(size(img,1)/10);
	dy = floor(size(img,2)/10);
	xmax = size(img,1)-1;
	ymax = size(img,2)-1;
	for x = dx:dx:xmax
		for y = dy:dy:ymax
			img(x:x+1,y:y+1) = 1.0;
		end
	end
end
		