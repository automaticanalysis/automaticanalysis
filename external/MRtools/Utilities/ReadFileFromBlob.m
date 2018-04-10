function [m, h] = ReadFileFromBlob(bin)
% This function is for stripping image data of binary streams
%
% https://brainder.org/2012/09/23/the-nifti-file-format/
% https://brainder.org/2015/04/03/the-nifti-2-file-format/

if iscell(bin)
    bin = bin{1};
end

hb = typecast(bin(1:4),'int32');
if hb == 348
    DIM = typecast(bin(41:56),'int16');
    dim = DIM(2:DIM(1)+1)';
    
    bitdepth = typecast(bin(71:72),'int16');
    
    ras = [typecast(bin(281:296), 'single') ...
           typecast(bin(297:312), 'single') ...
           typecast(bin(313:328), 'single') ...
           [0 0 0 1]']';

    tmp = [1 1 1 1]*ras;
    ras(1:3,4) = ras(1:3,4)-tmp(1:3)';
       
    h = struct(...
        'fname',   'aFile.nii',...
        'dim',     double(dim(1:3)),...
        'dt',      double([bitdepth 0]),...
        'mat',     double(ras),...
        'pinfo',   [double(typecast(bin(113:116),'single')) double(typecast(bin(117:120),'single')) 352]',...
        'descrip', char(bin(149:228)'));
    if h.pinfo(1)==0;
        h.pinfo(1) = 1;
    end
    
    h = spm_create_hdr(h);

    switch bitdepth 
        case 2
            m = typecast(bin(353:end),'uint8');
            m = reshape(m,dim);
        case 4
            m = typecast(bin(353:end),'int16');
            m = reshape(m,dim);    
        case 8
            m = typecast(bin(353:end),'int32');
            m = reshape(m,dim);
        case 16
            m = typecast(bin(353:end),'single');
            m = reshape(m,dim);
        case 64
            m = typecast(bin(353:end),'double');
            m = reshape(m,dim);
        case 256
            m = typecast(bin(353:end),'int8');
            m = reshape(m,dim);    
        case 512
            m = typecast(bin(353:end),'uint16');
            m = reshape(m,dim);   
        case 768
            m = typecast(bin(353:end),'uint32');
            m = reshape(m,dim);   
        otherwise
            disp('this bitdepth has not been configured');
            m = NaN;
    end
end


if hb == 540
    error('this is a nifti-2 format file, and reading this file has not been configured'); 
end


%%% Full NIFTI-1 header info
% hdr = bin(1:348)
% 
% h = [];
% h.hb = typecast(hdr(1:4),'int32');
% h.data_type = char(hdr(5:14)');
% h.db_name = char(hdr(15:32)');
% h.extents = typecast(hdr(33:36),'int32');
% h.session_error = typecast(hdr(37:38),'int16');
% h.regular = char(hdr(39));
% h.dim_info = char(hdr(40));
% h.dim = typecast(hdr(41:56),'int16');
% h.intent_p1 = typecast(hdr(57:60),'single');
% h.intent_p2 = typecast(hdr(61:64),'single');
% h.intent_p3 = typecast(hdr(65:68),'single');
% h.intent_code = typecast(hdr(69:70),'int16');
% h.datatype = typecast(hdr(71:72),'int16');
% h.bitpix = typecast(hdr(73:74),'int16');
% h.slicestart = typecast(hdr(75:76),'int16');
% h.pixdim = typecast(hdr(77:108),'single');
% h.vox_offset = typecast(hdr(109:112),'single');
% h.scl_slope =  typecast(hdr(113:116),'single');
% h.scl_inter =  typecast(hdr(117:120),'single');
% h.sliceend = typecast(hdr(121:122),'int16');
% h.slicecode = char(hdr(123:123));
% h.xyzt_units = char(hdr(124:124));
% h.cal_max = typecast(hdr(125:128),'single');
% h.cal_min = typecast(hdr(129:132),'single');
% h.slice_duration = typecast(hdr(133:136),'single');
% h.toffset = typecast(hdr(137:140),'single');
% h.glmax = typecast(hdr(141:144),'int32');
% h.glmin = typecast(hdr(145:148),'int32');
% h.descrip = char(hdr(149:228)');
% h.aux_file = char(hdr(229:252)');
% h.qform_code = typecast(hdr(253:254),'int16');
% h.sform_code = typecast(hdr(255:256),'int16');
% h.quatern_b = typecast(hdr(257:260),'single');
% h.quatern_c = typecast(hdr(261:264),'single');
% h.quatern_d = typecast(hdr(265:268),'single');
% h.quatern_x = typecast(hdr(269:272),'single');
% h.quatern_y = typecast(hdr(273:276),'single');
% h.quatern_z = typecast(hdr(277:280),'single');
% h.srow_x = typecast(hdr(281:296),'single');
% h.srow_y = typecast(hdr(297:312),'single');
% h.srow_z = typecast(hdr(313:328),'single');
% h.intent_name = char(hdr(329:344)');
% h.magic = char(hdr(345:348)');
% 
% h2 = h;
% [h2.srow_x h2.srow_y h2.srow_z [0 0 0 1]']'
