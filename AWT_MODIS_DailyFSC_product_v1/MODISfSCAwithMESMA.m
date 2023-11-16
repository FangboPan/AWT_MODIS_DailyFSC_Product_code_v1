%% ********************************************************************* %%
% AWT MODIS FSC Cloud Gap Filled product code
% 0-100:FSC, 237:Water£¬250:Cloud£¬253/255:Nodata
% Code by Fangbo Pan and Gongxue Wang, November 15, 2023
% jiang@bnu.edu.cn,panfb@mail.bnu.edu.cn
%% ********************************************************************* %%
function [ ] = MODISfSCAwithMESMA(fileMOD09GA,fileMODfSCA)

if ischar(fileMOD09GA),fileMOD09GA = cellstr(fileMOD09GA);end
if ischar(fileMODfSCA),fileMODfSCA = cellstr(fileMODfSCA);end

tf = checkSize(fileMOD09GA,fileMODfSCA) ||...
    ( checkSize(fileMOD09GA,fileMODfSCA));
if nargin>1
    assert(tf,'input and output filenames must have same size')
end

for i=1:length(fileMOD09GA)
    disp(['estimating fSCA from ',fileMOD09GA{i}])
    getSCA(fileMOD09GA{i},fileMODfSCA{i});
end


function [] = getSCA(fileRefl,fileSCA)

dtype = 'uint8';
cTable = defaultSnowColormap;
TiffTags = struct('Compression',Tiff.Compression.LZW);

[pathStr,fileName,ext] = fileparts(fileSCA);
fielSCA2 = [pathStr,'\',fileName,'.cloudMasked.tif'];

state = GetMOD09GA(fileRefl, 'state');
reflCube = GetMOD09GA(fileRefl, 'allbands');
solarZen = GetMOD09GA(fileRefl,'SolarZenith');

s = size(reflCube);

temp =  fillmissing(reflCube(:,:,6),'nearest');
reflCube(:,:,6)= temp;
clear temp

cloudMask = state.cloud ==1 | state.cloud ==2 | state.cirrus ==3;
lakeMask = state.landwater ==3 | state.landwater ==4 | state.landwater ==5;
oceanMask = state.landwater ==6 | state.landwater ==7;

cloudMask = imresize(cloudMask,[s(1),s(2)],'nearest');
lakeMask = imresize(lakeMask,[s(1),s(2)],'nearest');
oceanMask = imresize(oceanMask,[s(1),s(2)],'nearest');
%%
% -------------------------------------------------------------- %
% extract endmembers from the current scene

% [Ctrol-Endmember-threshold,SnowNDSI-threshold,SnowNDVI-threshold,VegNDSI-threshold,VegNDVI-threshold,
%  SoilNDSI-threshold,Soil-minNDVI-threshold,Soil-maxNDVI-threshold,WaterNDWI-threshold]
ctrThresh = [0.75,0.75,-0.035,-0.4,0.7,-0.4,0.0,0.15,0.2];
%  [Snow-0.55um,Snow-0.55um_threshold,Water-0.86um,Water-0.86um_threshold]
bandThresh = [1,0.7,2,0.2];
%%
[fSCA,fVCA,unmixError] = retrieveSCAforSensor(reflCube,'MODIS',ctrThresh,bandThresh);

fSCA(lakeMask) = 237;
fSCA(oceanMask) = 239;

fSCA = cast(fSCA,dtype);


maskedSCA = fSCA;
maskedSCA(cloudMask) = 250;
if strcmpi(ext,'.tif')
    info = hdfinfo(fileRefl,'eos');
    if strcmpi(info.Grid(2).Projection.ProjCode,'geo')
        s = [info.Grid(2).Rows,info.Grid(2).Columns];
        [latlim,lonlim] = EOSHDFcoord2Lim(info.Grid(2).UpperLeft,info.Grid(2).LowerRight);
        ref = georefcells(latlim,lonlim,s,'ColumnsStartFrom','north');
        geotiffwrite(fileSCA,fSCA,cTable,ref,'TiffTags',TiffTags)
        geotiffwrite(fielSCA2,maskedSCA,cTable,ref,'TiffTags',TiffTags)
    elseif  strcmpi(info.Grid(2).Projection.ProjCode,'snsoid')
        [~,fileName,~] = fileparts(fileRefl);
        tile = fileName(18:18+5);
        [ RefMatrix,ProjectionStructure,RasterReference,GeoKeyDirectoryTag] = ...
            sinusoidProjMODtile(tile);  %#ok<ASGLU>
        geotiffwrite(fileSCA,fSCA,cTable,RasterReference.RasterReference_500m, ...
            'GeoKeyDirectoryTag', GeoKeyDirectoryTag,'TiffTags',TiffTags)
        geotiffwrite(fielSCA2,maskedSCA,cTable,RasterReference.RasterReference_500m, ...
            'GeoKeyDirectoryTag', GeoKeyDirectoryTag, 'TiffTags',TiffTags)
    end
elseif strcmpi(ext,'.hdf')
    gridName = 'MODIS_Grid_500m_2D';
    fieldName = 'QC_500m_1';
    
    gridName2 = 'MODIS_Grid_1km_2D';
    fieldName2 = 'SensorZenith_1';
    fieldName3 = 'state_1km_1';
    
    fileNameOut = 'Fractional_Snow_Cover';
    fileNameOut2 = 'Fractional_Snow_Cover_Masked';
    
    fileNameOut4 = 'Model_RMSE';
    
    fill_value = uint8(255);
    QC = hdfread(fileRefl,fieldName);
    SensorZen = hdfread(fileRefl,fieldName2);
    state1km =  hdfread(fileRefl,fieldName3);
    
    import matlab.io.hdfeos.*
    gfid = gd.open(fileRefl);
    gridID = gd.attach(gfid,gridName);
    [xdim,ydim,upLeft,lowRight] = gd.gridInfo(gridID);
    originCode = gd.originInfo(gridID);
    pixRegCode = gd.pixRegInfo(gridID);
    [projCode,zoneCode,sphereCode,projParm] = gd.projInfo(gridID);
    [compCode,compParms] = gd.compInfo(gridID,fieldName);
    [dims,ntype,dimlist] = gd.fieldInfo(gridID,fieldName);
    gd.detach(gridID);
    
     gridID = gd.attach(gfid,gridName2);
    [xdim2,ydim2,~,~] = gd.gridInfo(gridID);
    [~,~,dimlist2] = gd.fieldInfo(gridID,fieldName2);
    gd.detach(gridID);
    gd.close(gfid);
    
    % the defProj function claims:
    % If projCode is 'geo', then zoneCode, sphereCode, and projParm should be specified as [].
    if strcmpi(projCode,'geo')
        zoneCode = [];
        sphereCode = [];
        projParm = [];
    end
    
    gfid = gd.open(fileSCA,'create');
    gridID = gd.create(gfid,gridName,xdim,ydim,upLeft,lowRight);
    % note that defProj must be used before defField
    gd.defProj(gridID,projCode,zoneCode,sphereCode,projParm)
    gd.defOrigin(gridID,originCode)
    gd.defPixReg(gridID,pixRegCode)
    gd.defComp(gridID,compCode,compParms)
    
    gd.defField(gridID,fileNameOut,dimlist,dtype)
    gd.setFillValue(gridID,fileNameOut,fill_value);
    gd.writeField(gridID,fileNameOut,fSCA');
    
    gd.defField(gridID,'Fractional_Snow_Cover_NDSI',dimlist,dtype)
    gd.setFillValue(gridID,'Fractional_Snow_Cover_NDSI',fill_value);
    
    gd.defField(gridID,fileNameOut3,dimlist,dtype);
    gd.setFillValue(gridID,fileNameOut3,fill_value);
    gd.writeField(gridID,fileNameOut3,fVCA');
    
    gd.defField(gridID,fileNameOut4,dimlist,dtype);
    gd.setFillValue(gridID,fileNameOut4,fill_value);
    gd.writeField(gridID,fileNameOut4,unmixError');
    
    gd.defField(gridID,fieldName,dimlist,class(QC))
    gd.setFillValue(gridID,fieldName,cast(787410671,class(QC)));
    gd.writeField(gridID,fieldName,QC');
    
    gd.writeAttr(gridID,'modification_date',datestr(now));
    
    gd.detach(gridID);
    
    gridID = gd.create(gfid,gridName2,xdim2,ydim2,upLeft,lowRight);
    % note that defProj must be used prior to defField
    gd.defProj(gridID,projCode,zoneCode,sphereCode,projParm)
    gd.defOrigin(gridID,originCode)
    gd.defPixReg(gridID,pixRegCode)
    gd.defComp(gridID,compCode,compParms)
    
    gd.defField(gridID,fieldName2,dimlist2,class(SensorZen))
    gd.setFillValue(gridID,fieldName2,cast(787410671,class(SensorZen)));
    gd.writeField(gridID,fieldName2,SensorZen');

    gd.defField(gridID,fieldName3,dimlist2,class(state1km))
    gd.setFillValue(gridID,fieldName3,intmax(class(state1km)));
    gd.writeField(gridID,fieldName3,state1km');
    
    gd.detach(gridID);
    
    gd.close(gfid);
    
end



function tf = checkSize(A,B,varargin)
fun = @(x,y) isempty(A) || isempty(B) || isequal(size(A),size(B));
tf = fun(A,B);

if ~tf,return,end

n = nargin-2;
i = 0;
if n>0
    while tf && i<n
        i = i + 1;
        tf = tf && fun(A,varargin{i});
    end
end

function [latlim,lonlim] = EOSHDFcoord2Lim(upleft,lowright)
%We need to readjust the limits of latitude and longitude. HDF-EOS is using DMS(DDDMMMSSS.SS) format to represent degrees.
%So to calculate the lat and lon in degree, one needs to convert minutes and seconds into degrees. 
%The following is the detailed description on how to calculate the latitude and longitude range based on lowright and upleft.
% One should observe the fact that 1 minute is 60 seconds and 1 degree is 60 minutes. 

upleft_positive = ones(size(upleft));
upleft_positive(upleft<0)=-1;
lowright_positive = ones(size(lowright));
lowright_positive(lowright<0)=-1;
upleft = abs(upleft);
lowright = abs(lowright);
% First, calculate the difference of .SS between lowright and upleft.
lowright_ss= lowright*100-floor(lowright)*100;
upleft_ss = upleft*100-floor(upleft)*100;

% Then, calculate the difference of SSS between lowright and upleft.
lowright_s = mod(floor(lowright),1000);
upleft_s = mod(floor(upleft),1000);

% Then, calculate the difference of MMM between lowright and upleft.
lowright_m = mod(floor(lowright/1000),1000);
upleft_m = mod(floor(upleft/1000),1000);

% Then, calculate the difference of DDD between lowright and upleft.
lowright_d = floor(lowright/1000000);
upleft_d = floor(upleft/1000000);


ul = ((upleft_ss/100+upleft_s)/60+upleft_m)/60+upleft_d ;
lr = ((lowright_ss/100+lowright_s)/60+lowright_m)/60+lowright_d ;

ul = ul.*upleft_positive;
lr =lr.*lowright_positive;
latlim  = [lr(2),ul(2)];
lonlim = [ul(1),lr(1)];
