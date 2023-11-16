%% ********************************************************************* %%
% AWT MODIS FSC Cloud Gap Filled product code
% 0-100:FSC, 237:Water，250:Cloud，253/255:Nodata
% Code by Fangbo Pan and Gongxue Wang, November 15, 2023
% jiang@bnu.edu.cn,panfb@mail.bnu.edu.cn
%% ********************************************************************* %%
function [ ] = blendMOYD(...
    filelistMOD,filelistMYD,filelistMCD,productVersion,productID)

narginchk(3,5)
if nargin<4
    productVersion = 'C6';
end
if ~iscell(filelistMOD),filelistMOD = cellstr(filelistMOD);end
if ~iscell(filelistMYD),filelistMYD = cellstr(filelistMYD);end
if ~iscell(filelistMCD),filelistMCD = cellstr(filelistMCD);end

if nargin<5, productID = 'MOD10A1';end
productID = validatestring(productID,...
    {'MODAGE','MOD10A1'},mfilename,'productID',2);

assert(isequal(size(filelistMOD),size(filelistMYD), size(filelistMCD)),...
    'input and output filelist must have same size')

import matlab.io.hdfeos.*
import matlab.io.hdf4.*

txtStr = '0=Terra MODIS, 1=Aqua MODIS';
n = 0;
for i=1:length(filelistMOD)
    disp(['making file: ',filelistMCD{i}])
    
    structMOD = getMOD10A1(filelistMOD{i},productVersion);
    structMYD = getMOD10A1(filelistMYD{i},productVersion);
    
    n = n + 1;
    if n==1
        info = hdfinfo(filelistMOD{i},'eos');
        fields = fieldnames(structMOD);
    end
    
    switch productID
        case 'MOD10A1'
            index = indBlendWithBestQA(structMOD.(fields{1}),structMOD.(fields{2}),...
                structMYD.(fields{1}),structMYD.(fields{2}));
        case 'MODAGE'
            fieldNames = {'Fractional_Snow_Cover','Fractional_Snow_Cover_Masked',...
                'Fractional_Vegetation_Cover','Model_RMSE',...
                'SensorZenith_1','QC_500m_1','state_1km_1'};
            tf = isfield(structMOD,fieldNames);
            if ~tf(2)
                warning('the fractional snow cover layer not found')
                continue
            else
                ind = find(strcmpi(fieldNames(2),fields));
                TSnow = structMOD.(fields{ind(1)});
                ASnow = structMYD.(fields{ind(1)});
            end
            if ~tf(5)
                TSAZ500m  = [];
                ASAZ500m  = [];
            else
                ind = find(strcmpi(fieldNames(5),fields));
                TSAZ500m = imresize(structMOD.(fields{ind(1)}),size(TSnow),'nearest');
                ASAZ500m = imresize(structMYD.(fields{ind(1)}),size(ASnow),'nearest');
            end
            [~,index] = MCDSnow(TSnow, TSAZ500m , ASnow, ASAZ500m);
    end
    
    % writing
    [pathstr,~,~] = fileparts(filelistMCD{i});
    
    if ~exist(pathstr,'dir'),mkdir(pathstr);end
    if exist(filelistMCD{i},'file'),delete(filelistMCD{i});end
    [status,message,messageId] = copyfile(filelistMOD{i},filelistMCD{i}, 'f'); %#ok<ASGLU>
    [status,msg,msgID] = fileattrib(filelistMCD{i},'+w'); %#ok<ASGLU>
    
    gfid = gd.open(filelistMCD{i},'rdwr');
    tf = true;
    for g = 1:numel(info.Grid)
        gridID = gd.attach(gfid,info.Grid(g).Name);
        gridFieldNames = cellfun(@(x) x,{info.Grid(g).DataFields.Name},'un',false);
        try
            for f = 1:numel(fields)
                if ~any(strcmpi(fields{f},gridFieldNames)),continue,end
                temp = structMOD.(fields{f});
                temp0 = structMYD.(fields{f});
                ind = imresize(index,size(temp),'nearest');
                temp(ind) = temp0(ind);
                gd.writeField(gridID,fields{f},temp');
                clear temp temp0
            end
            s = size(index);
            if tf && s(1) == info.Grid(g).Rows && s(2) == info.Grid(g).Columns
                dimlist = {'XDim','YDim'};
                dtype = 'uint8';
                gd.defField(gridID,'input_ID',dimlist,dtype)
                gd.writeField(gridID,'input_ID',cast(index,dtype)');
                tf = false;
            end
            gd.detach(gridID);
        catch
            warning('failed to write file')
            gd.detach(gridID);
        end
    end
    gd.close(gfid);
    clear structMOD structMYD
    try
        sdID = sd.start(filelistMCD{i},'rdwr');
        idx = sd.nameToIndex(sdID,'input_ID');
        sdsID = sd.select(sdID,idx);
        sd.setAttr(sdsID,'description',txtStr);
        sd.endAccess(sdsID);
        sd.close(sdID);
    catch
        warning('failed to write file')
        sd.close(sdID);
    end
    [status,msg,msgID] = fileattrib(filelistMCD{i},'-w'); %#ok<ASGLU>
end


function index = indBlendWithBestQA(MOD,MODQA,MYD,MYDQA)
% indBlendWithBestQA(MOD,MODQA,MYD,MYDQA)
% it returns the index of MYD to replace data in MOD
keyCloud = 250;
keySnow = 0:100;
keyMiss = 200;
keyNoDec = 201;
keyNight = 211;
keyLake = 237;
keyOcean = 239;
keySat = 254;
keyFill = 255;

% where Terra does not has valid detection but Aqua has
ix1 = MOD>keySnow(end) & MOD~=keyLake & MOD~=keyOcean;
ix2 = MYD<=keySnow(end) | MYD==keyLake | MYD==keyOcean | MYD==keyCloud;

% where both detect snow but Terra has worse QA
ix3 = MOD<=keySnow(end) & MYD<=keySnow(end);
ix4 = MODQA > MYDQA;

index = (ix1&ix2) | (ix3&ix4);

%
% % where both detect snow and have same QA
% ix1 = MOD<=keySnow(end) & MYD<=keySnow(end);
% ix2 = MODQA==MYDQA;
% temp = round(double(MOD(ix1&ix2)+MYD(ix1&ix2))/2);
% MCD(ix1&ix2) = cast(temp,class(MCD));
function dSnow = MCDSnow( TSnow, TSAZ, ASnow, ASAZ)
% dSnow = MCDSnow( TSnow, TSAZ, ASnow, ASAZ)
% to blend fractional snow cover from Terra MODIS and Aqua MODIS
% 
% 

scanAng = 55;
thr = 40;

cloudflag = 250;
snowflag=0:100;
gapflag = 254;

func = @(x) isnumeric(x) && isreal(x);

assert(func(TSnow),'1st input must be numeric and real')
assert(func(TSAZ),'2nd input must be numeric and real')
assert(func(ASnow),'3rd input must be numeric and real')
assert(func(ASAZ),'4th input must be numeric and real')

valSize = isequal(size(TSnow),size(TSAZ),size(ASnow),size(ASAZ));
assert(valSize,'all inputs must have same size')

% zenith angle should be in the range [0 180]
% and be smaller than threshold 55 degree.
ix = TSAZ <0 | TSAZ > scanAng;
TSAZ(ix)=NaN;

ix = ASAZ <0 | ASAZ > scanAng;
ASAZ(ix)=NaN;


dSnow = TSnow;

ix = (TSnow == cloudflag | TSnow ==gapflag) &...
    ASnow ~= cloudflag & ASnow ~= gapflag;
dSnow(ix) = ASnow(ix);

ix = TSnow >=snowflag(1) & TSnow <=snowflag(end)  & ...
    ASnow >=snowflag(1) & ASnow <=snowflag(end)  & ...
    TSAZ > thr  & TSAZ > ASAZ ;

dSnow(ix) = ASnow(ix);