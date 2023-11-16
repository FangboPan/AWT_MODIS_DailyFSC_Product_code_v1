%% ********************************************************************* %%
% AWT MODIS FSC Cloud Gap Filled product code
% 0-100:FSC, 237:Water，250:Cloud，253/255:Nodata
% Code by Fangbo Pan and Gongxue Wang, November 15, 2023
% jiang@bnu.edu.cn,panfb@mail.bnu.edu.cn
%% ********************************************************************* %%
addpath(genpath('/forMODIS/CodeShare'),'-end')

dates = compose('%03d',1:365);
tiles = {'h23v03','h23v04','h23v05','h24v04','h24v05','h24v06','h25v04','h25v05','h25v06',...
   'h26v05','h26v06','h27v06' };

idir = {'/MODIS_Daily_Ref/MOD09GA_C6/2000/','/MODIS_Daily_Ref/MYD09GA_C6/2000/'};
odir = {'/MODIS_Daily_Ref/MOD_result/MODAGE_C6/2000/','/MODIS_Daily_Ref/MOD_result/MYDAGE_C6/2000/'};

MODID = {'MOD09GA', 'MYD09GA','MCD09GA'};
snowID = {'MODAGE', 'MYDAGE','MCDAGE'};
%% main
for d=1:366
% M0D fSCA
    odir1 = [odir{1},dates{d},filesep];
    mkdir(odir1);
    if ~exist(odir{1},'dir'),mkdir(odir1);end
    temp = [idir{1},dates{d},filesep];
    files = dir([temp,'MOD09GA*.hdf']);
    fileMOD09GA = cellfun(@fullfile,{files.folder},{files.name},'un',false)';
    fileMODfSCA = fileMOD09GA;
    for i = 1:numel(fileMOD09GA)
        [~,name,~] = fileparts(fileMOD09GA{i});
        fileMODfSCA{i} = [odir1,name,'.fSCA.hdf'];
    end

    parfor i=1:numel(fileMOD09GA)
        if exist(fileMODfSCA{i},'file'),delete(fileMODfSCA{i});end
        try
            MODISfSCAwithMESMA(fileMOD09GA{i},fileMODfSCA{i})
        catch
            continue
        end
    end
    
%% MYD fSCA
    odir2 = [odir{2},dates{d},filesep];
    mkdir(odir2);
    if ~exist(odir{2},'dir'),mkdir(odir2);end
    temp2 = [idir{2},dates{d},filesep];
    files2 = dir([temp2,'MYD09GA*.hdf']);
    fileMYD09GA = cellfun(@fullfile,{files2.folder},{files2.name},'un',false)';
    fileMYDfSCA = fileMYD09GA;
    for j = 1:numel(fileMYD09GA)
        [~,name,~] = fileparts(fileMYD09GA{j});
        fileMYDfSCA{j} = [odir2,name,'.fSCA.hdf'];
    end
    fileMYD09A1 ={''};
    parfor j=1:numel(fileMYD09GA)
        if exist(fileMYDfSCA{j},'file'),delete(fileMYDfSCA{j});end
        try
            MODISfSCAwithMESMA(fileMYD09GA{j},fileMYD09A1,fileMYDfSCA{j})
        catch
            continue
        end
    end
%% MOD+MYD fSCA
   for m=1:numel(tiles)
    MCDAWfile = strrep(fileMOD09GA,'MOD','MCD');
    [pathStr,~,~]=fileparts(MCDAWfile{m});
    mkdir(pathStr);
     files = dir([odir1,MODID{1},'*',tiles{m},'*.fSCA.hdf']);
     filesTemp = cellfun(@fullfile,{files.folder},{files.name},'un',false)';
     files2 = dir([odir2,MODID{2},'*',tiles{m},'*.fSCA.hdf']);
     filesTemp2 = cellfun(@fullfile,{files2.folder},{files2.name},'un',false)';
       if ~isempty(files)
            if ~isempty(files2)
                 blendMOYD(filesTemp,filesTemp2,MCDAWfile{m},'C6','MODAGE') 
             else
                filesTemp2 = filesTemp;
                blendMOYD(filesTemp,filesTemp2,MCDAWfile{m},'C6','MODAGE') 
                str2=[MODID{2} '.' '2000' dates{d} '.' tiles{m}];
                fid=fopen('/mnt/lustre/users/jiang/MODIS_Daily_Ref/MOD_result/MCDAGE_C6/Filled.txt','a');
                fprintf(fid,'%s\n',str2);
                fclose(fid);
             end
       else
           if ~isempty(files2)
                 filesTemp = filesTemp2;
                 blendMOYD(filesTemp,filesTemp2,MCDAWfile{m},'C6','MODAGE')
                 str=[MODID{1} '.' '2000' dates{d} '.' tiles{m}];
                fid=fopen('/MODIS_Daily_Ref/MOD_result/MCDAGE_C6/Filled.txt','a');
                fprintf(fid,'%s\n',str);
                fclose(fid);
           else
                str=[MODID{1} '.' '2000' dates{d} '.' tiles{m}];
                str2=[MODID{2} '.' '2000' dates{d} '.' tiles{m}];
                fid=fopen('/MODIS_Daily_Ref/MOD_result/MCDAGE_C6/Filled.txt','a');
                fprintf(fid,'%s\n',str);
                fprintf(fid,'%s\n',str2);
                fclose(fid);
           end
       end   
   end
end