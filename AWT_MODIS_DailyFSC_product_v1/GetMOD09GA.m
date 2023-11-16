function X = GetMOD09GA( file, whichVariable )
% X = GetMOD09GA( file, whichVariable )
%extract variable from MOD09GA input file (or MYD09GA)
%
% Input
%   file - name of MOD09GA or MYD09GA .hdf file
%   whichVariable - choose among (case insensitive)
%       1 km spacing
%           'num_obs_1km', 'state', 'SensorZenith', 'SolarZenith','SolarAzimuth'
%           (others could be added but these are ones we use)
%       500 m spacing
%           'num_obs_500m', 'band1', 'band2', 'band3', ... 'band7',
%           'QC'
%           or 'allbands' to get bands 1 thru 7
%
% Output
%   variable scaled to appropriate values and units (the angle and band
%   variables are returned as floating point single precision, state and QC
%   are returned as structures, others as native size integers)

if iscell(file)
    file = char(file);
end
angle = false;
refl = false;
switch lower(whichVariable)
    case 'num_obs_1km'
        var = 'num_observations_1km';
    case 'state'
        var = 'state_1km_1';
    case lower('SensorZenith')
        var = 'SensorZenith_1';
        angle = true;
    case lower('SolarZenith')
        var = 'SolarZenith_1';
        angle = true;
    case lower('SolarAzimuth')
        var = 'SolarAzimuth_1';
        angle = true;
    case 'num_obs_500m'
        var = 'num_observations_500m';
    case lower('QC')
        var = 'QC_500m_1';
    case lower('allbands')
        X = GetMOD09GA(file,'band1');
        for b=2:7
            X = cat(3,X,GetMOD09GA(file,['band' num2str(b)]));
        end
        return
    otherwise
        if strncmpi(whichVariable,'band',4)
            bandNo = str2double(whichVariable(5:end));
            var = ['sur_refl_b' num2str(bandNo,'%02d') '_1'];
            refl = true;
        else
            error('Variable %s unrecognized',whichVariable)
        end
end

angleScaleDivisor = 100;
reflScaleDivisor = 10000;

iX = hdfread(file,var);
if angle
    X = single(iX)/angleScaleDivisor;
    X(X<-180) = NaN;
elseif refl
    X = single(iX)/reflScaleDivisor;
%     X(X<0) = 0;
    X(X<0) = NaN;
elseif strncmpi(whichVariable,'state',5)
    X = unpackMOD09state(iX);
elseif strncmpi(whichVariable,'qc',2)
    X = unpackMOD09QC(iX);
else
    X = iX;
end
end

function [ S ] = unpackMOD09QC( QC_500m )
%unpack the bit fields in the MOD09
%   from MOD09 User's Guide v1.3, Table 15
%   http://modis-sr.ltdri.org/products/MOD09_UserGuide_v1_3.pdf
%
% Input
%   QC_500m - QC variable read from GetMOD09GA
%
% Output
%   S - structure with the following fields
%       modland QA - 0 good, 1 not so good, 2 or 3 bad
%       band QA (vector for 7 bands) - 0 best, 7 noisy detector, 8 dead
%           detector, 9 solar Z >= 86, 10 solar Z>=85 && <86, 11 missing
%           input, 12 climatology for atmosphere, 13 out of bounds, 14 data
%           faulty, 15 not processed, deep ocean or cloud
%       atmospheric correction - true or false
%       adjacency correction - true or false

flagname = {'modland','bandQA','atmosCorr','adjCorr'};
datatype = {'uint8','uint8','logical','logical'};
nbits = [2 4 1 1];

nbands = 7;
N = QC_500m; % original 32-bit integers
for k=1:length(flagname)
    I = 2^nbits(k)-1; % fill the lower nbits(k), leave others 0
    if k~=2
        if strcmp(datatype{k},'logical')
            S.(flagname{k}) = false(size(QC_500m));
            S.(flagname{k}) = bitand(N,I)>0;
        else
            S.(flagname{k}) = zeros(size(QC_500m),datatype{k});
            S.(flagname{k}) = cast(bitand(N,I),datatype{k});
        end
        N = bitshift(N,-nbits(k));
    else
        bq = zeros([size(QC_500m) nbands],datatype{k});
        for b=1:nbands
            bq(:,:,b) = cast(bitand(N,I),datatype{k});
            N = bitshift(N,-nbits(k));
        end
        S.(flagname{k}) = bq;
    end
end

end
function [ S ] = unpackMOD09state( state_1km )
%unpack the bit fields in the MOD09 state_1km variable
%   from MOD09 User's Guide v1.4, Table 13
%   https://landweb.modaps.eosdis.nasa.gov/QA_WWW/forPage/user_guide/MOD09_UserGuide_v1.4.pdf
%
% Input
%   state_1km - state variable read from GetMOD09GA
% 
% Output
%   S - structure with the following fields
%       cloud - 0 clear, 1 cloudy, 2 mixed, 3 not sure
%       cloud shadow - true or false
%       landwater - 0 shallow ocean, 1 land, 2 coastline, 3 shallow lake,
%           4 ephemeral water, 5 deep lake, 6 continental ocean, 7 ocean
%       aerosol - 0 climatology, 1 low, 2 average, 3 high
%       cirrus - 0 none, 1 small, 2 average, 3 high
%       Intcloudflag - internal alg. cloud flag, 1 true
%       Intfireflag - internal alg. fire flag, 1 true
%       MOD35snow - MOD35 snow/ice flag, 1 true
%       cloudadj - pixel is adjacent to cloud, 1 true
%       saltpan - pixel is a salt pan, 1 true
%       intsnow - internal alg. snow mask, 1 true


flagname = {'cloud','cloudshadow','landwater','aerosol','cirrus','Intcloudflag',...
    'Intfireflag','MOD35snow','cloudadj','saltpan','Intsnow'};
datatype = {'uint8','logical','uint8','uint8','uint8','logical','logical',...
    'logical','logical','logical','logical'};

nbits = [2 1 3 2 2 1 1 1 1 1 1];

for k=1:length(flagname)
    if strcmp(datatype{k},'logical')
        S.(flagname{k}) = false(size(state_1km));
    else
        S.(flagname{k}) = zeros(size(state_1km),datatype{k});
    end
end

N = state_1km; % original 16-bit integers
for k=1:length(flagname)
    I = 2^nbits(k)-1; % fill the lower nbits(k), leave others 0
    if strcmp(datatype{k},'logical')
        S.(flagname{k}) = false(size(state_1km));
        switch flagname{k}
            case 'cloudshadow'
                pos=2;
            case 'Intcloudflag'
                pos=10;    
            case 'Intfireflag'
                pos=11;
            case 'MOD35snow'
                pos=12;
            case 'cloudadj'
                pos=13;
            case 'saltpan'
                pos=14;
            case 'Intsnow'
                pos=15;
            otherwise
            error('flagname %s not recognized as logical',flagname{k})
        end
                S.(flagname{k}) = bitget(state_1km,pos+1)>0;
%     end
                %         if strcmp(flagname{k},'cloudshadow')
%             S.(flagname{k}) = bitand(N,I)>0;
%         elseif strcmp(flagname{k},'MOD35snow') % MOD35 snow flag
%             S.(flagname{k}) = bitget(c,13)>0;
%         end
    else
        S.(flagname{k}) = zeros(size(state_1km),datatype{k});
        S.(flagname{k}) = cast(bitand(N,I),datatype{k});
    end
    N = bitshift(N,-nbits(k));
end
end