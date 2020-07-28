% Author: Tobias Gerken, tug15@psu.edu
% Version: 1.0 (2019_05_08)
% Compile final merged MetaData files

% This code reads in metadata flag ICART files compiled by Lily Cambell (ljc42@psu.edu) and
% 1) add Stop_UTC, and Mid_UTC timestamp colums
% 2) adds ManeuverFlag and Manuever Flag QC column
% 3) Replaces missing values (-9999) with -99999 for consistency.

clear all; close all;

startup % sets path and basedirectory

% Loop over files 

ICARTDir = './Data/ACTAMERICA_MetaDataFlags/4Campaigns/Final_R0/' ;
ICARTSDir = {'Summer16/', 'Winter17/', 'Fall17/', 'Spring18/'} ;

OutDir = './Data/ACTAMERICA_MetaDataFlags/4Campaigns/Final_R0_New/' 
MFDir = './Data/ACTAMERICA_ManeuverFlag/NC/'

MFNme = 'ManeuverFlag_' ;
MFExt = '.csv' ;

nanVal  = -999999;

%loop over days in campaign
for i = ICARTSDir
    fIC = dir([ICARTDir i{1} '*.ict']);
    
    % Make Output Dir
    mkdir([OutDir i{1}])
    
    for f = 1:length(fIC)
        fid = fopen([fIC(f).folder '/' fIC(f).name]);
        
        disp(['File:'  fIC(f).name]);
        
        % get number of lines in header
        clear line
        line{1} = fgetl(fid);
        C = textscan(line{1}, '%d%d', 'delim',',');
        
        % read and write header line by line and write to file
        
        Outname = fIC(f).name;
        Outname=strrep(Outname, 'ACTAMERICA_flag_metadata', 'ACTAMERICA-flag-metadata');
        fOut = fopen([OutDir i{1} Outname], 'w');
        %fprintf(fOut,[line '\n']);
        for l = 1: C{1}-1
            line{l+1} = fgetl(fid);
            %fprintf(fOut, [line '\n']);
        end
        
        % Modify lines: 
        % Line 2
        line{2} = strrep(line{2}, 'Gerkin', 'Gerken') ;
        % Line 10-12
        line{10} = num2str((str2num(line{10})-1)) ;
        dummy = strsplit(line{11},', ');
        line{11} = strjoin(dummy(1:end-1),', ');
        dummy = strsplit(line{12},', ');
        line{12} = strjoin(dummy(1:end-1),', ');
        % Select header line and replace Height with GPS_ALT
        ll = find(contains(line, 'Start_UTC,'));
        line{ll(end)} =  strrep(line{ll(end)}, 'Height', 'GPS_ALT');
        % Insert Special Comment Line 
        line = [line(1:23), {'0'}, line(24:end)] ;
        
        
        C = textscan(line{1}, '%d%d', 'delim',',');
        
        line{1} = sprintf('%d, %d', C{1}+1, C{2});
        
        for ll = 1:length(line)
            fprintf(fOut, [line{ll} '\n']);
        end
        % now read data portion of file
        D = textscan(fid, '%f%f%f%f%f%f%f%f', 'delim',',');      
        fclose(fid);
       
        % read ManeuverFlag File 
        p = lower(fIC(f).name(26:29)); % get plane 
        d = fIC(f).name(31:38);
        
        disp(['Plane: ' p ', Day: ' d]);
        
        fMF = dir([MFDir MFNme p '_' d '*' MFExt]);
        
        MF = csvread([fMF.folder '/' fMF.name] ,1) ;
        
        % Create Data Out Array
        DataOut = nan(length(D{1}),12);
        DataOut(:,1) = D{1};
        DataOut(:,2) = D{1}+5;
        DataOut(:,3) = D{1}+2.5;
        DataOut(:,4) = D{2};
        DataOut(:,5) = D{3};
        DataOut(:,6) = D{4};
        DataOut(:,7) = D{5};
        DataOut(:,8) = D{6};
        DataOut(:,9) = D{7};
        DataOut(:,10) = D{8};
        
        % match timestamps between MF and Metadata
        [~, iD, iMF] = intersect(squeeze(DataOut(:,1)), squeeze(MF(:,1))) ;
        if iMF ~=length(DataOut)
            warning('Data not same length for %s: %s', p, d)
            if length(iMF) ~= length(DataOut)
                warning('WARNING - MF Data missing for %s:', fIC(f).name(26:end-4))
            else
                warning('OK - MF Data complete:', fIC(f).name(26:end-4))
            end
        end
        DataOut(iD,11) = MF(iMF,6);
        DataOut(iD,12) = MF(iMF,7);
        
        % Replace missing values (-99999) with -999999 for consistency in
        % our data 
        ind = false(size(DataOut));
        ind(:,8:12) = (DataOut(:,8:12) < -8999);        
        DataOut(ind) = nanVal;
        % Replace NaN with -99999 for consistency
        DataOut(isnan(DataOut)) = nanVal;
        
        % now write to file 
        for li = 1:length(DataOut)
            outstr = sprintf('%d, %d, %7.1f, %d, %11.8f, %12.8f, %6.4f, %d, %d, %d, %d, %d \n', DataOut(li,:));
            outstr = strrep(outstr,'-9999.0000,','-999999');
            outstr = strrep(outstr,'-999999.00000000,','-999999');
            fprintf(fOut, outstr) ;
        end
        fclose(fOut);
    end
    
end