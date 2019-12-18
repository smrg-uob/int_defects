function [ outputStruct ] = read_dat( inputStruct )
%read_dat.m
%Harry Coules 2013
%   04/2014 - Minor modifications for improved string find performance.
%
%DESCRIPTION
%This function scans an Abaqus/Standard .dat file and extracts output
%variables. It is possible to extract most Abaqus output variables from the
%.dat file, but only to a relatively low precision (5 S.F.). For more
%precise output, it is necessary to use the .fil file.
%
%This function largely depends on the .dat file being in a specific format,
%i.e. for appropriate output requests to have been made in the Abaqus .inp
%file prior to executing the calculation.
%
%INPUT ARGUMENTS (Input structure can contain:)
%   inputStruct.filename - Name of the .dat file.
%   inputStruct.stepsToRead - Vector to define model step(s) to read output
%       data for.
%   inputStruct.incsToRead - Vector to define model increment(s) to read 
%       output data for.
%
%OUTPUT ARGUMENTS (Output structure can contain:)
%   outputStruct.alpha - Backstress tensor (element-wise at centroid)
%   outputStruct.ee - Elastic strain (element-wise at centroid)
%   outputStruct.ePl - Plastic strain (element-wise at centroid)
%   outputStruct.eL - Logarithmic (total) strain (element-wise at centroid)
%   outputStruct.eN - Nominal (total) strain (element-wise at centroid)
%   outputStruct.eThe - Thermal strain (element-wise at centroid)
%   outputStruct.evol - Element volumes (element-wise at centroid)
%   outputStruct.peeq - Equivalent plastic strain (element-wise at centroid)
%   outputStruct.s - Stress (element-wise at centroid)
%
%   outputStruct.sIP - Stress (at element integration points)
%
%   outputStruct.pener - Energy dissipated by plasticity per unit volume (element-wise)
%   outputStruct.sener - Elastic strain energy density (element-wise)
%
%   outputStruct.rf - Nodal forces (node-wise)
%   outputStruct.u - Nodal displacements (node-wise)
%   outputStruct.T - Nodal temperature (node-wise)
%
%   outputStruct.k - K Factor information (evaluated on requested contours)
%   outputStruct.jInt - J-integral (evaluated on requested contours)
%   outputStruct.tStress - T-stress (evaluated on requested contours)
%
%NOTES
%   1. This function always returns a structure containing arrays, which
%       themselves contain the extracted data. There is a single array for
%       each type of output (eg. stress, nodal force etc.).
%   2. When multiple steps and increments are requested using
%       inputStruct.stepsToRead and inputStruct.incsToRead, the results for
%       all steps and increments are read, undistinguished, into the single
%       output array for that data type. read_dat cannot currently
%       distinguish between data read from different increments in a single
%       call. If you do need to distinguidh between results from multiple
%       different increments, it is necessary to call read_dat multiple
%       times, each time specifying a different (single) step and/or
%       increment.
%   3. **NEED TO EXPLAIN BEHAVIOUR WHEN MULTIPLE ELEMENT TYPES ARE PRESENT**
%   4. Preallocation of output arrays is used. Currently, the array
%       preallocation length is 5e6. If the amount of ouput data for a
%       single data type is more than 5e6 then this function is likely
%       to ehibit a severe drop in performance. No warnings are given.
%
%EXAMPLE ABAQUS OUTPUT REQUESTS - note that this list of examples is not exhaustive.
%Elastic strain (EE)
%   *El Print, Elset = ExtractSetComplete, Position = Centroidal, Frequency=99999999
%   coord, EE
%
%Logarithmic (total) strain (LE)
%   *El Print, Elset = ExtractSetComplete, Position = Centroidal, Frequency=99999999
%   coord, LE
%
%Nominal (total) strain (NE)
%   *El Print, Elset = ExtractSetComplete, Position = Centroidal, Frequency=99999999
%   coord, NE
%
%Plastic strain (PE)
%   *El Print, Elset = ExtractSetComplete, Position = Centroidal, Frequency=99999999
%   PE
%
%Stress (S)
%   *El Print, Elset = ExtractSetComplete, Position = Centroidal, Frequency=99999999
%   coord, S
%
%% Preliminary checks, and set default values if input struct is incomplete
%Which step(s) to read output for. Sanitise this input.
if isfield(inputStruct, 'stepsToRead')
    if isnumeric(inputStruct.stepsToRead)
        if ~any(inputStruct.stepsToRead<0)
            if ~any(mod(inputStruct.stepsToRead,1)~=0)
                stepsToRead = inputStruct.stepsToRead;
            else
                error('inputStruct.stepsToRead contains non-integer value(s)')
            end
        else
            error('inputStruct.stepsToRead contains negative value(s).')
        end
    else
        error('inputStruct.stepsToRead should be a numeric.')
    end
else
    warning('inputStruct.stepsToRead not set in input to read_dat.m. Reading output for all steps present in the .dat file.')
    stepsToRead = 0;
end

%Which increment(s) to read output for. Sanitise this input.
if isfield(inputStruct, 'incsToRead')
    if isnumeric(inputStruct.incsToRead)
        if ~any(inputStruct.incsToRead<0)
            if ~any(mod(inputStruct.incsToRead,1)~=0)
                incsToRead = inputStruct.incsToRead;
            else
                error('inputStruct.incsToRead contains non-integer value(s)')
            end
        else
            error('inputStruct.incsToRead contains negative value(s).')
        end
    else
        error('inputStruct.incsToRead should be a numeric.')
    end
else
    warning('inputStruct.incsToRead not set in input to read_dat.m. Reading output for all increments in specified the step(s).')
    incsToRead = 0;
end

%If filename is empty in input structure, give an error.
if ~isfield(inputStruct,'filename')
    error('No filename defined for the .dat file to be read.');
else
    filename = inputStruct.filename;
end

%Create cell arrays of strings from stepsToRead and incsToRead - used to
%determine whether output for each step/inc should be read
f1 = @(x) sprintf(' %.f ',x);
stepsToReadCell = cellfun(f1, num2cell(stepsToRead),'UniformOutput', false);
incsToReadCell = cellfun(f1, num2cell(incsToRead),'UniformOutput', false);

%% Read in the lines of the .dat file
fid = fopen(filename,'r','n','UTF-8');
C = textscan(fid, '%s', Inf, 'Delimiter','\n'); C = C{1};
fclose(fid);

alpha = []; %Elemental back-stress tensor
ee = [];  %Elastic strain
ePl = []; %Plastic strain
eL = [];    %Logarithmic (total) strain
eN = [];    %Nominal (total) strain
eThe = [];  %Thermal strain
evol = [];  %Elemental volume
peeq = []; %Elemental PEEQ values
pener = []; %Elemental energy dissipation by plasticity w.r.t volume
rf = []; %Nodal forces
s = []; %Elemental stresses
sIP = []; %Stresses at element integration points
sener = []; %Elemental strain energy density
u = []; %Nodal displacements
T = []; %Nodal temperature

jInt_1 = [];    %J-integral (column 1 - contour no.s)
jInt_2 = [];    %J-integral (column 2 - values)

tStress_1 = [];    %T-stress (column 1 - contour no.s)
tStress_2 = [];    %T-stress (column 2 - values)

kContours = []; %Contour no.s over which K is evaluated
k1 = []; %Mode I SIF
k2 = []; %Mode II SIF
k3 = []; %Mode III SIF
kMTSDir = [];   %Propagation direction (via maximum tangential stress criterion)
jFromK = [];  %J calculated from SIFs

outputFlag = 'NONE';    %When header lines corresponding to an output type are read, this will be set.
outputType = '';
elemOutputType = '';

stepFlag = false;
incFlag = false;

%Find strings indicating step and increment headers and footers
stepStrFlag = ~cellfun('isempty',strfind(C,'S T E P'));                  
analysisStrFlag = ~cellfun('isempty',strfind(C,'A N A L Y S I S'));
incrementStrFlag = ~cellfun('isempty',strfind(C,'INCREMENT'));
summaryStrFlag = ~cellfun('isempty',strfind(C,'SUMMARY'));

%Scan the .dat file line-wise
for k = 1:length(C)
    %Scan lines
    %For models with multiple steps, make sure you're getting the output from the correct step:
    if stepStrFlag(k) && analysisStrFlag(k)
        if stepsToRead ~= 0
            stepFlag = false;
            for k1 = 1:length(stepsToReadCell)
                if ~isempty(strfind(C{k},stepsToReadCell{k1}))
                    stepFlag = true;
                end
            end
        else
            stepFlag = true;
        end
    end
    
    if incrementStrFlag(k) && summaryStrFlag(k)
        if incsToRead ~= 0
            incFlag = false;
            for k1 = 1:length(incsToReadCell)
                if ~isempty(strfind(C{k},incsToReadCell{k1}))
                    incFlag = true;
                end
            end
        else
            incFlag = true;
        end
    end
    
    %Note that these lines will only be used on the correct step and inc
    if stepFlag && incFlag
        %Determine what type of output will be given in subsequent lines by
        %searching for header strings.
        if strcmpi(C{k},'E L E M E N T   O U T P U T');
            outputFlag = 'ELEM';
        elseif strcmpi(C{k},'N O D E   O U T P U T');
            outputFlag = 'NODE';
        elseif strcmpi(C{k},'J - I N T E G R A L   E S T I M A T E S');
            outputFlag = 'J';
        elseif strcmpi(C{k},'K   F A C T O R       E S T I M A T E S');
            outputFlag = 'K';
        elseif strcmpi(C{k},'T - S T R E S S   E S T I M A T E S');
            outputFlag = 'T';
        elseif strcmpi(C{k},'THE ANALYSIS HAS BEEN COMPLETED');
            outputFlag = 'END';
        elseif ~isempty(strfind(C{k}, 'THE FOLLOWING TABLE IS PRINTED AT THE CENTROID OF THE ELEMENT'))
            elemOutputType = 'CENTROID';
        elseif ~isempty(strfind(C{k}, 'THE FOLLOWING TABLE IS PRINTED AT THE INTEGRATION POINTS'))
            elemOutputType = 'IP';
        elseif ~isempty(strfind(C{k}, 'MAXIMUM')) || ~isempty(strfind(C{k}, 'MINIMUM')) || ~isempty(strfind(C{k}, 'PROCESS'))
            outputType = '';
            elemOutputType = '';    %Reset elemOutputType if you've got to the end of the table
        elseif ~isempty(strfind(C{k}, 'ALPHA11'));
            outputType = 'ALPHA';
        elseif ~isempty(strfind(C{k}, 'EVOL'));
            outputType = 'EVOL';
        elseif ~isempty(strfind(C{k}, 'EE11'));
            outputType = 'EE';
        elseif ~isempty(strfind(C{k}, 'LE11'));
            outputType = 'LE';
        elseif ~isempty(strfind(C{k}, 'NE11'));
            outputType = 'NE';
        elseif ~isempty(strfind(C{k}, 'PE11'));
            outputType = 'PE';
        elseif ~isempty(strfind(C{k}, 'PEEQ'));
            outputType = 'PEEQ';
        elseif ~isempty(strfind(C{k}, 'PENER'));
            outputType = 'PENER';
        elseif ~isempty(strfind(C{k}, 'RF1'));
            outputType = 'RF';
        elseif ~isempty(strfind(C{k}, 'S11'));
            outputType = 'S';
        elseif ~isempty(strfind(C{k}, 'SENER'));
            outputType = 'SENER';
        elseif ~isempty(strfind(C{k}, 'THE11'));
            outputType = 'THE';
        elseif ~isempty(strfind(C{k}, 'U1'));
            outputType = 'U';
        elseif ~isempty(strfind(C{k}, 'NT')) && isempty(strfind(C{k}, 'INCREMENT')) && isempty(strfind(C{k}, 'CONTROL'));
            outputType = 'NT';
        else
            %Extract values from lines
            if ~isempty(strfind(C{k}, 'OR'))
                if strcmpi(elemOutputType,'IP')
                    line = sscanf(C{k},'%u %u OR %f %f %f %f %f %f %f %f %f')';
                else
                    line = sscanf(C{k},'%u OR %f %f %f %f %f %f %f %f %f')';
                end
            else
                line = sscanf(C{k},'%u %f %f %f %f %f %f %f %f %f')';
            end
            %Now see if these values are anything we're wanting to scan for
            if strcmp(outputFlag,'ELEM')   %Element output
                if strcmpi(elemOutputType,'CENTROID')     %IF THE OUTPUT IS GIVEN AT ELEMENT CENTROIDS:
                    if length(line) == 2    %Will have written two only if the line has only the element no. and PEEQ/PENER/SENER
                        if strcmp(outputType,'PEEQ');   %PEEQ only
                            if isempty(peeq)
                                peeq = zeros(5e6,2);    %Preallocate array
                                peeqCount = 1;
                            else
                                peeqCount = peeqCount+1;
                            end
                            peeq(peeqCount,:) = line;
                        elseif strcmp(outputType,'PENER');  %PENER only
                            if isempty(pener)
                                pener = zeros(5e6,2);    %Preallocate array
                                penerCount = 1;
                            else
                                penerCount = penerCount+1;
                            end
                            pener(penerCount,:) = line;
                        elseif strcmp(outputType,'SENER');  %SENER only
                            if isempty(sener)
                                sener = zeros(5e6,2);    %Preallocate array
                                senerCount = 1;
                            else
                                senerCount = senerCount+1;
                            end
                            sener(senerCount,:) = line;
                        end
                        
                    elseif length(line) == 6    %Will have written six if there are stresses/elastic or thermal strains and coordinates (2D plane stress case)
                        if strcmp(outputType,'ALPHA');  %2D, coords and ALPHA
                            if isempty(alpha)
                                alpha = zeros(5e6,6);    %Preallocate array
                                alphaCount = 1;
                                outputStruct.dimensions = 2;
                            else
                                alphaCount = alphaCount+1;
                            end
                            alpha(alphaCount,:) = line;
                        elseif strcmp(outputType,'S');  %2D, coords and S
                            if isempty(s)
                                s = zeros(5e6,6);    %Preallocate array
                                sCount = 1;
                                outputStruct.dimensions = 2;
                            else
                                sCount = sCount+1;
                            end
                            s(sCount,:) = line;
                        elseif strcmp(outputType,'EE'); %2D, coords and EE
                            if isempty(ee)
                                ee = zeros(5e6,6);    %Preallocate array
                                eeCount = 1;
                                outputStruct.dimensions = 2;
                            else
                                eeCount = eeCount+1;
                            end
                            ee(eeCount,:) = line;
                        elseif strcmp(outputType,'LE'); %2D, coords and LE
                            if isempty(eL)
                                eL = zeros(5e6,6);    %Preallocate array
                                eLCount = 1;
                                outputStruct.dimensions = 2;
                            else
                                eLCount = eLCount+1;
                            end
                            eL(eLCount,:) = line;
                        elseif strcmp(outputType,'NE'); %2D, coords and NE
                            if isempty(eN)
                                eN = zeros(5e6,6);    %Preallocate array
                                eNCount = 1;
                                outputStruct.dimensions = 2;
                            else
                                eNCount = eNCount+1;
                            end
                            eN(eNCount,:) = line;
                        elseif strcmp(outputType,'PE'); %2D, PE only
                            if isempty(ePl)
                                ePl = zeros(5e6,4);    %Preallocate array
                                ePlCount = 1;
                                outputStruct.dimensions = 2;
                            else
                                ePlCount = ePlCount+1;
                            end
                            ePl(ePlCount,:) = line(1:4); %Note that only the element no.s and individual plastic strain components are written to ePl. PEEQ and PEMAG are ignored.
                        elseif strcmp(outputType,'THE'); %2D, coords and THE
                            if isempty(eThe)
                                eThe = zeros(5e6,6);    %Preallocate array
                                eTheCount = 1;
                                outputStruct.dimensions = 2;
                            else
                                eTheCount = eTheCount+1;
                            end
                            eThe(eTheCount,:) = line;
                        end
                        
                    elseif length(line) == 7    %Will have written seven if there are stresses/elastic or thermal strains and coordinates (2D plane strain case)
                        if strcmp(outputType,'ALPHA');  %2D plane strain, coords and ALPHA
                            if isempty(alpha)
                                alpha = zeros(5e6,7);    %Preallocate array
                                alphaCount = 1;
                                outputStruct.dimensions = 2;
                            else
                                alphaCount = alphaCount+1;
                            end
                            alpha(alphaCount,:) = line;
                        elseif strcmp(outputType,'S');  %2D plane strain, coords and S
                            if isempty(s)
                                s = zeros(5e6,7);    %Preallocate array
                                sCount = 1;
                                outputStruct.dimensions = 2;
                            else
                                sCount = sCount+1;
                            end
                            s(sCount,:) = line;
                        elseif strcmp(outputType,'EE'); %2D plane strain, coords and EE
                            if isempty(ee)
                                ee = zeros(5e6,7);    %Preallocate array
                                eeCount = 1;
                                outputStruct.dimensions = 2;
                            else
                                eeCount = eeCount+1;
                            end
                            ee(eeCount,:) = line;
                        elseif strcmp(outputType,'LE'); %2D plane strain, coords and LE
                            if isempty(eL)
                                eL = zeros(5e6,7);    %Preallocate array
                                eLCount = 1;
                                outputStruct.dimensions = 2;
                            else
                                eLCount = eLCount+1;
                            end
                            eL(eLCount,:) = line;
                        elseif strcmp(outputType,'NE'); %2D plane strain, coords and NE
                            if isempty(eN)
                                eN = zeros(5e6,7);    %Preallocate array
                                eNCount = 1;
                                outputStruct.dimensions = 2;
                            else
                                eNCount = eNCount+1;
                            end
                            eN(eNCount,:) = line;
                        elseif strcmp(outputType,'PE'); %2D plane strain, PE only
                            if isempty(ePl)
                                ePl = zeros(5e6,5);    %Preallocate array
                                ePlCount = 1;
                                outputStruct.dimensions = 2;
                            else
                                ePlCount = ePlCount+1;
                            end
                            ePl(ePlCount,:) = line(1:5); %Note that only the element no.s and individual plastic strain components are written to ePl. PEEQ and PEMAG are ignored.
                        elseif strcmp(outputType,'THE'); %2D plane strain, coords and THE
                            if isempty(eThe)
                                eThe = zeros(5e6,7);    %Preallocate array
                                eTheCount = 1;
                                outputStruct.dimensions = 2;
                            else
                                eTheCount = eTheCount+1;
                            end
                            eThe(eTheCount,:) = line;
                        end
                        
                    elseif length(line) == 8
                        if strcmp(outputType,'PE'); %2D plane stress, coords and PE
                            if isempty(ePl)
                                ePl = zeros(5e6,8);    %Preallocate array
                                ePlCount = 1;
                                outputStruct.dimensions = 2;
                            else
                                ePlCount = ePlCount+1;
                            end
                            ePl(ePlCount,:) = line;
                        end
                        
                    elseif length(line) == 9
                        if strcmp(outputType,'PE'); %3D, PE only
                            if isempty(ePl)
                                ePl = zeros(5e6,7);    %Preallocate array
                                ePlCount = 1;
                                outputStruct.dimensions = 3;
                            else
                                ePlCount = ePlCount+1;
                            end
                            ePl(ePlCount,:) = line(1:7); %Note that only the element no.s and individual plastic strain components are written to ePl. PEEQ and PEMAG are ignored.
                        end
                        
                    elseif length(line) == 10    %Will have written ten if there are stresses/elastic or thermal strains and coordinates (3D case)
                        if strcmp(outputType,'ALPHA');  %3D, coords and ALPHA
                            if isempty(alpha)
                                alpha = zeros(5e6,10);    %Preallocate array
                                alphaCount = 1;
                                outputStruct.dimensions = 3;
                            else
                                alphaCount = alphaCount+1;
                            end
                            alpha(alphaCount,:) = line;
                        elseif strcmp(outputType,'S');  %3D, coords and S
                            if isempty(s)
                                s = zeros(5e6,10);    %Preallocate array
                                sCount = 1;
                                outputStruct.dimensions = 3;
                            else
                                sCount = sCount+1;
                            end
                            s(sCount,:) = line;
                        elseif strcmp(outputType,'EE'); %3D, coords and EE
                            if isempty(ee)
                                ee = zeros(5e6,10);    %Preallocate array
                                eeCount = 1;
                                outputStruct.dimensions = 3;
                            else
                                eeCount = eeCount+1;
                            end
                            ee(eeCount,:) = line;
                        elseif strcmp(outputType,'LE'); %3D, coords and LE
                            if isempty(eL)
                                eL = zeros(5e6,10);    %Preallocate array
                                eLCount = 1;
                                outputStruct.dimensions = 3;
                            else
                                eLCount = eLCount+1;
                            end
                            eL(eLCount,:) = line;
                        elseif strcmp(outputType,'NE'); %3D, coords and NE
                            if isempty(eN)
                                eN = zeros(5e6,10);    %Preallocate array
                                eNCount = 1;
                                outputStruct.dimensions = 3;
                            else
                                eNCount = eNCount+1;
                            end
                            eN(eNCount,:) = line;
                        elseif strcmp(outputType,'THE'); %3D, coords and THE
                            if isempty(eThe)
                                eThe = zeros(5e6,10);    %Preallocate array
                                eTheCount = 1;
                                outputStruct.dimensions = 3;
                            else
                                eTheCount = eTheCount+1;
                            end
                            eThe(eTheCount,:) = line;
                        end
                    end
                    
                elseif strcmpi(elemOutputType,'IP')     %IF THE OUTPUT IS GIVEN AT ELEMENT INTEGRATION POINTS:
                    %%%DEVELOPMENT%%% Note that this part currently only
                    %%%supports output of stresses at the integration points.
                    if length(line) == 7    %Will have written seven if there are stresses and coordinates (2D plane stress case)
                        if strcmp(outputType,'S');  %2D, coords and S
                            if isempty(sIP)
                                sIP = zeros(5e6,7);    %Preallocate array
                                sIPCount = 1;
                                outputStruct.dimensions = 2;
                            else
                                sIPCount = sIPCount+1;
                            end
                            sIP(sIPCount,:) = line;
                        end
                        
                    elseif length(line) == 8    %Will have written eight if there are stresses and coordinates (2D plane strain case)
                        if strcmp(outputType,'S');  %2D, coords and S
                            if isempty(sIP)
                                sIP = zeros(5e6,8);    %Preallocate array
                                sIPCount = 1;
                                outputStruct.dimensions = 2;
                            else
                                sIPCount = sIPCount+1;
                            end
                            sIP(sIPCount,:) = line;
                        end
                        
                    elseif length(line) == 11    %Will have written eleven if there are stresses and coordinates (3D case)
                        if strcmp(outputType,'S');  %3D, coords and S
                            if isempty(sIP)
                                sIP = zeros(5e6,11);    %Preallocate array
                                sIPCount = 1;
                                outputStruct.dimensions = 3;
                            else
                                sIPCount = sIPCount+1;
                            end
                            sIP(sIPCount,:) = line;
                        end
                    end
                    
                elseif strcmpi(elemOutputType,'')   %IF THE OUTPUT IS NOT SPECIFIED TO BE AT EITHER CENTROID OR INTEGRATION POINTS
                    if length(line) == 2    %Will have written two only if the line has only the element no. and EVOL
                        if strcmp(outputType,'EVOL');   %EVOL only
                            if isempty(evol)
                                evol = zeros(5e6,2);    %Preallocate array
                                evolCount = 1;
                            else
                                evolCount = evolCount+1;
                            end
                            evol(evolCount,:) = line;
                        end
                    end
                end
                
            elseif strcmp(outputFlag,'NODE')    %Nodal output
                if length(line) == 4    %Will have written two only if the line has only the node no., coord and NT for a 2D model, or node no. and RF.
                    if strcmp(outputType,'NT');
                        if isempty(T)
                            T = zeros(5e6,5);    %Preallocate array
                            TCount = 1;
                        else
                            TCount = TCount+1;
                        end
                        T(TCount,:) = line;
                    elseif strcmp(outputType,'RF');
                        if isempty(rf)
                            rf = zeros(5e6,4);    %Preallocate array
                            rfCount = 1;
                        else
                            rfCount = rfCount+1;
                        end
                        rf(rfCount,:) = line;
                    end
                    
                elseif length(line) == 5    %Will have written five only if it's node no., coords and NT for a 3D model or node no., coords and U for a 2D model.
                    if strcmp(outputType,'NT');
                        if isempty(T)
                            T = zeros(5e6,5);    %Preallocate array
                            TCount = 1;
                        else
                            TCount = TCount+1;
                        end
                        T(TCount,:) = line;
                    elseif strcmp(outputType,'U');
                        if isempty(u)
                            u = zeros(5e6,5);    %Preallocate array
                            uCount = 1;
                        else
                            uCount = uCount+1;
                        end
                        u(uCount,:) = line;
                    end
                    
                elseif length(line) == 7    %Will have written five only if it's node no., coords and U for a 2D model
                    if strcmp(outputType,'U');
                        if isempty(u)
                            u = zeros(5e6,7);    %Preallocate array
                            uCount = 1;
                        else
                            uCount = uCount+1;
                        end
                        u(uCount,:) = line;
                    end
                end
                
            elseif strcmp(outputFlag,'J')   %J-Integral output
                if strcmpi(C{k},'LABELS REFERENCED IN THE ABOVE TABLE');
                    outputFlag = 'NONE';   %This prevents the code from reading data below the final data line
                else
                    %Special regexp for J-integral output format
                    if ~isempty(strfind(C{k}, '- '))
                        jLine = sscanf(C{k},' -%u- %f %f %f %f %f')';
                        if length(jLine) > 1
                            jLine = jLine(2:end);
                            dataSwitch = 1; %If there's a -%u- in the line, and at least one other number, we've moved on to the J-Integral data lines
                        end
                    else
                        jLine = sscanf(C{k},'%f %f %f %f %f')';
                    end
                    
                    if length(jLine) ~=0
                        if ~exist('prevJLineLength','var')
                            prevJLineLength = Inf;
                            dataSwitch = 0;
                        end
                        
                        if dataSwitch == 0;
                            jInt_1 = [jInt_1;jLine'];
                        else
                            jInt_2 = [jInt_2;jLine'];
                        end
                        
                        prevJLineLength = length(jLine);
                    end
                end
                
            elseif strcmp(outputFlag,'T')    %T-stress output
                if strcmpi(C{k},'LABELS REFERENCED IN THE ABOVE TABLE');
                    outputFlag = 'NONE';   %This prevents the code from reading data below the final data line
                else
                    %Special regexp for T-stress output format (note it's the same as for the J-integral)
                    if ~isempty(strfind(C{k}, '- '))
                        tLine = sscanf(C{k},' -%u- %f %f %f %f %f')';
                        if length(tLine) > 1
                            tLine = tLine(2:end);
                            dataSwitch = 1; %If there's a -%u- in the line, and at least one other number, we've moved on to the T-stress data lines
                        end
                    else
                        tLine = sscanf(C{k},'%f %f %f %f %f')';
                    end
                    
                    if length(tLine) ~=0
                        if ~exist('prevTLineLength','var')
                            prevTLineLength = Inf;
                            dataSwitch = 0;
                        end
                        
                        if dataSwitch == 0;
                            tStress_1 = [tStress_1;tLine'];
                        else
                            tStress_2 = [tStress_2;tLine'];
                        end
                        
                        prevTLineLength = length(tLine);
                    end
                end
                
            elseif strcmp(outputFlag,'K')    %K factors output
                if strcmpi(C{k},'LABELS REFERENCED IN THE ABOVE TABLE');
                    outputFlag = 'NONE';   %This prevents the code from reading data below the final data line
                else
                    %Special regexp for K factor output
                    if strfind(C{k}, 'K1:')
                        k1Line = sscanf(C{k},' K1: %f %f %f %f %f');
                        if isempty(k1Line)
                            k1Line = sscanf(C{k},' -%u- K1: %f %f %f %f %f');
                            k1Line(1) = [];
                        end
                        k1 = [k1;k1Line];
                    elseif strfind(C{k}, 'K2:')
                        k2Line = sscanf(C{k},' K2: %f %f %f %f %f');
                        k2 = [k2;k2Line];
                    elseif strfind(C{k}, 'K3:')
                        k3Line = sscanf(C{k},' K3: %f %f %f %f %f');
                        k3 = [k3;k3Line];
                    elseif strfind(C{k}, 'MTS   DIRECTION (DEG):')
                        kMTSDirLine = sscanf(C{k},' MTS   DIRECTION (DEG): %f %f %f %f %f');
                        kMTSDir = [kMTSDir;kMTSDirLine];
                    elseif strfind(C{k}, 'J from Ks:')
                        jFromKLine = sscanf(C{k},' J from Ks: %f %f %f %f %f');
                        jFromK = [jFromK;jFromKLine];
                    else kContoursLine = sscanf(C{k},'%f %f %f %f %f');
                        kContours = [kContours;kContoursLine];
                    end
                end
            end
        end
    end
end

%% Organise output into the output struct
%Ensure all output data is sorted by element or node number (i.e. the 1st
%column in each array)
if ~isempty(alpha)
    alpha = alpha(any(alpha~=0,2),:);
    alpha = sortrows(alpha,1);
end

if ~isempty(ee)
    ee = ee(any(ee~=0,2),:);
    ee = sortrows(ee,1);
end

if ~isempty(evol)
    evol = evol(any(ee~=0,2),:);
    evol = sortrows(evol,1);
end

if ~isempty(eL)
    eL = eL(any(eL~=0,2),:);
    eL = sortrows(eL,1);
end

if ~isempty(eN)
    eN = eN(any(eN~=0,2),:);
    eN = sortrows(eN,1);
end

if ~isempty(ePl)
    ePl = ePl(any(ePl~=0,2),:);
    ePl = sortrows(ePl,1);
end

if ~isempty(eThe)
    eThe = eThe(any(eThe~=0,2),:);
    eThe = sortrows(eThe,1);
end

if ~isempty(peeq)
    peeq = peeq(any(peeq~=0,2),:);
    peeq = sortrows(peeq,1);
end

if ~isempty(pener)
    pener = pener(any(pener~=0,2),:);
    pener = sortrows(pener,1);
end

if ~isempty(rf)
    rf = rf(any(rf~=0,2),:);
    rf = sortrows(rf,1);
end

if ~isempty(sener)
    sener = sener(any(sener~=0,2),:);
    sener = sortrows(sener,1);
end

if ~isempty(s)
    s = s(any(s~=0,2),:);
    s = sortrows(s,1);
end

if ~isempty(sIP)
    sIP = sIP(any(sIP~=0,2),:);
    sIP = sortrows(sIP,[1,2]);
end

if ~isempty(T)
    T = T(any(T~=0,2),:);
    T = sortrows(T,[1,2]);
end

if ~isempty(u)
    u = u(any(u~=0,2),:);
    u = sortrows(u,1);
end

%Flag for problems with contour integral output. This may occur if eg.
%there are multiple cracks in the model that can't be distinguished.
contourOutputErrorFlag = 0;
if length(jInt_1)~=length(jInt_2)
    contourOutputErrorFlag = 1;
elseif length(kContours)~=length(k1)
    contourOutputErrorFlag = 1;
end

%J-integral and T-stress
if contourOutputErrorFlag
    jInt = [];
    tStress = [];
else
    jInt = [jInt_1,jInt_2];
    tStress = [tStress_1,tStress_2];
end

%Stress intensity factors
if isempty(k3);     %For 2D problems, there will be no Mode III SIF calculated. Create a dummy array full of zeros.
    k3 = zeros(size(kContours));
end
if contourOutputErrorFlag
    kFactors = [];
else
    kFactors = [kContours,k1,k2,k3,kMTSDir,jFromK];
end

%Construct a structure containing the extracted data
outputStruct.alpha = alpha;
outputStruct.ee = ee;
outputStruct.eL = eL;
outputStruct.eN = eN;
outputStruct.evol = evol;
outputStruct.ePl = ePl;
outputStruct.eThe = eThe;
outputStruct.jInt = jInt;
outputStruct.k = kFactors;
outputStruct.peeq = peeq;
outputStruct.pener = pener;
outputStruct.rf = rf;
outputStruct.sener = sener;
outputStruct.s = s;
outputStruct.sIP = sIP;
outputStruct.tStress = tStress;
outputStruct.T = T;
outputStruct.u = u;

end

