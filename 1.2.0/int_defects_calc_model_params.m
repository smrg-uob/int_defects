function [ modelParamStruct, varargout ] = int_defects_calc_model_params( naturalParamStruct )
%int_defects_calc_model_params.m
%Harry Coules 2015
%
%DESCRIPTION
%This function takes in a structure containing geometric information for a 
%pair of semi-elliptical edge cracks, given in terms of 'natural'
%parameters (mainly ratios between geometric lengths). It returns a
%structure containing geometric parameters in the format required for
%writing into an Abaqus python script.
%
%INPUT ARGUMENTS
%   naturalParamStruct - Structure of parameters in 'natural' form, i.e. in
%       terms of ratios between the various characteristic lengths.
%
%OUTPUT ARGUMENTS
%   modelParamStruct - Structure containing parameters in absolute form,
%       i.e. the actual lengths etc. to be provided to Abaqus in order to
%       build a .inp file.
%   OPTIONAL OUTPUT ARGUMENT
%   validityComment - A string containing information about the validity of
%       the model parameters.
%
%NOTES
% - If naturalParamStruct.aB1OveraB2 is set to NaN, the same aspect ratio
%   will be used for Crack B as for Crack A.
%
%% Set validFlag and valitityComment
%These will subsequently be changes if it is determined that the parameters requested would not result in a valid model.
validFlag = true;
validityComment = 'Okay';

%% Dynamically unpack naturalParamStruct
%Loop over main fields and sub-fields, creating variables from the
%fields of the structure. Note that this will throw an error if any fields
%have the same name.
fieldNames = fieldnames(naturalParamStruct);
for k1 = 1:length(fieldNames)
    if isstruct(eval(['naturalParamStruct.',fieldNames{k1},';']))
        subFieldNames = fieldnames(eval(['naturalParamStruct.',fieldNames{k1},';']));
        for k2 = 1:length(subFieldNames)
            if exist(subFieldNames{k2},'var')
                error(['Error unpacking naturalParamStruct. Duplicate field name: ',subFieldNames{k2}]);
            else
                eval([subFieldNames{k2},' = naturalParamStruct.',fieldNames{k1},'.',subFieldNames{k2},';']);
            end
        end
    else
        if exist(fieldNames{k1},'var')
            error(['Error unpacking naturalParamStruct. Duplicate field name: ',fieldNames{k1}]);
        else
            eval([fieldNames{k1},' = naturalParamStruct.',fieldNames{k1},';']);
        end
    end
end

%% If necessary, make initial modifications to crack geometry
%This is useful to, for example, create a set of models which represent
%assesment procedure re-characterisation of the defects run in an earlier
%set of models
if exist('modifyNaturalParamStructStr','var')
    if strcmpi(modifyNaturalParamStructStr,'')
        %If modifyNaturalParamStructStr is empty, do nothing.
    elseif strcmpi(modifyNaturalParamStructStr,'DNVGL-RP-F108')
        subsurfaceAFlag = false;
        
        %Determine new parameters
        aA1Overb = (2*naturalParamStruct.aA1Overb)+(naturalParamStruct.d2Overa*naturalParamStruct.aA1Overb);
        aA2Overb = ((2*naturalParamStruct.aA1Overb./naturalParamStruct.aA1OveraA2)+(naturalParamStruct.d2Overa*naturalParamStruct.aA1Overb))/2;
        aA1OveraA2 = aA1Overb./aA2Overb;
        clear d2Overa
        
        disp('Modified crack geometry parameters according to DNVGL-PR-F108 embedded defect re-characterisation.');
    else
        error('Unexpected string in: naturalParamStruct.modifyNaturalParamStructStr.')
    end
else
    modifyNaturalParamStructStr = '';
end

%% Check geometry definitions in naturalParamStruct
%Wall thickness
if exist('b','var')
    if exist('ri','var') && exist('riOverro','var')
        bFromPipeParams = (ri./riOverro)-ri;
        if bFromPipeParams ~= b
            error('Multiple conflicting definitions of wall thickness exist (1: b, 2: ri and riOverro).')
        end
    end
else
    if exist('ri','var') && exist('riOverro','var')
        b = (ri./riOverro)-ri;
    else
        error('Wall thickness is not full defined. Either b or ri and riOverro is required.');
    end
end

%Check definition of size (in depth direction) of Crack A
if ~exist('aA1Overb','var')
    error('Depth of Crack A is not defined. See naturalParamStruct.geometryParams.');
end

%Check definition of Crack A aspect ratio
if ~exist('aA1OveraA2','var')
    error('Aspect ratio of Crack A is not defined. See naturalParamStruct.geometryParams.');
end

%Check definition of size (in depth direction) of Crack B
if exist('aB1OveraA1','var') && ~exist('aB1Overb','var')
    %Do nothing
elseif exist('aB1Overb','var') && ~exist('aB1OveraA1','var')
    %Do nothing
else
    error('Size of Crack B in the depth direction not properly defined. See naturalParamStruct.geometryParams.');
end

%Check definition of Crack B aspect ratio
% - If naturalParamStruct.aB1OveraB2 equals NaN, then use the same aspect ratio for Crack B as Crack A has.
if length(aB1OveraB2)==1
    if isnan(aB1OveraB2)
        if ~singleCrackFlag
            disp('aB1OveraB2 is NaN. Using the same aspect ratio for Crack B as for Crack A.');
            aB1OveraB2 = aA1OveraA2;
        end
    end
end

%Check distance between the two ellipses.
if ~singleCrackFlag
    if sum(ismember({'dOverb','dOveraA1','dOverMax_aA1_2aB1','twoaA2OverdCents','overlapRatio','sOverb','normalisedOffset'},who))>1
        error('Two conflicting crack spacing definitions are present.');
    end
end

%Check depth of any subsurface crack.
if subsurfaceAFlag || subsurfaceBFlag
    if sum(ismember({'d2Overa','d2Overb','SOvera','SOverb','SOverc','SOverSqrtac','normalisedOffsetY'},who))>1
        error('Two conflicting embedded crack depth definitions are present.');
    end
    if exist('d2Overa','var')
        if d2Overa < 0
            error('d2Overa cannot be negative.')
        end
    end
    if exist('d2Overb','var')
        if d2Overb < 0
            error('d2Overb cannot be negative.')
        end
    end
end

%Check plate depth
if ~exist('plDepth','var')
    plDepth = 100*b;   %Default is 100 times the thickness
end

%Check plate half-width
if sum([exist('plHalfWidth','var'),exist('plHalfWidthOveraA2','var'),exist('plHalfWidthOveraA2Plusb','var')])>1
    error('Two conflicting plate width definitions are present.');
elseif ~exist('plHalfWidth','var') && ~exist('plHalfWidthOveraA2','var') && ~exist('plHalfWidthOveraA2Plusb','var')
    plHalfWidth = 50*b;   %Default is 50 times the thickness
elseif exist('plHalfWidthOveraA2','var')
    plHalfWidth = plHalfWidthOveraA2*(aA1Overb*b/aA1OveraA2);
elseif exist('plHalfWidthOveraA2Plusb','var')
    plHalfWidth = plHalfWidthOveraA2Plusb*((aA1Overb*b/aA1OveraA2)+b);
end

%Determine bounding box magnitude (defined conservatively)
boundBoxMag = 2*max([plDepth,plHalfWidth]);

%% Define model geometry parameters
%Crack ellipse radii and crack centre locations...
%Crack dimensions
aA1 = aA1Overb*b;
aA2 = aA1/aA1OveraA2;

if exist('aB1OveraA1','var')
    aB1 = aB1OveraA1*aA1;
elseif exist('aB1Overb','var')
    aB1 = aB1Overb*b;
end
aB2 = aB1/aB1OveraB2;

%Determine y coordinate of crack centres
if subsurfaceAFlag
    if exist('d2Overa','var')
        yA = (d2Overa*aA1)+aA1;
    elseif exist('d2Overb','var')
        yA = (d2Overb*b)+aA1;
    elseif exist('normalisedOffsetY','var')   %Normalised offset in y.
        yA = b*(0.5-normalisedOffsetY);
    else
        error('When subsurfaceAFlag = true, the depth of Crack A must be defined by d2Overa or d2Overb.')
    end
    yB = NaN;
elseif subsurfaceBFlag
    if exist('d2Overa','var')       %Distance from shallowest point to y=0 plane, normalised to aA1.
        yB = (d2Overa*aA1)+aB1;
    elseif exist('d2Overb','var')   %Distance from shallowest point to y=0 plane, normalised to b.
        yB = (d2Overb*b)+aB1;
    elseif exist('SOvera','var')    %Distance in y between cracks, normalised to aA1.
        yB = (SOvera*aA1)+aA1+aB1;
    elseif exist('SOverb','var')    %Distance in y between cracks, normalised to b.
        yB = (SOverb*b)+aA1+aB1;
    elseif exist('SOverc','var')    %Distance in y between cracks, normalised to aA2.
        yB = (SOverc*aA2)+aA1+aB1;
    elseif exist('SOverSqrtac','var')   %Distance in y between cracks, normalised to sqrt(aA1*aA2).
        yB = (SOverSqrtac*sqrt(aA1*aA2))+aA1+aB1;
    elseif exist('normalisedOffsetY','var')   %Normalised offset in y.
        yB = normalisedOffsetY*(aA1+aB1);
    end
    yA = 0;
else
    yA = 0;
    if oppositeFlag
        yB = b;
    else
        yB = 0;
    end
end

%Crack separation distance in x-direction
if exist('dOverb','var')    %Distance in x between crack edges, normalised to b.
    if ~oppositeFlag
        if dOverb < 0
            error('For cracks emanating from the same side of the plate, dOverb cannot be negative.')
        end
    end
    if isinf(dOverb)
        disp('dOverb is given as Inf in naturalParamStruct. Assuming cracks are remote from each other.');
        dOverb = max([aA2*10,aB2*10,10]);
    end
    d = dOverb.*b;
elseif exist('dOveraA1','var')    %Distance in x between crack edges, normalised to aA1.
    if ~oppositeFlag
        if dOveraA1 < 0
            error('For cracks emanating from the same side of the plate, dOveraA1 cannot be negative.')
        end
    end
    if isinf(dOveraA1)
        disp('dOveraA1 is given as Inf in naturalParamStruct. Assuming cracks are remote from each other.');
        dOveraA1 = max([aA2*10,aB2*10,10]);
    end
    d = dOveraA1.*aA1;
elseif exist('dOverMax_aA1_2aB1','var')    %Distance in x between crack edges, defined as a ratio of max([aA1,2*aB1])
    if dOverMax_aA1_2aB1 < 0
        error('dOverMax_aA1_2aB1 cannot be negative.')
    end
    if isinf(dOverMax_aA1_2aB1)
        error('dOverMax_aA1_2aB1 is given as Inf in naturalParamStruct. Invalid value.');
    end
    if subsurfaceAFlag || ~subsurfaceBFlag
        error('dOverMax_aA1_2aB1 may only be used when Crack A is at the surface and Crack B is sub-surface.')
    end
    d = dOverMax_aA1_2aB1.*max([aA1,2*aB1]);
elseif exist('twoaA2OverdCents','var')  %Used for comparison with Sethuraman et al.
    if ~oppositeFlag
        if twoaA2OverdCents > 1 && ~isinf(twoaA2OverdCents)
            error('For cracks emanating from the same side of the plate, twoaA2OverdCents  cannot be greater than 1.')
        end
    end
    if isinf(twoaA2OverdCents)
        disp('twoaA2OverdCents is given as Inf in naturalParamStruct. Assuming cracks are remote from each other.');
        twoaA2OverdCents = min([(2*aA2)/(10*aA2+aA2+aB2),(2*aA2)/(10*aB2+aA2+aB2),(2*aA2)/(10+aA2+aB2)]);
    end
    d = ((2*aA2)/twoaA2OverdCents)-(aA2+aB2);
elseif exist('overlapRatio','var')  %At overlapRatio = 0, near-side crack tips are at the same location in x. At overlapRatio = 1, the cracks are aligned
    if ~oppositeFlag
        error('overlapRatio should not be used to determine the distance between cracks on the same side of the plate.')
    end
    if isinf(overlapRatio)
        disp('overlapRatio is given as Inf in naturalParamStruct. Assuming cracks are remote from each other.');
        overlapRatio = min([-10,-10/(aA2+aB2)]);
    end
    d = -overlapRatio*(aA2+aB2);
elseif exist('sOverb','var')
    if oppositeFlag
        if subsurfaceBFlag
            %When a subsurface crack is considered, we need to modify the input
            %to account for this...
            [ d ] = int_defects_d_from_sOverb( aA1, aA2, aB1, aB2, yB, sOverb, 1e-4 );   %Determine d from sOverb
            warning(['Determining d from sOverb for the case where Crack B is a subsurface defect using\n',...
                'int_defects_d_from_sOverb.m. I have not yet validated this method.']);
        else
            %Otherwise proceed as normal.
            [ d ] = int_defects_d_from_sOverb( aA1, aA2, aB1, aB2, b, sOverb, 1e-4 );   %Determine d from sOverb
        end
        if isnan(d)
            modelParamStruct.validFlag = false;
            validityComment = 'Required separation for opposite-side cracks is too small.';
        end
    else
        error('sOverb should not be used to define the inter-crack spacing unless oppositeFlag = true.')
    end
elseif exist('normalisedOffset','var')
    if ~oppositeFlag
        if normalisedOffset < 1
            error('For cracks emanating from the same side of the plate, normalisedOffset cannot be less than 1.')
        end
    end
    if isinf('normalisedOffset')
        disp('normalisedOffset is given as Inf in naturalParamStruct. Assuming cracks are remote from each other.');
        normalisedOffset = 11;  %The distance between crack edges is 10 times the mean crack half-width.
    end
    d = (normalisedOffset-1)*(aA2+aB2);
else
    error('Unable to determine the distance between ellipses.')
end

%x coordinates of crack centres
if singleCrackFlag
    xA = 0; %x coordinate of centre of Crack A for a single crack model
else
    xA = aA2 + d/2;
end
xB = -(aB2 + d/2);

%Notes:
% - Crack A is in +ive x direction, Crack B in -ive x direction.
% - x=0 is half-way between the near-side extremes of the cracks.
% - In certain cases where cracks on opposite sides of the plate overlap in
%   x, the centres of both cracks may be at positive locations in x. This
%   occurs when d < 0 && aB2 < -d/2.

%% Define crack tip element zone sizes
%Number of elements along the crack tip line
if ~exist('noTipElemsA','var')
    warning('Number of elements along the Crack A tip line not defined. Defaulting to 80.');
    noTipElemsA = 80;
end
if ~exist('noTipElemsB','var')
    if singleCrackFlag
        noTipElemsB = NaN;
    else
        warning('Number of elements along the Crack B tip line not defined. Defaulting to 80.');
        noTipElemsB = 80;
    end
end

%Number of elements along the radial direction in the outer crack tip region
if ~exist('noTipElemsRadial2A','var')
    warning('Number of elements in radial direction in the Crack A outer region not defined, defaulting to 5.');
    noTipElemsRadial2A = 5;
end
if ~exist('noTipElemsRadial2B','var')
    if singleCrackFlag
        noTipElemsRadial2B = NaN;
    else
        warning('Number of elements in radial direction in the Crack B outer region not defined, defaulting to 5.');
        noTipElemsRadial2B = 5;
    end
end

%Crack tip element zone size - flag for automatic mode
if ~exist('varShellSizeFlag','var')
    varShellSizeFlag = true;
    warning('Field varShellSizeFlag not found in naturalParamStruct. Assuming varTubeShellFlag = true.');
end

%Define crack tip element zone size
if varShellSizeFlag
    %Notes:
    % - I have changed the default radius slightly to accomodate sOverb =
    %   0.2 when meshing opposite cracks. It is now aA1/11.
    % - For high aspect ratio cracks aA1OveraA2 > 2 or aB1OveraB2 > 2,
    %   TubeShellLarge is unlikely to sweep successfully with a
    %   normal-sized profile radius. The default profile radius for these
    %   high aspect ratio cases is therefore aA1/21 or aB1/21. This may
    %   come at a penalty in terms of accuracy.
    %Profile radius for TubeShellLargeA
    if aA1OveraA2 > 2
        rpA2 = aA1/21;  %Reduced default radius for high aspect ratio cracks
    else
        rpA2 = aA1/11;
    end
    
    %Profile radius for TubeShellLargeB
    if aB1OveraB2 > 2
        rpB2 = aB1/21;   %Reduced default radius for high aspect ratio cracks
    else
        rpB2 = aB1/11;
    end
else
    if ~exist('rpA2','var')
        error('Explicit crack tip zone sizing has been selected (varShellSizeFlag=false) but crack tip zone size rpA2 not defined.');
    end
    if singleCrackFlag
        if ~exist('rpB2','var')
            rpB2 = NaN;
        end
    end
end

%For non-linear modes where J-integral output is needed, more accurate
%results will be given with a larger crack tip mesh zone size.
if isfield(naturalParamStruct.outputRequests,'contourTypeA')
    if strcmpi(naturalParamStruct.outputRequests.contourTypeA,'J_INTEGRAL')
        rpA2 = 2*rpA2;
    end
end
if isfield(naturalParamStruct.outputRequests,'contourTypeB')
    if strcmpi(naturalParamStruct.outputRequests.contourTypeB,'J_INTEGRAL')
        rpB2 = 2*rpB2;
    end
end

rpA1 = rpA2/(noTipElemsRadial2A+1);    %Profile radius for TubeShellSmallA
rpB1 = rpB2/(noTipElemsRadial2B+1);    %Profile radius for TubeShellSmallB
%Note: Currently, this is set up so that the elements in the rpA1 are
%the same size as those in rpA2.

%% Check the crack tip zone sizes. PART 1 - Intersection with the other crack's crack tip zone
%nB. This code currently attempts to adjust the crack tip zone size only
%for intersection of the two crack tip zones, and only when you are using
%variable rather than statically-defined crack tip zone sizes.
crackTipElemSizeWarning = false;
if oppositeFlag
    %When the cracks are opposite each other, need to use
    %calc_ellipses_approach to check for crack tip element zone
    %intersection.
    %NOTES
    % - This method should work even when subsurfaceBFlag = true. Note that
    % the coordinate yB is used in the input to calc_ellipses_approach.
    % - This method is only used when Crack B is subsurface. Crack A is 
    % only ever subsurface then there is only a single crack.
    [ intersectionFlag,~,~,~,sMin,~,~] = calc_ellipses_approach( aA1+rpA2, aA2+rpA2, aB1+rpB2, aB2+rpB2, d-rpA2-rpB2, yB );   %Check for crack tip element zone intersection. Note modification to d for larger ellipses.
    if intersectionFlag || sMin<rpA1 || sMin<rpB1
        if varShellSizeFlag
            warning('The crack tip element zones of the two cracks will intersect or approach each other too closely! Attempting adjustment...')
            %Halve all crack tip zone radii
            rpA1 = rpA1/2;   %Profile radius for TubeShellSmallA
            rpB1 = rpB1/2;   %Profile radius for TubeShellSmallB
            rpA2 = rpA2/2;   %Profile radius for TubeShellLargeA
            rpB2 = rpB2/2;   %Profile radius for TubeShellLargeB
            [ intersectionFlag,~,~,~,sMin,~,~] = calc_ellipses_approach( aA1+rpA2, aA2+rpA2, aB1+rpB2, aB2+rpB2, d-rpA2-rpB2, yB );   %Re-check for intersection
            if intersectionFlag || sMin<rpA1 || sMin<rpB1
                warning('Crack tip element zone size adjustment unsuccessful, the cracks are too close together.')
                validFlag = false;
                validityComment = 'Cracks are too close to accomodate (variable) crack tip element zones.';
            else
                warning('Adjustment successful. Crack tip element sizes halved. Results should be treated with caution.')
                validFlag = true;
                crackTipElemSizeWarning = true;
            end
        else
            warning('The crack tip element zones of the two cracks will intersect.')
            validFlag = false;
            validityComment = 'Cracks are too close to accomodate (static) crack tip element zones.';
        end
    end
else
    %When the cracks are on the same side, it is simpler to check for
    %crack tip element zone intersection
    if d <= (rpA2 + rpB2)
        if varShellSizeFlag
            warning('The crack tip element zones of the two cracks will intersect! Attempting adjustment...')
            %Halve all crack tip zone radii
            rpA1 = rpA1/2;   %Profile radius for TubeShellSmallA
            rpB1 = rpB1/2;   %Profile radius for TubeShellSmallB
            rpA2 = rpA2/2;   %Profile radius for TubeShellLargeA
            rpB2 = rpB2/2;   %Profile radius for TubeShellLargeB
            if d <= (rpA2 + rpB2)
                warning('Crack tip element zone size adjustment unsuccessful, the cracks are too close together.')
                validFlag = false;
                validityComment = 'Cracks are too close to accomodate (variable) crack tip element zones.';
            else
                warning('Adjustment successful. Crack tip element sizes halved. Results should be treated with caution.')
                validFlag = true;
                crackTipElemSizeWarning = true;
            end
        else
            warning('The crack tip element zones of the two cracks will intersect.')
            validFlag = false;
            validityComment = 'Cracks are too close to accomodate (static) crack tip element zones.';
        end
    end
end

%% Check the crack tip zone sizes. PART 2 - Intersection with edges of the plate
%nB. No attempt is made to adjust crack tip zone sizes if there is this
%type of intersection problem.

%Crack A
if subsurfaceAFlag
    %Criteria for use when Crack A is a subsurface crack.
    if yA <= (aA1 + rpA2)
        warning('The crack tip element zone of Crack A will extend beyond the front face of the plate.')
        validFlag = false;
        validityComment = 'Crack tip element zone of Crack A will extend beyond the front face of the plate.';
    end
    if (b-yA) <= (aA1 + rpA2)
        warning('The crack tip element zone of Crack A will extend beyond the back face of the plate.')
        validFlag = false;
        validityComment = 'Crack tip element zone of Crack A will extend beyond the back face of the plate.';
    end
else
    %Criteria for use when Crack A is a surface crack.
    if aA1 <= rpA2
        warning('The crack tip element zone of Crack A will extend beyond the front face of the plate.')
        validFlag = false;
        validityComment = 'Crack tip element zone of Crack A will extend beyond the front face of the plate.';
    end
    if b <= (aA1 + rpA2)
        warning('The crack tip element zone of Crack A will extend beyond the back face of the plate.')
        validFlag = false;
        validityComment = 'Crack tip element zone of Crack A will extend beyond the back face of the plate.';
    end
end
%Crack B
if subsurfaceBFlag
    %Criteria for use when Crack B is a subsurface crack.
    if yB <= (aB1 + rpB2)
        warning('The crack tip element zone of Crack B will extend beyond the front face of the plate.')
        validFlag = false;
        validityComment = 'Crack tip element zone of Crack B will extend beyond the front face of the plate.';
    end
    if (b-yB) <= (aB1 + rpB2)
        warning('The crack tip element zone of Crack B will extend beyond the back face of the plate.')
        validFlag = false;
        validityComment = 'Crack tip element zone of Crack B will extend beyond the back face of the plate.';
    end
else
    %Criteria for use when Crack B is a surface crack.
    if aB1 <= rpB2
        warning('The crack tip element zone of Crack B will extend beyond the front face of the plate.')
        validFlag = false;
        validityComment = 'Crack tip element zone of Crack B will extend beyond the front face of the plate.';
    end
    if b <= (aB1 + rpB2)
        warning('The crack tip element zone of Crack B will extend beyond the back face of the plate.')
        validFlag = false;
        validityComment = 'Crack tip element zone of Crack B will extend beyond the back face of the plate.';
    end
end

%% Define element sizes on the edges on x-y plane close to cracks
%These define the minimum of a linear gradient up to 0.2*b, so the min
%ensures that they are never larger than 0.2*b (which would cause Abaqus CAE
%to error). Also elSize2 should not be too small relative to the size of
%the crack tip zones.

%Side edges
elSizeASide = min([rpA2/3,0.19*b]);       %Edge to the side of the cracks on the 'A' side (minimum of a gradient up to 0.2*b)
if singleCrackFlag                        %Edge to the side of the cracks on the 'B' side - note must be modified in the case of a single crack model
    elSizeBSide = elSizeASide;
else
    if subsurfaceBFlag
        %Surface defect (A) and a subsurface one (B). The elements along the surface where the Crack B approaches must be sufficiently refined.
        elSizeBSide = min([rpA2/3,0.19*b]);
        %If near side of B is closer to the surface than the apex of A, smaller element size needed:
        if (yB-aB1)<=aA1
            crackTipElemEdgeLength = calc_ellipse_circ(aB1+rpB2,aB2+rpB2,100)/(noTipElemsB*2);  %Length of one of the elements in the outer crack tip zone of Crack B
            elSizeBSide = min([elSizeBSide, 0.1*b, crackTipElemEdgeLength/2, (yB-aB1-rpB2)/4]);
        end
    else
        %Two surface defects.
        elSizeBSide = min([rpB2/3,0.19*b]);
    end
end

%Edge opposite the cracks on the other side of the plate (minimum of a gradient up to 0.2*b). nB. Not used when cracks are on opposite edges.                         
if singleCrackFlag
    elSize3 = min([0.1*b,(b-yA-aA1-rpA2)/2]);
else
    elSize3 = min([0.1*b,(b-yA-aA1-rpA2)/2,(b-yB-aB1-rpB2)/2]);
end

%Special case of a single subsurface defect. The elements along the surface
%where the crack approaches must be sufficiently refined.
if singleCrackFlag && subsurfaceAFlag
    elSize4 = min([0.1*b, (yA-aA1-rpA2)/2]);    %Minimum of 0.1b and half distance from crack tip zone edge to plate face
    
    crackTipElemEdgeLength = calc_ellipse_circ(aA1+rpA2,aA2+rpA2,100)/(noTipElemsA*2);  %Length of one of the elements in the outer crack tip zone
    elSize4 = max([elSize4,crackTipElemEdgeLength/2]);
else
    elSize4 = NaN;
end

%Along the x-direction edge at the mouth of each crack. When the ratio of
%aA2 to rpA2 is too low, or aA2<aA1, set as NaN and no seeds will be generated.
if aA2/8 <= elSizeASide
    elSizeAMid1 = NaN;
    elSizeAMid2 = NaN;
else
    elSizeAMid1 = aA2/8;
    elSizeAMid2 = elSizeASide;
end

if subsurfaceBFlag
    elSizeBMid1 = NaN;
    elSizeBMid2 = NaN;
else
    if aB2/8 <= elSizeBSide
        elSizeBMid1 = NaN;
        elSizeBMid2 = NaN;
    else
        elSizeBMid1 = aB2/8;
        elSizeBMid2 = elSizeBSide;
    end
end

%Edge in-between the cracks - this is only used when cracks are on the same
%edge. The central edge is only divided when there is a certain distance
%in-between the cracks.
if oppositeFlag
    %Note that if subsurfaceBFlag = true, oppositeFlag = true as well. So
    %this condition will apply.
    divideFlag = false;
    elSizeMid1 = NaN;
    elSizeMid2 = NaN;
else
    if (d-rpA2-rpB2) > mean([aA2,aB2])  %Not used if the length to seed is less than the average ellipse axis along x
        divideFlag = true;
        elSizeMidDivA1 = (d-rpA2-rpB2)/8;
        elSizeMidDivA2 = elSizeASide;
        elSizeMidDivB1 = (d-rpA2-rpB2)/8;
        elSizeMidDivB2 = elSizeBSide;
    else
        divideFlag = false;
        elSizeMid1 = NaN;
        elSizeMid2 = NaN;
    end
end

%% Define mesh transition plane locations in x direction
meshTransDistA = max([1.5*aA2, 1.5*b]);   %Distance from the edge of the crack to the mesh transition
if singleCrackFlag
    meshTransDistB = max([1.5*aA2, 1.5*b]);
else
    meshTransDistB = max([1.5*aB2, 1.5*b]);
end

meshTransA = xA + aA2 + meshTransDistA;  %On the side with Crack A
if singleCrackFlag
    meshTransB = xA - aA2 - meshTransDistA;  %Other side of Crack A for a single-crack model
else
    meshTransB = xB - aB2 - meshTransDistB;  %On the side with Crack B
end

%If the mesh transition locations are outside of the plate, set them so
%they're just inside. This can occur for narrow plates (or short pipes).
if meshTransA > 0.9*plHalfWidth || meshTransB < -0.9*plHalfWidth
    meshTransA = 0.9*plHalfWidth;
    meshTransB = -0.9*plHalfWidth;
end

meshTransSketchA = -meshTransA; %Due to the way the sketch was drawn, these need to be flipped
meshTransSketchB = -meshTransB;

%Checks on the mesh transition plane locations - note no attampt is made to
%adjust either parameter in case of a warning.
if meshTransDistA <= rpA2
    warning('The crack tip element zone of Crack A will extend beyond the mesh transition plane.')
    validFlag = false;
    validityComment = 'Crack tip element zone of Crack A extend beyond the mesh transition plane.';
end
if meshTransDistB <= rpB2
    warning('The crack tip element zone of Crack B will extend beyond the mesh transition plane.')
    validFlag = false;
    validityComment = 'Crack tip element zone of Crack B extend beyond the mesh transition plane.';
end

%Z direction mesh transition plane location
meshTransZ = b*3;   %This is statically defined to be three times the plate thickness

%% Define crack plane normal vectors
%All allowable combinations of crack positions. Crack plane normal vectors for each crack.
if singleCrackFlag
    %Regardless of subsurfaceAFlag:
    CPNVA = [0,0,-1];
    CPNVB = NaN;
else
    if oppositeFlag
        if subsurfaceBFlag
            if isfield(naturalParamStruct.outputRequests,'contourTypeA')
                %if strcmpi(naturalParamStruct.outputRequests.contourTypeA,'J_INTEGRAL')
                %    CPNVA = [0,0,1];    %One surface, one subsurface. Different direction for J-integral output.
                %    CPNVB = [0,0,1];
                %else
                    CPNVA = [0,0,-1];   %One surface, one subsurface
                    CPNVB = [0,0,1];
                %end
            else
                CPNVA = [0,0,-1];   %One surface, one subsurface
                CPNVB = [0,0,1];
            end
        else
            CPNVA = [0,0,-1];   %Surface cracks on opposite sides
            CPNVB = [0,0,-1];
        end
    else
        CPNVA = [0,0,-1];    %Surface cracks on the same side
        CPNVB = [0,0,-1];
    end
end

%% Construct modelParamStruct - include loading parameters, material parameters and output requests
%Include loading parameters
modelParamStruct.loadingParams = naturalParamStruct.loadingParams;  %Assign loadingParams in modelParamStruct
if isfield(naturalParamStruct,'loadingParams')
    %Biaxiality
    if strcmpi(geometryType,'pipe')
        if biaxialFlag
            if isnan(biaxialStress)
                %Determine pipe axial stress
                disp('biaxialStress given as NaN. Using the axial stress for a closed-ended cylinder.');
                pipeEndForce = pi.*(ri.^2).*pressureMagnitude;                         %Total axial force for a complete pipe
                modelParamStruct.loadingParams.biaxialStress = pipeEndForce./((pi.*((ri+b).^2)) - (pi.*(ri.^2)));     %Pipe axial stress
            end
        end
    else
        if biaxialFlag
            if isnan(biaxialStress)
                error('Biaxial stress is undefined (loadingParams.biaxialStress=NaN and geometryType=''plate'').');
            end
        end
    end
    
    %If necessary, convert per unit width load/moment to overall load/moment
    if isfield(naturalParamStruct.loadingParams,'loadMagnitudeUnitWidthFlag')
        if naturalParamStruct.loadingParams.loadMagnitudeUnitWidthFlag
            modelParamStruct.loadingParams.loadMagnitude = modelParamStruct.loadingParams.loadMagnitude.*2.*plHalfWidth;    %Convert per unit width load/moment to overall
        end
        modelParamStruct.loadingParams = rmfield(modelParamStruct.loadingParams,'loadMagnitudeUnitWidthFlag');
    end
else
    error('No loading parameters specified in naturalParamStruct.')
end

%Include material parameters
if isfield(naturalParamStruct,'materialParams')
    modelParamStruct.materialParams = naturalParamStruct.materialParams;
else
    error('Material parameters not defined in naturalParamStruct.');
end

%Include output requests
if isfield(naturalParamStruct,'outputRequests')
    modelParamStruct.outputRequests = naturalParamStruct.outputRequests;
else
    error('Output requests not defined in naturalParamStruct.');
end

%% Construct modelParamStruct - build the rest of the structure
modelParamStruct.validFlag = validFlag;
modelParamStruct.crackTipElemSizeWarning = crackTipElemSizeWarning;
modelParamStruct.modifyNaturalParamStructStr = modifyNaturalParamStructStr;

modelParamStruct.geometryType = geometryType;
modelParamStruct.pipeCrackA = pipeCrackA;

modelParamStruct.nlgeomFlag = nlgeomFlag;
modelParamStruct.singleCrackFlag = singleCrackFlag;

modelParamStruct.plateSizes.b = b;
modelParamStruct.plateSizes.plDepth = plDepth;
modelParamStruct.plateSizes.ri = ri;
modelParamStruct.plateSizes.plHalfWidth = plHalfWidth;
modelParamStruct.plateSizes.boundBoxMag = boundBoxMag;

modelParamStruct.crackSizes.aA1 = aA1;
modelParamStruct.crackSizes.aA2 = aA2;
modelParamStruct.crackSizes.aB1 = aB1;
modelParamStruct.crackSizes.aB2 = aB2;

modelParamStruct.crackPositions.xA = xA;
modelParamStruct.crackPositions.xB = xB;
modelParamStruct.crackPositions.yA = yA;
modelParamStruct.crackPositions.yB = yB;
modelParamStruct.crackPositions.oppositeFlag = oppositeFlag;
modelParamStruct.crackPositions.subsurfaceAFlag = subsurfaceAFlag;
modelParamStruct.crackPositions.subsurfaceBFlag = subsurfaceBFlag;

modelParamStruct.crackTipZoneSizes.rpA1 = rpA1;
modelParamStruct.crackTipZoneSizes.rpA2 = rpA2;
modelParamStruct.crackTipZoneSizes.rpB1 = rpB1;
modelParamStruct.crackTipZoneSizes.rpB2 = rpB2;

modelParamStruct.elementSizes.elSizeASide = elSizeASide;
modelParamStruct.elementSizes.elSizeBSide = elSizeBSide;
modelParamStruct.elementSizes.elSize3 = elSize3;
modelParamStruct.elementSizes.elSize4 = elSize4;
modelParamStruct.elementSizes.elSizeAMid1 = elSizeAMid1;
modelParamStruct.elementSizes.elSizeAMid2 = elSizeAMid2;
modelParamStruct.elementSizes.elSizeBMid1 = elSizeBMid1;
modelParamStruct.elementSizes.elSizeBMid2 = elSizeBMid2;
if divideFlag
    modelParamStruct.elementSizes.elSizeMidDivA1 = elSizeMidDivA1;
    modelParamStruct.elementSizes.elSizeMidDivA2 = elSizeMidDivA2;
    modelParamStruct.elementSizes.elSizeMidDivB1 = elSizeMidDivB1;
    modelParamStruct.elementSizes.elSizeMidDivB2 = elSizeMidDivB2;
else
    modelParamStruct.elementSizes.elSizeMid1 = elSizeMid1;
    modelParamStruct.elementSizes.elSizeMid2 = elSizeMid2;
end

modelParamStruct.meshTransitionPositions.meshTransA = meshTransA;
modelParamStruct.meshTransitionPositions.meshTransB = meshTransB;
modelParamStruct.meshTransitionPositions.meshTransSketchA = meshTransSketchA;
modelParamStruct.meshTransitionPositions.meshTransSketchB = meshTransSketchB;
modelParamStruct.meshTransitionPositions.divideFlag = divideFlag;
modelParamStruct.meshTransitionPositions.meshTransZ = meshTransZ;

modelParamStruct.meshDensities.noTipElemsA = noTipElemsA;
modelParamStruct.meshDensities.noTipElemsB = noTipElemsB;
modelParamStruct.meshDensities.noTipElemsRadial2A = noTipElemsRadial2A;
modelParamStruct.meshDensities.noTipElemsRadial2B = noTipElemsRadial2B;

modelParamStruct.crackPlaneNormalVectors.CPNVA = CPNVA;
modelParamStruct.crackPlaneNormalVectors.CPNVB = CPNVB;

modelParamStruct.stepStr = stepStr;

%% Calculate edge parameters
modelParamStruct.edgeCoords = int_defects_calc_model_edge_params(modelParamStruct);

%% Provide optional output argument (validityComment)
if nargout == 2
    varargout{1} = validityComment;
elseif nargout > 2
    error('Too many output arguments.')
end

end