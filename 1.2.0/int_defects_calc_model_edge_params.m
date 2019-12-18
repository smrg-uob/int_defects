function [ edgeCoords ] = int_defects_calc_model_edge_params( modelParamStructIn )
%int_defects_calc_model_edge_params.m
%Harry Coules 2015
%
%DESCRIPTION
%This function determines the coordinates of edges, surfaces and cells in a
%model of interacting semi-elliptical cracks in a plate. It is called by
%int_defects_calc_model_edge_params. The coordinates calculated by this
%function are used in automatic generation of .py files, for selecting
%parts of the geometry.
%
%INPUT ARGUMENTS
% modelParamStructIn - Structure containing model geometry information.
%
%OUTPUT ARGUMENTS
% edgeCoords - Structure containing coordinate inromation used for
%   selecting model edges.

%% Unpack structure
geometryType = modelParamStructIn.geometryType;
if isfield(modelParamStructIn,'pipeCrackA')
    pipeCrackA = modelParamStructIn.pipeCrackA;
end
singleCrackFlag = modelParamStructIn.singleCrackFlag;

b = modelParamStructIn.plateSizes.b;
plDepth = modelParamStructIn.plateSizes.plDepth;
plHalfWidth = modelParamStructIn.plateSizes.plHalfWidth;

aA1 = modelParamStructIn.crackSizes.aA1;
aA2 = modelParamStructIn.crackSizes.aA2;
aB1 = modelParamStructIn.crackSizes.aB1;
aB2 = modelParamStructIn.crackSizes.aB2;

xA = modelParamStructIn.crackPositions.xA;
xB = modelParamStructIn.crackPositions.xB;
yA = modelParamStructIn.crackPositions.yA;
yB = modelParamStructIn.crackPositions.yB;
oppositeFlag = modelParamStructIn.crackPositions.oppositeFlag;
subsurfaceAFlag = modelParamStructIn.crackPositions.subsurfaceAFlag;
subsurfaceBFlag = modelParamStructIn.crackPositions.subsurfaceBFlag;

rpA1 = modelParamStructIn.crackTipZoneSizes.rpA1;
rpA2 = modelParamStructIn.crackTipZoneSizes.rpA2;
rpB1 = modelParamStructIn.crackTipZoneSizes.rpB1;
rpB2 = modelParamStructIn.crackTipZoneSizes.rpB2;

meshTransA = modelParamStructIn.meshTransitionPositions.meshTransA;
meshTransB = modelParamStructIn.meshTransitionPositions.meshTransB;
meshTransSketchA = modelParamStructIn.meshTransitionPositions.meshTransSketchA;
meshTransSketchB = modelParamStructIn.meshTransitionPositions.meshTransSketchB;
divideFlag = modelParamStructIn.meshTransitionPositions.divideFlag;
meshTransZ = modelParamStructIn.meshTransitionPositions.meshTransZ;

smallDist = 2e-6; %An arbitrarily-small distance. Note this is defined statically, due to the need for a specific precision which arises from CAE's behaviour.
                  %See:  Abaqus scripting reference manual, Section 7.4.1: findAt(...)
%% Crack A
%Radial lines in fine crack tip zone
edgeCoords.A.radial1 = [xA+aA2+(0.5*rpA1),yA,0;
    xA+aA2-(0.5*rpA1),yA,0;
    xA+aA2,yA,0.5*rpA1;
    xA-aA2+(0.5*rpA1),yA,0;
    xA-aA2-(0.5*rpA1),yA,0;
    xA-aA2,yA,0.5*rpA1];
%Radial lines in coarse crack tip zone
edgeCoords.A.radial2 = [xA+aA2+(0.5*rpA2),yA,0;
    xA+aA2-(0.5*rpA2),yA,0;
    xA+aA2,yA,0.5*rpA2;
    xA-aA2+(0.5*rpA2),yA,0;
    xA-aA2-(0.5*rpA2),yA,0;
    xA-aA2,yA,0.5*rpA2];
%Circumferential lines
edgeCoords.A.circ = [xA+aA2+rpA1,yA,smallDist;
    xA+aA2-rpA1,yA,smallDist;
    xA+aA2+rpA2,yA,smallDist;
    xA+aA2-rpA2,yA,smallDist;
    xA-aA2+rpA1,yA,smallDist;
    xA-aA2-rpA1,yA,smallDist;
    xA-aA2+rpA2,yA,smallDist;
    xA-aA2-rpA2,yA,smallDist];
%Swept
%nB: Five lines on xy plane, then two below surface.
edgeCoords.A.sweep = [xA+aA2+rpA2,yA+smallDist,0;
    xA+aA2+rpA1,yA+smallDist,0;
    xA+aA2,yA+smallDist,0;
    xA+aA2-rpA1,yA+smallDist,0;
    xA+aA2-rpA2,yA+smallDist,0;
    xA+aA2,yA+smallDist,rpA1;
    xA+aA2,yA+smallDist,rpA2];
%Crack tip
edgeCoords.A.tip = [xA+aA2,yA+smallDist,0];

if subsurfaceAFlag
    %Add extra swept lines (for other half of the ellipse)
    edgeCoords.A.sweep = [edgeCoords.A.sweep;
        xA+aA2+rpA2,yA-smallDist,0;
        xA+aA2+rpA1,yA-smallDist,0;
        xA+aA2,yA-smallDist,0;
        xA+aA2-rpA1,yA-smallDist,0;
        xA+aA2-rpA2,yA-smallDist,0;
        xA+aA2,yA-smallDist,rpA1;
        xA+aA2,yA-smallDist,rpA2];
    %Add extra crack tip line (for other half of the ellipse)
    edgeCoords.A.tip = [edgeCoords.A.tip;
        xA+aA2,yA-smallDist,0];
else
    %Edge in the middle of the crack
    edgeCoords.A.mid = [xA,0,0];
    
    %Coordinates for auto-trim of Crack A semi-ellipse
    edgeCoords.A.autoTrim = [0,-aA1];
end

%% Crack B
if oppositeFlag
        %Radial lines in fine crack tip zone
    edgeCoords.B.radial1 = [xB+aB2+(0.5*rpB1),yB,0;
        xB+aB2-(0.5*rpB1),yB,0;
        xB+aB2,yB,0.5*rpB1;
        xB-aB2+(0.5*rpB1),yB,0;
        xB-aB2-(0.5*rpB1),yB,0;
        xB-aB2,yB,0.5*rpB1];
        %Radial lines in coarse crack tip zone
    edgeCoords.B.radial2 = [xB+aB2+(0.5*rpB2),yB,0;
        xB+aB2-(0.5*rpB2),yB,0;
        xB+aB2,yB,0.5*rpB2;
        xB-aB2+(0.5*rpB2),yB,0;
        xB-aB2-(0.5*rpB2),yB,0;
        xB-aB2,yB,0.5*rpB2];
        %Circumferential lines
    edgeCoords.B.circ = [xB+aB2+rpB1,yB,smallDist;
        xB+aB2-rpB1,yB,smallDist;
        xB+aB2+rpB2,yB,smallDist;
        xB+aB2-rpB2,yB,smallDist;
        xB-aB2+rpB1,yB,smallDist;
        xB-aB2-rpB1,yB,smallDist;
        xB-aB2+rpB2,yB,smallDist;
        xB-aB2-rpB2,yB,smallDist];
        %Swept
    %nB: Five lines on xy plane, then two below surface.
    edgeCoords.B.sweep = [xB+aB2+rpB2,yB-smallDist,0;
        xB+aB2+rpB1,yB-smallDist,0;
        xB+aB2,yB-smallDist,0;
        xB+aB2-rpB1,yB-smallDist,0;
        xB+aB2-rpB2,yB-smallDist,0;
        xB+aB2,yB-smallDist,rpB1;
        xB+aB2,yB-smallDist,rpB2];
    %Crack tip
    edgeCoords.B.tip = [xB+aB2,yB-smallDist,0];
    
    if subsurfaceBFlag
        %Add extra swept lines (for other half of the ellipse)
        edgeCoords.B.sweep = [edgeCoords.B.sweep;
        xB+aB2+rpB2,yB+smallDist,0;
        xB+aB2+rpB1,yB+smallDist,0;
        xB+aB2,yB+smallDist,0;
        xB+aB2-rpB1,yB+smallDist,0;
        xB+aB2-rpB2,yB+smallDist,0;
        xB+aB2,yB+smallDist,rpB1;
        xB+aB2,yB+smallDist,rpB2];
        %Add extra crack tip line (for other half of the ellipse)
        edgeCoords.B.tip = [edgeCoords.B.tip;
            xB+aB2,yB+smallDist,0];
    else
        %Edge in the middle of the crack
        edgeCoords.B.mid = [xB,yB,0];
        
        %Coordinates for auto-trim of Crack B semi-ellipse
        edgeCoords.B.autoTrim = [0,-aB1];
    end
else
    %Radial lines in fine crack tip zone
    edgeCoords.B.radial1 = [xB+aB2+(0.5*rpB1),0,0;
        xB+aB2-(0.5*rpB1),0,0;
        xB+aB2,0,0.5*rpB1;
        xB-aB2+(0.5*rpB1),0,0;
        xB-aB2-(0.5*rpB1),0,0;
        xB-aB2,0,0.5*rpB1];
    %Radial lines in coarse crack tip zone
    edgeCoords.B.radial2 = [xB+aB2+(0.5*rpB2),0,0;
        xB+aB2-(0.5*rpB2),0,0;
        xB+aB2,0,0.5*rpB2;
        xB-aB2+(0.5*rpB2),0,0;
        xB-aB2-(0.5*rpB2),0,0;
        xB-aB2,0,0.5*rpB2];
    %Circumferential lines
    edgeCoords.B.circ = [xB+aB2+rpB1,0,smallDist;
        xB+aB2-rpB1,0,smallDist;
        xB+aB2+rpB2,0,smallDist;
        xB+aB2-rpB2,0,smallDist;
        xB-aB2+rpB1,0,smallDist;
        xB-aB2-rpB1,0,smallDist;
        xB-aB2+rpB2,0,smallDist;
        xB-aB2-rpB2,0,smallDist];
    %Swept
    %nB: Five lines on xy plane, then two below surface.
    edgeCoords.B.sweep = [xB+aB2+rpB2,smallDist,0;
        xB+aB2+rpB1,smallDist,0;
        xB+aB2,smallDist,0;
        xB+aB2-rpB1,smallDist,0;
        xB+aB2-rpB2,smallDist,0;
        xB+aB2,smallDist,rpB1;
        xB+aB2,smallDist,rpB2];
    %Edge in the middle of the crack
    edgeCoords.B.mid = [xB,0,0];
    %Coordinates for auto-trim of Crack B semi-ellipse
    edgeCoords.B.autoTrim = [xB,-yB-aB1];
    %Crack tip
    edgeCoords.B.tip = [xB+aB2,smallDist,0];
end

%% Other geometry
%Around the cracks
if oppositeFlag
    edgeCoords.aroundCracks.sideCracksA.end1 = [meshTransB+rpA2,0,0];
    edgeCoords.aroundCracks.sideCracksA.end2 = [meshTransA-rpA2,0,0];
    edgeCoords.aroundCracks.sideCracksB.end1 = [meshTransA-rpB2,b,0];
    edgeCoords.aroundCracks.sideCracksB.end2 = [meshTransB+rpB2,b,0];
else
    edgeCoords.aroundCracks.sideCracksA.end2 = [meshTransA-rpA2,0,0];   %Lines either side of the cracks on the x axis
    if singleCrackFlag
        edgeCoords.aroundCracks.sideCracksB.end1 = [meshTransB+rpA2,0,0];   %Need this condition because rpB2=NaN for a single crack, so you can't use it here.
    else
        edgeCoords.aroundCracks.sideCracksB.end1 = [meshTransB+rpB2,0,0];
    end
end

if divideFlag
    edgeCoords.aroundCracks.inBetweenCracks.A.end1 = [rpA2/100,0,0];    %Line along x axis in between cracks, A side
    edgeCoords.aroundCracks.inBetweenCracks.B.end2 = [-rpB2/100,0,0];   %Line along x axis in between cracks, B side
    edgeCoords.aroundCracks.oppositeCracks.A.end1 = [rpA2/100,b,0];     %Line on opposite side of plate, A side
    edgeCoords.aroundCracks.oppositeCracks.B.end2 = [-rpB2/100,b,0];    %Line on opposite side of plate, B side
else
    edgeCoords.aroundCracks.inBetweenCracks = [0,0,0];                  %Line along x axis in between cracks
    edgeCoords.aroundCracks.oppositeCracks = [0,b,0];                   %Line on opposite side of plate, opposite cracks
    edgeCoords.aroundCracks.sameSideAsCracks = [0,0,0];                 %Line on y=0 side of plate, used when subsurfaceAFlag = true.
end

%Other edges aligned in the x-direction
if divideFlag
    edgeCoords.otherXdir.crackZone = [rpA2/100,0,meshTransZ;    %Edges along x direction in crack zone at meshTransZ
        -rpB2/100,0,meshTransZ;
        rpA2/100,b,meshTransZ;
        -rpB2/100,b,meshTransZ];
else
    edgeCoords.otherXdir.crackZone = [0,0,meshTransZ;
    0,b,meshTransZ];
end

edgeCoords.otherXdir.outer = [-0.99*plHalfWidth,0,0; %Edges along x direction at z=0 and z=meshTransZ away from crack zone
    -0.99*plHalfWidth,b,0;
    0.99*plHalfWidth,0,0;
    0.99*plHalfWidth,b,0;
    -0.99*plHalfWidth,0,meshTransZ;
    -0.99*plHalfWidth,b,meshTransZ;
    0.99*plHalfWidth,0,meshTransZ;
    0.99*plHalfWidth,b,meshTransZ];

edgeCoords.otherXdir.lower = [0,0,plDepth; %Edges along x direction at the other side of the plate
    0,b,plDepth];

%All edges aligned in the y-direction
edgeCoords.otherYdir.all = [-plHalfWidth,0.5*b,0; %All edges aligned with the y direction
    -plHalfWidth,0.5*b,meshTransZ;
    plHalfWidth,0.5*b,0;
    plHalfWidth,0.5*b,meshTransZ;
    meshTransA,0.5*b,0;
    meshTransA,0.5*b,meshTransZ;
    meshTransB,0.5*b,0;
    meshTransB,0.5*b,meshTransZ;
    -plHalfWidth,0.5*b,plDepth;
    plHalfWidth,0.5*b,plDepth];

%Other edges aligned in the z-direction
edgeCoords.otherZdir.crackZone.end1 = [-plHalfWidth,0,0.1*meshTransZ;    %Edges along the z direction at the same height in z as the crack zone
    -plHalfWidth,b,0.1*meshTransZ;
    plHalfWidth,0,0.1*meshTransZ;
    plHalfWidth,b,0.1*meshTransZ;
    meshTransA,b,0.1*meshTransZ;
    meshTransB,b,0.1*meshTransZ];
edgeCoords.otherZdir.crackZone.end2 = [ meshTransA,0,0.1*meshTransZ;
    meshTransB,0,0.1*meshTransZ];
edgeCoords.otherZdir.lower = [-plHalfWidth,0,0.5*plDepth; %Edges along the z direction lower down in z
    -plHalfWidth,b,0.5*plDepth;
    plHalfWidth,0,0.5*plDepth;
    plHalfWidth,b,0.5*plDepth];

%% Coordinates used for defining sets
%Set for zSymm boundary condition
if singleCrackFlag
    if oppositeFlag
        error('Cannot have singleCrackFlag=true and oppositeFlag=true.')
    elseif divideFlag
        error('Cannot have singleCrackFlag=true and divideFlag=true.')
    elseif subsurfaceBFlag
        error('Cannot have singleCrackFlag=true and subsurfaceBFlag=true.')
    else
        edgeCoords.sets.zSymm = [meshTransA-rpA2,0.5*b,0; %Select central region
            0.99*plHalfWidth,0.5*b,0;                 %Faces either side of central region
            -0.99*plHalfWidth,0.5*b,0;
            xA+aA2+(rpA1/2),yA+(rpA1/10),0;  %Faces of crack tip zones
            xA+aA2+(rpA2/2),yA+(rpA2/10),0];
        if subsurfaceAFlag
            %Extra faces of the other side of the subsurface Crack A
            %ellipse.
            edgeCoords.sets.zSymm = [edgeCoords.sets.zSymm;
                xA+aA2+(rpA1/2),yA-(rpA1/10),0;
                xA+aA2+(rpA2/2),yA-(rpA2/10),0];
        end
    end
else
    if oppositeFlag
        if divideFlag
            error('Cannot have oppositeFlag=true and divideFlag=true.')
        else
            edgeCoords.sets.zSymm = [meshTransA-rpA2,0.5*b,0;
                0.99*plHalfWidth,0.5*b,0;                 %Faces either side of central region
                -0.99*plHalfWidth,0.5*b,0;
                xA+aA2+(rpA1/2),rpA1/10,0;  %Faces of crack tip zones
                xA+aA2+(rpA2/2),rpA2/10,0;
                xB-aB2-(rpB1/2),yB-(rpB1/10),0;
                xB-aB2-(rpB2/2),yB-(rpB2/10),0];
            if subsurfaceBFlag
                %Extra faces of the other side of the subsurface Crack B
                %ellipse.
                edgeCoords.sets.zSymm = [edgeCoords.sets.zSymm;
                    xB-aB2-(rpB1/2),yB+(rpB1/10),0;
                    xB-aB2-(rpB2/2),yB+(rpB2/10),0];
            end
        end
    else
        edgeCoords.sets.zSymm = [meshTransA-rpA2,0.5*b,0;
            0.99*plHalfWidth,0.5*b,0;                 %Faces either side of central region
            -0.99*plHalfWidth,0.5*b,0;
            xA+aA2+(rpA1/2),rpA1/10,0;  %Faces of crack tip zones
            xA+aA2+(rpA2/2),rpA2/10,0;
            xB-aB2-(rpB1/2),rpB1/10,0;
            xB-aB2-(rpB2/2),rpB2/10,0];
        if divideFlag
            edgeCoords.sets.zSymm = [meshTransB+rpB2,0.5*b,0; %Crack region is divided, add right-hand side
                edgeCoords.sets.zSymm];
        end
    end
end

if strcmpi(geometryType,'pipe')             %For a pipe, the opposite face also need to have zSymm applied
    edgeCoords.sets.zSymm = [edgeCoords.sets.zSymm;
        0,0.5*b,plDepth];
end

%Set for crack face pressure boundary condition
if singleCrackFlag
    edgeCoords.sets.crackFacePressureA = [xA,yA+(rpA2/10),0;   %Main surfaces of Crack A
        xA+aA2-(rpA1/2),yA+(rpA1/10),0;                       %Faces of crack tip zone for Crack A
        xA+aA2-(rpA2/2),yA+(rpA2/10),0];
    if subsurfaceAFlag
        %Add the extra faces of the crack tip zones which occur when
        %Crack A is a subsurface crack.
        edgeCoords.sets.crackFacePressureA = [edgeCoords.sets.crackFacePressureA;
            xA+aA2-(rpA1/2),yA-(rpA1/10),0;
            xA+aA2-(rpA2/2),yA-(rpA2/10),0];
    end
    edgeCoords.sets.crackFacePressure = edgeCoords.sets.crackFacePressureA;
else
    if oppositeFlag
        edgeCoords.sets.crackFacePressureA = [xA,rpA2/10,0;   %Crack A
            xA+aA2-(rpA1/2),rpA1/10,0;                        %Faces of crack tip zone
            xA+aA2-(rpA2/2),rpA2/10,0];
        edgeCoords.sets.crackFacePressureB = [xB,yB-(rpB2/10),0;    %Crack B
            xB-aB2+(rpB1/2),yB-(rpB1/10),0;                         %Faces of crack tip zone
            xB-aB2+(rpB2/2),yB-(rpB2/10),0];
        if subsurfaceBFlag
            %Add the extra faces of the crack tip zones which occur when
            %Crack B is a subsurface crack.
            edgeCoords.sets.crackFacePressureB = [edgeCoords.sets.crackFacePressureB;
                xB-aB2+(rpB1/2),yB+(rpB1/10),0;
                xB-aB2+(rpB2/2),yB+(rpB2/10),0];
        end
    else
        edgeCoords.sets.crackFacePressureA = [xA,rpA2/10,0;   %Crack A
            xA+aA2-(rpA1/2),rpA1/10,0;                        %Faces of crack tip zone
            xA+aA2-(rpA2/2),rpA2/10,0];
        edgeCoords.sets.crackFacePressureB = [xB,rpB2/10,0;   %Crack B
            xB-aB2+(rpB1/2),rpB1/10,0;                        %Faces of crack tip zone
            xB-aB2+(rpB2/2),rpB2/10,0];
    end
    edgeCoords.sets.crackFacePressure = [edgeCoords.sets.crackFacePressureA; edgeCoords.sets.crackFacePressureB];
end

%Set for pipe internal pressure boundary condition - excluding crack face
%pressure loading
if strcmpi(geometryType,'pipe')
    if strcmpi(pipeCrackA,'internal')
        if singleCrackFlag
            edgeCoords.sets.pipePressure = [xA,0,0.5*meshTransZ;    %Central remote faces
                xA,0,meshTransZ+rpA2;
                meshTransA+rpA2,0,0.5*meshTransZ;                   %Left and right remote faces
                meshTransB-rpA2,0,0.5*meshTransZ];
            if ~subsurfaceAFlag
                edgeCoords.sets.pipePressure = [edgeCoords.sets.pipePressure;
                    xA-aA2-(rpA1/2),0,rpA1/10;                          %Crack A -ive x side
                    xA-aA2-(rpA2/2),0,rpA1/10;
                    xA-aA2+(rpA1/2),0,rpA1/10;
                    xA-aA2+(rpA2/2),0,rpA1/10;
                    xA+aA2-(rpA1/2),0,rpA1/10;                          %Crack A +ive x side
                    xA+aA2-(rpA2/2),0,rpA1/10;
                    xA+aA2+(rpA1/2),0,rpA1/10;
                    xA+aA2+(rpA2/2),0,rpA1/10];
            end
        else
            edgeCoords.sets.pipePressure = [xA,0,0.5*meshTransZ;    %Central remote faces
                xA,0,meshTransZ+rpA2;
                meshTransA+rpA2,0,0.5*meshTransZ;                   %Left and right remote faces
                meshTransB-rpA2,0,0.5*meshTransZ;
                xA-aA2-(rpA1/2),0,rpA1/10;                          %Crack A -ive x side
                xA-aA2-(rpA2/2),0,rpA1/10;
                xA-aA2+(rpA1/2),0,rpA1/10;
                xA-aA2+(rpA2/2),0,rpA1/10;
                xA+aA2-(rpA1/2),0,rpA1/10;                          %Crack A +ive x side
                xA+aA2-(rpA2/2),0,rpA1/10;
                xA+aA2+(rpA1/2),0,rpA1/10;
                xA+aA2+(rpA2/2),0,rpA1/10];
            if ~oppositeFlag
                edgeCoords.sets.pipePressure = [edgeCoords.sets.pipePressure;
                    xB-aB2-(rpB1/2),0,rpB1/10;                          %Crack B -ive x side
                    xB-aB2-(rpB2/2),0,rpB1/10;
                    xB-aB2+(rpB1/2),0,rpB1/10;
                    xB-aB2+(rpB2/2),0,rpB1/10;
                    xB+aB2-(rpB1/2),0,rpB1/10;                          %Crack B +ive x side
                    xB+aB2-(rpB2/2),0,rpB1/10;
                    xB+aB2+(rpB1/2),0,rpB1/10;
                    xB+aB2+(rpB2/2),0,rpB1/10];
            end
        end
    elseif strcmpi(pipeCrackA,'external')
        edgeCoords.sets.pipePressure = [xA,b,0.5*meshTransZ;    %Central remote faces
            xA,b,meshTransZ+rpA2;
            meshTransA+rpA2,b,0.5*meshTransZ;                   %Left and right remote faces
            meshTransB-rpA2,b,0.5*meshTransZ];
        if ~singleCrackFlag 
            if oppositeFlag
                if ~subsurfaceBFlag
                    edgeCoords.sets.pipePressure = [edgeCoords.sets.pipePressure;
                        xB-aB2-(rpB1/2),b,rpB1/10;                          %Crack B -ive x side
                        xB-aB2-(rpB2/2),b,rpB1/10;
                        xB-aB2+(rpB1/2),b,rpB1/10;
                        xB-aB2+(rpB2/2),b,rpB1/10;
                        xB+aB2-(rpB1/2),b,rpB1/10;                          %Crack B +ive x side
                        xB+aB2-(rpB2/2),b,rpB1/10;
                        xB+aB2+(rpB1/2),b,rpB1/10;
                        xB+aB2+(rpB2/2),b,rpB1/10];
                end
            end
        end
    else
        error('Unexpected string in pipeCrackA. Should be ''internal'' or ''external''.');
    end
else
    edgeCoords.sets.pipePressure = [];
end

%If loadingParams.loadingType is 'pipePressureInclCrack', the set for pipe
%internal pressure boundary condition must also include any crack faces
%which link with the interior.
if strcmpi(modelParamStructIn.loadingParams.loadingType, 'pipePressureInclCrack')
    if strcmpi(pipeCrackA,'internal')
        if ~subsurfaceAFlag
            edgeCoords.sets.pipePressure = [edgeCoords.sets.pipePressure;
                edgeCoords.sets.crackFacePressureA];
        end
        if ~singleCrackFlag && ~oppositeFlag && ~subsurfaceBFlag
            edgeCoords.sets.pipePressure = [edgeCoords.sets.pipePressure;
                edgeCoords.sets.crackFacePressureB];
        end
    elseif strcmpi(pipeCrackA,'external')
        if ~singleCrackFlag && oppositeFlag && ~subsurfaceBFlag
            edgeCoords.sets.pipePressure = [edgeCoords.sets.pipePressure;
                edgeCoords.sets.crackFacePressureB];
        end
    else
        error('Unexpected string in pipeCrackA. Should be ''internal'' or ''external''.');
    end
end

%Other sets for loading and BCs
edgeCoords.sets.fixed = [plHalfWidth,0,0];                              %First point for fixed displacement BCs
edgeCoords.sets.fixed2 = [-plHalfWidth,0,0];                            %Second point for fixed displacement BCs
edgeCoords.sets.coupling = [0,0.5*b,plDepth];                           %Coupling constraint surface
edgeCoords.sets.biaxialLoadSurf = [-plHalfWidth,0.5*b,0.5*meshTransZ;   %Loading surface for biaxial plate load (pipe axial load)
    -plHalfWidth,0.5*b,plDepth-smallDist];
edgeCoords.sets.xSymm = [plHalfWidth,0.5*b,0.5*meshTransZ;              %xSymm BC surface for biaxial loading
    plHalfWidth,0.5*b,plDepth-smallDist];

%% Cells for assigning mesh controls
if oppositeFlag
    edgeCoords.meshControls.tipZones1 = [xA+aA2+(rpA1/2),rpA1/10,rpA1/10;
        xA+aA2-(rpA1/2),rpA1/10,rpA1/10;
        xB-aB2+(rpB1/2),yB-(rpB1/10),rpB1/10;
        xB-aB2-(rpB1/2),yB-(rpB1/10),rpB1/10];
    edgeCoords.meshControls.tipZones2 = [xA+aA2+(rpA2/2),rpA2/10,rpA2/10;
        xA+aA2-(rpA2/2),rpA2/10,rpA2/10;
        xB-aB2+(rpB2/2),yB-(rpB2/10),rpB2/10;
        xB-aB2-(rpB2/2),yB-(rpB2/10),rpB2/10];
    if subsurfaceBFlag
        %Extra set of cells to be selected on Crack B if subsurface
        edgeCoords.meshControls.tipZones1 = [edgeCoords.meshControls.tipZones1;
            xB-aB2+(rpB1/2),yB+(rpB1/10),rpB1/10;
            xB-aB2-(rpB1/2),yB+(rpB1/10),rpB1/10];
        edgeCoords.meshControls.tipZones2 = [edgeCoords.meshControls.tipZones2;
            xB-aB2+(rpB2/2),yB+(rpB2/10),rpB2/10;
            xB-aB2-(rpB2/2),yB+(rpB2/10),rpB2/10];
    end
else
    if singleCrackFlag
        edgeCoords.meshControls.tipZones1 = [xA+aA2+(rpA1/2),yA+(rpA1/10),rpA1/10;
            xA+aA2-(rpA1/2),yA+(rpA1/10),rpA1/10;];
        edgeCoords.meshControls.tipZones2 = [xA+aA2+(rpA2/2),yA+(rpA2/10),rpA2/10;
            xA+aA2-(rpA2/2),yA+(rpA2/10),rpA2/10;];
        if subsurfaceAFlag
            %Extra set of cells to be selected on Crack B if subsurface
            edgeCoords.meshControls.tipZones1 = [edgeCoords.meshControls.tipZones1;
                xA+aA2+(rpA1/2),yA-(rpA1/10),rpA1/10;
                xA+aA2-(rpA1/2),yA-(rpA1/10),rpA1/10;];
            edgeCoords.meshControls.tipZones2 = [edgeCoords.meshControls.tipZones2;
                xA+aA2+(rpA2/2),yA-(rpA2/10),rpA2/10;
                xA+aA2-(rpA2/2),yA-(rpA2/10),rpA2/10;];
        end
    else
        edgeCoords.meshControls.tipZones1 = [xA+aA2+(rpA1/2),rpA1/10,rpA1/10;
            xA+aA2-(rpA1/2),rpA1/10,rpA1/10;
            xB-aB2+(rpB1/2),rpB1/10,rpB1/10;
            xB-aB2-(rpB1/2),rpB1/10,rpB1/10];
        edgeCoords.meshControls.tipZones2 = [xA+aA2+(rpA2/2),rpA2/10,rpA2/10;
            xA+aA2-(rpA2/2),rpA2/10,rpA2/10;
            xB-aB2+(rpB2/2),rpB2/10,rpB2/10;
            xB-aB2-(rpB2/2),rpB2/10,rpB2/10];
    end
end

%% Partitioning of model cells
if divideFlag
    edgeCoords.meshControls.aroundCracks = [meshTransA-rpA2,0.5*b,1;  %Crack region is divided, add right-hand side
        meshTransB+rpB2,0.5*b,1];
else
    edgeCoords.meshControls.aroundCracks = [meshTransA-rpA2,0.5*b,1];
end

if subsurfaceAFlag
    %Used for partitioning Crack A tip cells using a datum plane if it is a
    %subsurface crack.
    edgeCoords.ellipseToDivide.A = [xA-aA2+(rpA1/2),yA,rpA1/10;
        xA-aA2-(rpA1/2),yA,rpA1/10;
        xA-aA2+(rpA2/2),yA,rpA2/10;
        xA-aA2-(rpA2/2),yA,rpA2/10];
end

if subsurfaceBFlag
    %Used for partitioning Crack B tip cells using a datum plane if it is a
    %subsurface crack.
    edgeCoords.ellipseToDivide.B = [xB-aB2+(rpB1/2),yB,rpB1/10;
        xB-aB2-(rpB1/2),yB,rpB1/10;
        xB-aB2+(rpB2/2),yB,rpB2/10;
        xB-aB2-(rpB2/2),yB,rpB2/10];
end

