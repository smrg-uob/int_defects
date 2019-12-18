function [ lineOut ] = int_defects_py_line_edges( lineIn, modelParamStruct)
%int_defects_py_line_edges.m
%Harry Coules 2015
%
%DESCRIPTION
%This function is used to write the lines of a Python script used
%generating an Abaqus model. Specifically, this function is used for the
%later parts of the .py file where it is necessary to reselect existing
%geometry in order to assign mesh seeds etc. For earlier parts of the
%Python script, where this geometry is being created,
%in_defects_py_line_basic is used. Both functions are called by
%int_defects_write_py.
%
%INPUT ARGUMENTS
%lineIn - String containing the current line in the master .py file.
%modelParamStruct - Structure containing parameters for the current model.
%
%OUTPUT ARGUMENTS
%lineOut - A string or cell array of strings which will be used to replace
%   the current line (i.e. lineIn) in the generated .py file.
%%
lineOut = lineIn;
if ~isfield(modelParamStruct,'relaxMeshFlag1')
    modelParamStruct.relaxMeshFlag1 = false;
end
if ~isfield(modelParamStruct,'relaxMeshFlag2')
    modelParamStruct.relaxMeshFlag2 = false;
end

%formatSpec = '%8.8e';   %Specify the number format to be written.
formatSpec = '%.7f';   %Specify the number format to be written.

%% Material
if any(strfind(lineIn,'#material#'))    
    if strcmpi(modelParamStruct.materialParams.type,'linear elasticity')
        materialParamsTable = [modelParamStruct.materialParams.E, modelParamStruct.materialParams.nu];
        line1 = 'mdb.models[''Model-1''].Material(name=''Material-1'')';
        line2 = 'mdb.models[''Model-1''].materials[''Material-1''].Elastic(table=';
        line3 = int_defects_py_line_builder(materialParamsTable, formatSpec);
        lineOut = {line1;line2;line3};
    elseif strcmpi(modelParamStruct.materialParams.type,'deformation plasticity')
        materialParamsTable = [modelParamStruct.materialParams.E,...
            modelParamStruct.materialParams.nu,...
            modelParamStruct.materialParams.s0,...
            modelParamStruct.materialParams.n,...
            modelParamStruct.materialParams.alpha]';
        line1 = 'mdb.models[''Model-1''].Material(name=''Material-1'')';
        line2 = 'mdb.models[''Model-1''].materials[''Material-1''].DeformationPlasticity(table=';
        line3 = int_defects_py_line_builder(materialParamsTable, formatSpec, 'defplast');
        lineOut = {line1;line2;line3};
    elseif strcmpi(modelParamStruct.materialParams.type,'incremental plasticity')
        %Note that elastic properties but be defined also, whenever
        %incremental plasticity is used.
        materialParamsTable = [modelParamStruct.materialParams.E, modelParamStruct.materialParams.nu];
        line1 = 'mdb.models[''Model-1''].Material(name=''Material-1'')';
        line2 = 'mdb.models[''Model-1''].materials[''Material-1''].Elastic(table=';
        line3 = int_defects_py_line_builder(materialParamsTable, formatSpec);
        line4 = 'mdb.models[''Model-1''].materials[''Material-1''].Plastic(table=';
        line5 = int_defects_py_line_builder(modelParamStruct.materialParams.plasticTable, formatSpec, 'incplast');
        lineOut = {line1;line2;line3;line4;line5};
    end
end

%% Contour integrals
%Crack A
if any(strfind(lineIn,'#contourIntegralA#'))
    if ~strcmpi(modelParamStruct.outputRequests.contourTypeA,'none')
        line1 = 'mdb.models[''Model-1''].HistoryOutputRequest(contourIntegral=''CrackA'',';
        line2 = ['contourType=', modelParamStruct.outputRequests.contourTypeA,', createStepName=''ApplyLoad'', frequency=', num2str(modelParamStruct.outputRequests.contourFreqA),];
        line3 = [', name=''', modelParamStruct.outputRequests.contourNameA, ''', numberOfContours=',num2str(modelParamStruct.outputRequests.noContoursA), ', rebar=EXCLUDE, sectionPoints=DEFAULT)'];
        lineOut = {line1;line2;line3};
    end
end

%Crack B
if any(strfind(lineIn,'#contourIntegralB#'))
    if ~strcmpi(modelParamStruct.outputRequests.contourTypeB,'none')
        line1 = 'mdb.models[''Model-1''].HistoryOutputRequest(contourIntegral=''CrackB'',';
        line2 = ['contourType=', modelParamStruct.outputRequests.contourTypeB,', createStepName=''ApplyLoad'', frequency=', num2str(modelParamStruct.outputRequests.contourFreqB),];
        line3 = [', name=''', modelParamStruct.outputRequests.contourNameB, ''', numberOfContours=',num2str(modelParamStruct.outputRequests.noContoursB), ', rebar=EXCLUDE, sectionPoints=DEFAULT)'];
        lineOut = {line1;line2;line3};
    end
end

%% Step parameters
%modelParamStruct.stepStr can be directly written into the py file if it
%exists.
if any(strfind(lineIn,'#step#'))
    if isempty(modelParamStruct.stepStr)
        lineOut = '#NOTE: This step uses default parameters.';
    else
        lineOut = modelParamStruct.stepStr;
    end
end

%% Crack A
%Radial lines in inner crack tip zone
if any(strfind(lineIn,'#edgeCoords.A.radial1#'))
    line1 = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByNumber(constraint=FINER, edges=';
    line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.A.radial1, formatSpec);
    line4 = ', number=1)';
    lineOut = {line1;line2;line3;line4};
end

%Radial lines in outer crack tip zone
if any(strfind(lineIn,'#edgeCoords.A.radial2#'))
    line1 = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByNumber(constraint=FINER, edges=';
    line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.A.radial2, formatSpec);
    line4 = [', number=',num2str(modelParamStruct.meshDensities.noTipElemsRadial2A),')'];
    lineOut = {line1;line2;line3;line4};
end

%Circumferential lines
if any(strfind(lineIn,'#edgeCoords.A.circ#'))
    line1 = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByNumber(constraint=FINER, edges=';
    line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.A.circ, formatSpec);
    line4 = ', number=5)';
    lineOut = {line1;line2;line3;line4};
end

%Swept edges - Note FIXED number of crack tip elements
if any(strfind(lineIn,'#edgeCoords.A.sweep#'))
    line1 = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByNumber(constraint=FIXED, edges=';
    line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.A.sweep, formatSpec);
    line4 = [', number=',num2str(modelParamStruct.meshDensities.noTipElemsA),')'];
    lineOut = {line1;line2;line3;line4};
end

%Line along x axis at crack mouth
if any(strfind(lineIn,'#edgeCoords.A.mid#'))
    if modelParamStruct.crackPositions.subsurfaceAFlag
        lineOut = '#NOTE: No code written for edgeCoords.A.mid. Crack A is a complete ellipse and has no crack mouth to be seeded.';
    elseif isnan(modelParamStruct.elementSizes.elSizeAMid1)
        lineOut = '#NOTE: No code written here. The ratio between aA2 and rpA2 is too small for the crack mouth to be seeded.';
    elseif modelParamStruct.relaxMeshFlag2
        lineOut = '#NOTE: No code written here. modeParamStruct.relaxMeshFlag2 = true. Mesh constraints have been relaxed.';
    elseif modelParamStruct.relaxMeshFlag1
        lineOut = '#NOTE: No code written here. modeParamStruct.relaxMeshFlag1 = true. Mesh constraints have been relaxed.';
    else
        line1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByBias(biasMethod=DOUBLE,';
        line1b = 'endEdges=';
        line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
        line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.A.mid, formatSpec);
        line4 = [', constraint=FINER, maxSize=',num2str(modelParamStruct.elementSizes.elSizeAMid1),', minSize=',num2str(modelParamStruct.elementSizes.elSizeAMid2),')'];
        lineOut = {line1a;line1b;line2;line3;line4};
    end
end

%% Crack B
%Radial lines in inner crack tip zone
if any(strfind(lineIn,'#edgeCoords.B.radial1#'))
    line1 = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByNumber(constraint=FINER, edges=';
    line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.B.radial1, formatSpec);
    line4 = ', number=1)';
    lineOut = {line1;line2;line3;line4};
end

%Radial lines in outer crack tip zone
if any(strfind(lineIn,'#edgeCoords.B.radial2#'))
    line1 = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByNumber(constraint=FINER, edges=';
    line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.B.radial2, formatSpec);
    line4 = [', number=',num2str(modelParamStruct.meshDensities.noTipElemsRadial2B),')'];
    lineOut = {line1;line2;line3;line4};
end

%Circumferential lines
if any(strfind(lineIn,'#edgeCoords.B.circ#'))
    line1 = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByNumber(constraint=FINER, edges=';
    line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.B.circ, formatSpec);
    line4 = ', number=5)';
    lineOut = {line1;line2;line3;line4};
end

%Swept edges - Note FIXED number of crack tip elements
if any(strfind(lineIn,'#edgeCoords.B.sweep#'))
    line1 = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByNumber(constraint=FIXED, edges=';
    line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.B.sweep, formatSpec);
    line4 = [', number=',num2str(modelParamStruct.meshDensities.noTipElemsB),')'];
    lineOut = {line1;line2;line3;line4};
end

%Line along x axis at crack mouth
if any(strfind(lineIn,'#edgeCoords.B.mid#'))
    if modelParamStruct.crackPositions.subsurfaceBFlag
        lineOut = '#NOTE: No code written for edgeCoords.B.mid. Crack B is a complete ellipse and has no crack mouth to be seeded.';
    elseif isnan(modelParamStruct.elementSizes.elSizeBMid1)
        lineOut = '#NOTE: No code written here. The ratio between aB2 and rpB2 is too small for the crack mouth to be seeded.';
    elseif modelParamStruct.relaxMeshFlag2
        lineOut = '#NOTE: No code written here. modeParamStruct.relaxMeshFlag2 = true. Mesh constraints have been relaxed.';
    elseif modelParamStruct.relaxMeshFlag1
        lineOut = '#NOTE: No code written here. modeParamStruct.relaxMeshFlag1 = true. Mesh constraints have been relaxed.';
    else
        line1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByBias(biasMethod=DOUBLE,';
        line1b = 'endEdges=';
        line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
        line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.B.mid, formatSpec);
        line4 = [', constraint=FINER, maxSize=',num2str(modelParamStruct.elementSizes.elSizeBMid1),', minSize=',num2str(modelParamStruct.elementSizes.elSizeBMid2),')'];
        lineOut = {line1a;line1b;line2;line3;line4};
    end
end

%% Other geometry - around the cracks
%Line along x axis in between cracks - note this is not used for cracks on
%opposite sides of the plate, or when there is only a single crack.
if any(strfind(lineIn,'#edgeCoords.aroundCracks.inBetweenCracks#'))
    if modelParamStruct.singleCrackFlag
        lineOut = '#NOTE: No code has been written here - none needed when there is only a single crack.';
    elseif modelParamStruct.crackPositions.oppositeFlag
        lineOut = '#NOTE: No code has been written here - none needed when cracks on opposite sides.';
    elseif modelParamStruct.relaxMeshFlag2
        lineOut = '#NOTE: No code written here. modeParamStruct.relaxMeshFlag2 = true. Mesh constraints have been relaxed.';
    else
        if modelParamStruct.meshTransitionPositions.divideFlag
            lineA1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByBias(biasMethod=SINGLE,';
            lineA1b = 'constraint=FINER, end1Edges=';
            lineA1c = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
            lineA2 = int_defects_py_line_builder(modelParamStruct.edgeCoords.aroundCracks.inBetweenCracks.A.end1, formatSpec);
            lineA3 = [',maxSize=',num2str(modelParamStruct.elementSizes.elSizeMidDivA1),', minSize=',num2str(modelParamStruct.elementSizes.elSizeMidDivA2),')'];
            lineB1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByBias(biasMethod=SINGLE,';
            lineB1b = 'constraint=FINER, end2Edges=';
            lineB1c = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
            lineB2 = int_defects_py_line_builder(modelParamStruct.edgeCoords.aroundCracks.inBetweenCracks.B.end2, formatSpec);
            lineB3 = [',maxSize=',num2str(modelParamStruct.elementSizes.elSizeMidDivB1),', minSize=',num2str(modelParamStruct.elementSizes.elSizeMidDivB2),')'];
            lineOut = {lineA1a;lineA1b;lineA1c;lineA2;lineA3;lineB1a;lineB1b;lineB1c;lineB2;lineB3};
        else
            if ~isnan(modelParamStruct.elementSizes.elSizeMid1)
                line1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByBias(biasMethod=DOUBLE,';
                line1b = 'endEdges=';
                line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
                line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.aroundCracks.inBetweenCracks, formatSpec);
                line4 = [', constraint=FINER, maxSize=',num2str(modelParamStruct.elementSizes.elSizeMid1),', minSize=',num2str(modelParamStruct.elementSizes.elSizeMid2),')'];
                lineOut = {line1a;line1b;line2;line3;line4};
            else
                lineOut = '#NOTE: Seeds for edge in between cracks not written - cracks are too close together.';
            end
        end
    end
end

%Line on other side of plate, opposite cracks - note this is not used for
%cracks on opposite sides of the plate.
if any(strfind(lineIn,'#edgeCoords.aroundCracks.oppositeCracks#'))
    if modelParamStruct.crackPositions.oppositeFlag && ~modelParamStruct.crackPositions.subsurfaceBFlag
        lineOut = '#NOTE: No code has been written here. None needed when surface cracks emanate from opposite sides.';
    else
        if modelParamStruct.meshTransitionPositions.divideFlag
            lineA1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByBias(biasMethod=SINGLE,';
            lineA1b = 'constraint=FINER, end1Edges=';
            lineA1c = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
            lineA2 = int_defects_py_line_builder(modelParamStruct.edgeCoords.aroundCracks.oppositeCracks.A.end1, formatSpec);
            lineA3 = [', maxSize=',num2str(0.2*modelParamStruct.plateSizes.b),', minSize=',num2str(modelParamStruct.elementSizes.elSize3),')'];
            lineB1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByBias(biasMethod=SINGLE,';
            lineB1b = 'constraint=FINER, end2Edges=';
            lineB1c = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
            lineB2 = int_defects_py_line_builder(modelParamStruct.edgeCoords.aroundCracks.oppositeCracks.B.end2, formatSpec);
            lineB3 = [', maxSize=',num2str(0.2*modelParamStruct.plateSizes.b),', minSize=',num2str(modelParamStruct.elementSizes.elSize3),')'];
            lineOut = {lineA1a;lineA1b;lineA1c;lineA2;lineA3;lineB1a;lineB1b;lineB1c;lineB2;lineB3};
        else
            %Note that this will cover the condition when the crack or
            %cracks are only on one side, and also when Crack B is a
            %subsurface crack.
            line1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByBias(biasMethod=DOUBLE,';
            line1b = 'centerEdges=';
            line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
            line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.aroundCracks.oppositeCracks, formatSpec);
            line4 = [', constraint=FINER, maxSize=',num2str(0.2*modelParamStruct.plateSizes.b),', minSize=',num2str(modelParamStruct.elementSizes.elSize3),')'];
            lineOut = {line1a;line1b;line2;line3;line4};
        end
    end
end

%Line on the y=0 side of the plate. It is only used when there is a single
%crack and it is a subsurface defect.
if any(strfind(lineIn,'#edgeCoords.aroundCracks.sameSideAsCracks#'))
    if modelParamStruct.crackPositions.subsurfaceAFlag
        line1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByBias(biasMethod=DOUBLE,';
        line1b = 'centerEdges=';
        line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
        line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.aroundCracks.sameSideAsCracks, formatSpec);
        line4 = [', constraint=FINER, maxSize=',num2str(0.2*modelParamStruct.plateSizes.b),', minSize=',num2str(modelParamStruct.elementSizes.elSize4),')'];    %Note use of elSize4, which is only used in the case of a single subsurface crack.
        lineOut = {line1a;line1b;line2;line3;line4};
    else
        lineOut = '#NOTE: No code has been written here. None needed unless Crack A is a subsurface defect.';
    end
end

%Lines either side of the cracks on the x axis
if any(strfind(lineIn,'#edgeCoords.aroundCracks.eitherSideCracks#'))
    if modelParamStruct.crackPositions.subsurfaceAFlag
        lineOut = '#NOTE: No code written here. Crack A is a subsurface defect, and so the y=0 edge seeds are handled by edgeCoords.aroundCracks.sameSideAsCracks.';
    elseif modelParamStruct.relaxMeshFlag2
        lineOut = '#NOTE: No code written here. modeParamStruct.relaxMeshFlag2 = true. Mesh constraints have been relaxed.';
    elseif modelParamStruct.crackPositions.oppositeFlag
        %Either side of Crack A
        %+ive x side ("A side")
        lineAA1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByBias(biasMethod=SINGLE,';
        lineAA1b = 'constraint=FINER, end2Edges=';
        lineAA1c = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
        lineAA2 = int_defects_py_line_builder(modelParamStruct.edgeCoords.aroundCracks.sideCracksA.end2, formatSpec);
        lineAA3 = [', maxSize=',num2str(0.2*modelParamStruct.plateSizes.b),', minSize=',num2str(modelParamStruct.elementSizes.elSizeASide),')'];
        %-ive x side ("B side")
        lineAB1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByBias(biasMethod=SINGLE,';
        lineAB1b = 'constraint=FINER, end1Edges=';
        lineAB1c = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
        lineAB2 = int_defects_py_line_builder(modelParamStruct.edgeCoords.aroundCracks.sideCracksA.end1, formatSpec);
        if modelParamStruct.crackSizes.aA1>(modelParamStruct.crackPositions.yB - modelParamStruct.crackSizes.aB1)   %Smaller max. element size when Crack B is close to the x-z plane.
            lineAB3 = [', maxSize=',num2str(0.1*modelParamStruct.plateSizes.b),', minSize=',num2str(modelParamStruct.elementSizes.elSizeBSide),')'];
        else
            lineAB3 = [', maxSize=',num2str(0.2*modelParamStruct.plateSizes.b),', minSize=',num2str(modelParamStruct.elementSizes.elSizeBSide),')'];
        end
        
        %If B is a subsurface crack, we only need either side of Crack A
        if modelParamStruct.crackPositions.subsurfaceBFlag
            lineOut = {lineAA1a;lineAA1b;lineAA1c;lineAA2;lineAA3;...
                lineAB1a;lineAB1b;lineAB1c;lineAB2;lineAB3};
        else
            %Either side of Crack B
            %+ive x side
            lineBA1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByBias(biasMethod=SINGLE,';
            lineBA1b = 'constraint=FINER, end2Edges=';
            lineBA1c = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
            lineBA2 = int_defects_py_line_builder(modelParamStruct.edgeCoords.aroundCracks.sideCracksB.end2, formatSpec);
            lineBA3 = [', maxSize=',num2str(0.2*modelParamStruct.plateSizes.b),', minSize=',num2str(modelParamStruct.elementSizes.elSizeBSide),')'];
            %-ive x side
            lineBB1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByBias(biasMethod=SINGLE,';
            lineBB1b = 'constraint=FINER, end1Edges=';
            lineBB1c = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
            lineBB2 = int_defects_py_line_builder(modelParamStruct.edgeCoords.aroundCracks.sideCracksB.end1, formatSpec);
            lineBB3 = [', maxSize=',num2str(0.2*modelParamStruct.plateSizes.b),', minSize=',num2str(modelParamStruct.elementSizes.elSizeBSide),')'];
            lineOut = {lineAA1a;lineAA1b;lineAA1c;lineAA2;lineAA3;...
                lineAB1a;lineAB1b;lineAB1c;lineAB2;lineAB3;
                lineBA1a;lineBA1b;lineBA1c;lineBA2;lineBA3;
                lineBB1a;lineBB1b;lineBB1c;lineBB2;lineBB3};
        end
    else
        %Crack A side
        lineA1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByBias(biasMethod=SINGLE,';
        lineA1b = 'constraint=FINER, end2Edges=';
        lineA1c = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
        lineA2 = int_defects_py_line_builder(modelParamStruct.edgeCoords.aroundCracks.sideCracksA.end2, formatSpec);
        lineA3 = [', maxSize=',num2str(0.2*modelParamStruct.plateSizes.b),', minSize=',num2str(modelParamStruct.elementSizes.elSizeASide),')'];
        %Crack B side
        lineB1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByBias(biasMethod=SINGLE,';
        lineB1b = 'constraint=FINER, end1Edges=';
        lineB1c = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
        lineB2 = int_defects_py_line_builder(modelParamStruct.edgeCoords.aroundCracks.sideCracksB.end1, formatSpec);
        lineB3 = [', maxSize=',num2str(0.2*modelParamStruct.plateSizes.b),', minSize=',num2str(modelParamStruct.elementSizes.elSizeBSide),')'];
        lineOut = {lineA1a;lineA1b;lineA1c;lineA2;lineA3;lineB1a;lineB1b;lineB1c;lineB2;lineB3};
    end
end

%% Other geometry aligned with x-direction
%Edges along x direction in crack zone at z=3
if any(strfind(lineIn,'#edgeCoords.otherXdir.crackZone#'))
    line1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeBySize(constraint=FINER,';
    line1b = 'deviationFactor=0.1, edges=';
    line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.otherXdir.crackZone, formatSpec);
    line4 = [', size=',num2str(2*modelParamStruct.plateSizes.b),')'];
    lineOut = {line1a;line1b;line2;line3;line4};
end

%Edges along x direction at z=0 and z=3 away from crack zone
if any(strfind(lineIn,'#edgeCoords.otherXdir.outer#'))
    line1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeBySize(constraint=FINER,';
    line1b = 'deviationFactor=0.1, edges=';
    line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.otherXdir.outer, formatSpec);
    line4 = [', size=',num2str(modelParamStruct.plateSizes.plHalfWidth/20),')'];
    lineOut = {line1a;line1b;line2;line3;line4};
end

%Edges along x direction at the other side of the plate
if any(strfind(lineIn,'#edgeCoords.otherXdir.lower#'))
    line1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeBySize(constraint=FINER,';
    line1b = 'deviationFactor=0.1, edges=';
    line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.otherXdir.lower, formatSpec);
    line4 = [', size=',num2str(modelParamStruct.plateSizes.plHalfWidth/40),')'];
    lineOut = {line1a;line1b;line2;line3;line4};
end

%% Other geometry aligned with y-direction
if any(strfind(lineIn,'#edgeCoords.otherYdir.all#'))
    line1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeBySize(constraint=FINER,';
    line1b = 'deviationFactor=0.1, edges=';
    line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.otherYdir.all, formatSpec);
    line4 = [', size=',num2str(0.2*modelParamStruct.plateSizes.b),')'];
    lineOut = {line1a;line1b;line2;line3;line4};
end

%% Other geometry aligned with z-direction
%Edges along the z direction at the same height in z as the crack zone
if any(strfind(lineIn,'#edgeCoords.otherZdir.crackZone#'))
    if strcmpi(modelParamStruct.geometryType,'pipe')
        maxNearElemSizeZ = min([modelParamStruct.plateSizes.b,modelParamStruct.plateSizes.plDepth/10]); %Need a lower value here to account for any strongly-curved pipes
    else
        maxNearElemSizeZ = modelParamStruct.plateSizes.b;
    end
    line1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeByBias(biasMethod=SINGLE,';
    line1b = 'constraint=FINER, end1Edges=';
    lineE1a = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    lineE1b = int_defects_py_line_builder(modelParamStruct.edgeCoords.otherZdir.crackZone.end1, formatSpec);
    lineE1c = ', end2Edges=';
    lineE2a = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    lineE2b = int_defects_py_line_builder(modelParamStruct.edgeCoords.otherZdir.crackZone.end2, formatSpec);
    line3 = [',maxSize=',num2str(maxNearElemSizeZ),', minSize=',num2str(0.2*modelParamStruct.plateSizes.b),')'];
    lineOut = {line1a;line1b;lineE1a;lineE1b;lineE1c;lineE2a;lineE2b;line3};
end

%Edges along the z direction lower down in z
if any(strfind(lineIn,'#edgeCoords.otherZdir.lower#'))
    line1a = 'mdb.models[''Model-1''].rootAssembly.seedEdgeBySize(constraint=FINER,';
    line1b = 'deviationFactor=0.1, edges=';
    line2 = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    line3 = int_defects_py_line_builder(modelParamStruct.edgeCoords.otherZdir.lower, formatSpec);
    line4 = [', size=',num2str(modelParamStruct.plateSizes.plDepth/20),')'];
    lineOut = {line1a;line1b;line2;line3;line4};
end

%% Coordinates used for defining sets and surfaces
%zSymm condition (Set-5)
if any(strfind(lineIn,'#edgeCoords.sets.zSymm#'))
    line1a = 'mdb.models[''Model-1''].rootAssembly.Set(faces=';
    line1b = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].faces.findAt(';
    line2 = int_defects_py_line_builder(modelParamStruct.edgeCoords.sets.zSymm, formatSpec);
    line3 = ', name=''Set-5'')';
    line4 = 'mdb.models[''Model-1''].ZsymmBC(createStepName=''Initial'', localCsys=None, name=';
    line5 = '''zSymm'', region=mdb.models[''Model-1''].rootAssembly.sets[''Set-5''])';
    lineOut = {line1a;line1b;line2;line3;line4;line5};
end

%xSymm condition (Set-8) - used for biaxial loading only
if any(strfind(lineIn,'#edgeCoords.sets.xSymm#'))
    if modelParamStruct.loadingParams.biaxialFlag
        line1a = 'mdb.models[''Model-1''].rootAssembly.Set(faces=';
        line1b = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].faces.findAt(';
        line2 = int_defects_py_line_builder(modelParamStruct.edgeCoords.sets.xSymm, formatSpec);
        line3 = ', name=''Set-8'')';
        line4 = 'mdb.models[''Model-1''].XsymmBC(createStepName=''Initial'', localCsys=None, name=';
        line5 = '''xSymm'', region=mdb.models[''Model-1''].rootAssembly.sets[''Set-8''])';
        lineOut = {line1a;line1b;line2;line3;line4;line5};
    else
        lineOut = '#NOTE: No code written here. modelParamStruct.loadingParams.biaxialFlag = false.';
    end
end

%Set-6 and Set7: used for defining fixed displacement at single points
if any(strfind(lineIn,'#edgeCoords.sets.fixed#'))
    line1a = 'mdb.models[''Model-1''].rootAssembly.Set(vertices=';
    line1b = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].vertices.findAt(';
    line1c = int_defects_py_line_builder(modelParamStruct.edgeCoords.sets.fixed, formatSpec);
    line1d = ', name=''Set-6'')';
    
    line2a = 'mdb.models[''Model-1''].rootAssembly.Set(vertices=';
    line2b = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].vertices.findAt(';
    line2c = int_defects_py_line_builder(modelParamStruct.edgeCoords.sets.fixed2, formatSpec);
    line2d = ', name=''Set-7'')';
    
    line3a = 'mdb.models[''Model-1''].DisplacementBC(amplitude=UNSET, createStepName=''Initial'',';     %yFixed BC
    line3b = 'distributionType=UNIFORM, fieldName='''', localCsys=None, name=''yFixed'',';
    line3c = 'region=mdb.models[''Model-1''].rootAssembly.sets[''Set-6''], u1=UNSET, u2=SET,'; 
    line3d = 'u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)';
    
    lineOut = {line1a;line1b;line1c;line1d;line2a;line2b;line2c;line2d;line3a;line3b;line3c;line3d};
    
    if ~modelParamStruct.loadingParams.biaxialFlag
        line4a = 'mdb.models[''Model-1''].DisplacementBC(amplitude=UNSET, createStepName=''Initial'',';     %xFixed BC
        line4b = 'distributionType=UNIFORM, fieldName='''', localCsys=None, name=''xFixed'',';
        line4c = 'region=mdb.models[''Model-1''].rootAssembly.sets[''Set-6''], u1=SET, u2=UNSET,';
        line4d = 'u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)';
        
        line5a = 'mdb.models[''Model-1''].DisplacementBC(amplitude=UNSET, createStepName=''Initial'',';     %Second yFixed BC
        line5b = 'distributionType=UNIFORM, fieldName='''', localCsys=None, name=''yFixed2'',';
        line5c = 'region=mdb.models[''Model-1''].rootAssembly.sets[''Set-7''], u1=UNSET, u2=SET,';
        line5d = 'u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)';
        
       lineOut = [lineOut;line4a;line4b;line4c;line4d;line5a;line5b;line5c;line5d];
    end
end

%s_Surf-1 - used for defining coupling constraint
if any(strfind(lineIn,'#edgeCoords.sets.coupling#'))
    if any(strcmpi(modelParamStruct.loadingParams.loadingType,{'load','moment','combined'}))
        line1a = 'RP1 = mdb.models[''Model-1''].rootAssembly.ReferencePoint(point=(';
        line1b = [num2str(modelParamStruct.edgeCoords.sets.coupling(1),formatSpec),',',...
            num2str(modelParamStruct.edgeCoords.sets.coupling(2),formatSpec),',',...
            num2str(modelParamStruct.edgeCoords.sets.coupling(3),formatSpec),'))'];
        line1c = 'mdb.models[''Model-1''].rootAssembly.Set(name=''m_Set-9'', referencePoints=(';
        line1d = 'mdb.models[''Model-1''].rootAssembly.referencePoints[RP1.id], ))';
        
        line2a = 'mdb.models[''Model-1''].rootAssembly.Surface(name=''s_Surf-1'', side1Faces=';
        line2b = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].faces.findAt(';
        line2c = [int_defects_py_line_builder(modelParamStruct.edgeCoords.sets.coupling, formatSpec),')'];
               
        line4a = 'mdb.models[''Model-1''].rootAssembly.Set(name=''Set-10'', referencePoints=(';
        line4b = 'mdb.models[''Model-1''].rootAssembly.referencePoints[RP1.id], ))';
        
        %Assign different kinematic constraints and BCs depending on whether
        %we're performing limit load analysis or another analysis type. I
        %have found that while minimal BCs work for most models, for limit
        %load models it is necessary to add additional kinematic
        %constraints to avoid solver issues.
        line3a = 'mdb.models[''Model-1''].Coupling(controlPoint=';
        line3b = 'mdb.models[''Model-1''].rootAssembly.sets[''m_Set-9''], couplingType=KINEMATIC,';
        line3c = 'influenceRadius=WHOLE_SURFACE, localCsys=None, name=''Constraint-1'',';
        
        line5a = 'mdb.models[''Model-1''].DisplacementBC(amplitude=UNSET, createStepName=''Initial'','; %...which necessitates an extra BC on RP1 for the x and y directions.
        line5b = 'distributionType=UNIFORM, fieldName='''', localCsys=None, name=''xyFixed'',';
        
        if strcmpi(modelParamStruct.materialParams.type,'incremental plasticity')
            if size(modelParamStruct.materialParams.plasticTable,1)==1  %We assume that incremental plasticity models with one entry in the plasticTable are limit load models
                %Limit load models only
                line3d = 'surface=mdb.models[''Model-1''].rootAssembly.surfaces[''s_Surf-1''], u1=OFF, u2=';    %Kinematic constraint on U2, U3, UR1
                line3e = 'ON, u3=ON, ur1=ON, ur2=OFF, ur3=OFF)';
                
                line5c = 'region=mdb.models[''Model-1''].rootAssembly.sets[''Set-10''], u1=SET, u2=UNSET,';     %RP BC: U1, UR2, UR3 fixed
                line5d = 'u3=UNSET, ur1=UNSET, ur2=SET, ur3=SET)';
            else
                %Other models using incremental plasticity
                line3d = 'surface=mdb.models[''Model-1''].rootAssembly.surfaces[''s_Surf-1''], u1=OFF, u2=';    %Kinematic constraint on U3, UR1
                line3e = 'OFF, u3=ON, ur1=ON, ur2=OFF, ur3=OFF)';
                
                line5c = 'region=mdb.models[''Model-1''].rootAssembly.sets[''Set-10''], u1=SET, u2=SET,';       %RP BC: U1, U2 fixed
                line5d = 'u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)';
            end
        else
            %All other models with loadingType 'load', 'moment' or 'combined'.
            line3d = 'surface=mdb.models[''Model-1''].rootAssembly.surfaces[''s_Surf-1''], u1=OFF, u2=';    %Kinematic constraint on U3, UR1
            line3e = 'OFF, u3=ON, ur1=ON, ur2=OFF, ur3=OFF)';
            
            line5c = 'region=mdb.models[''Model-1''].rootAssembly.sets[''Set-10''], u1=SET, u2=SET,';       %RP BC: U1, U2 fixed
            line5d = 'u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)';
        end
        
        %Create loads/moments
        if strcmpi(modelParamStruct.loadingParams.loadingType,'load')
            line6a = ['mdb.models[''Model-1''].ConcentratedForce(cf3=',num2str(modelParamStruct.loadingParams.loadMagnitude,formatSpec),', createStepName=''ApplyLoad'','];
            line6b = 'distributionType=UNIFORM, field='''', localCsys=None, name=''Load-1'', region=';
            line6c = 'mdb.models[''Model-1''].rootAssembly.sets[''Set-10''])';
            lineOut = {line1a;line1b;line1c;line1d;line2a;line2b;line2c;line3a;line3b;line3c;line3d;line3e;...
                line4a;line4b;line5a;line5b;line5c;line5d;line6a;line6b;line6c};
        elseif strcmpi(modelParamStruct.loadingParams.loadingType,'moment')
            line6a = ['mdb.models[''Model-1''].Moment(cm1=',num2str(modelParamStruct.loadingParams.loadMagnitude,formatSpec),', createStepName=''ApplyLoad'','];
            line6b = 'distributionType=UNIFORM, field='''', localCsys=None, name=''Moment-1'',';
            line6c = 'region=mdb.models[''Model-1''].rootAssembly.sets[''Set-10''])';
            lineOut = {line1a;line1b;line1c;line1d;line2a;line2b;line2c;line3a;line3b;line3c;line3d;line3e;...
                line4a;line4b;line5a;line5b;line5c;line5d;line6a;line6b;line6c};
        elseif strcmpi(modelParamStruct.loadingParams.loadingType,'combined')
            line6a = ['mdb.models[''Model-1''].ConcentratedForce(cf3=',num2str(modelParamStruct.loadingParams.loadMagnitude(1),formatSpec),', createStepName=''ApplyLoad'','];
            line6b = 'distributionType=UNIFORM, field='''', localCsys=None, name=''Load-1'', region=';
            line6c = 'mdb.models[''Model-1''].rootAssembly.sets[''Set-10''])';
            line7a = ['mdb.models[''Model-1''].Moment(cm1=',num2str(modelParamStruct.loadingParams.loadMagnitude(2),formatSpec),', createStepName=''ApplyLoad'','];
            line7b = 'distributionType=UNIFORM, field='''', localCsys=None, name=''Moment-1'',';
            line7c = 'region=mdb.models[''Model-1''].rootAssembly.sets[''Set-10''])';
            lineOut = {line1a;line1b;line1c;line1d;line2a;line2b;line2c;line3a;line3b;line3c;line3d;line3e;...
                line4a;line4b;line5a;line5b;line5c;line5d;line6a;line6b;line6c;line7a;line7b;line7c};
        end
    else
        lineOut = '#NOTE: No code written here. modelParamStruct.loadingParams.loadingType is not ''load'', ''moment'' or ''combined''.';
    end
end

%Surf-2 - used for defining crack face pressure loads
if any(strfind(lineIn,'#edgeCoords.sets.crackFacePressure#'))
    if strcmpi(modelParamStruct.loadingParams.loadingType,'crackFacePressure')
        line1a = ['mdb.models[''Model-1''].ExpressionField(description='''', expression=''',modelParamStruct.loadingParams.CFPstring,''','];
        line1b = 'localCsys=None, name=''AnalyticalField-1'')';
        line2a = 'mdb.models[''Model-1''].rootAssembly.Surface(name=''Surf-2'', side1Faces=';
        line2b = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].faces.findAt(';
        line3 = [int_defects_py_line_builder(modelParamStruct.edgeCoords.sets.crackFacePressure, formatSpec),')'];
        line4 = 'mdb.models[''Model-1''].Pressure(amplitude=UNSET, createStepName=''ApplyLoad'',';
        line5 = ['distributionType=FIELD, field=''AnalyticalField-1'', magnitude=1, name='];
        line6 = '''Load-1'', region=mdb.models[''Model-1''].rootAssembly.surfaces[''Surf-2''])';
        lineOut = {line1a;line1b;line2a;line2b;line3;line4;line5;line6};
    else
        lineOut = '#NOTE: No code written here. modelParamStruct.loadingParams.loadingType is not ''crackFacePressure''.';
    end
end

%Surf-3 - used for defining pipe pressure
if any(strfind(lineIn,'#edgeCoords.sets.pipePressure#'))
    if any(strcmpi(modelParamStruct.loadingParams.loadingType,{'pipePressure','pipePressureInclCrack'}))
        line1a = 'mdb.models[''Model-1''].rootAssembly.Surface(name=''Surf-3'', side1Faces=';
        line1b = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].faces.findAt(';
        line2 = [int_defects_py_line_builder(modelParamStruct.edgeCoords.sets.pipePressure, formatSpec),')'];
        line3 = 'mdb.models[''Model-1''].Pressure(amplitude=UNSET, createStepName=''ApplyLoad'',';
        line4 = ['distributionType=UNIFORM, field='''', magnitude=',num2str(modelParamStruct.loadingParams.pressureMagnitude,formatSpec),', name='];
        line5 = '''Load-2'', region=mdb.models[''Model-1''].rootAssembly.surfaces[''Surf-3''])';
        lineOut = {line1a;line1b;line2;line3;line4;line5};
    else
        lineOut = '#NOTE: No code written here. modelParamStruct.loadingParams.loadingType is not ''pipePressure'' or ''pipePressureInclCrack''.';
    end
end

%Surf-4 - used for defining stress in 'biaxial' (i.e. x) direction
if any(strfind(lineIn,'#edgeCoords.sets.biaxialStress#'))
    if modelParamStruct.loadingParams.biaxialFlag
        line1a = 'mdb.models[''Model-1''].rootAssembly.Surface(name=''Surf-4'', side1Faces=';
        line1b = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].faces.findAt(';
        line2 = [int_defects_py_line_builder(modelParamStruct.edgeCoords.sets.biaxialLoadSurf, formatSpec),')'];
        line3 = 'mdb.models[''Model-1''].Pressure(amplitude=UNSET, createStepName=''ApplyLoad'',';
        line4 = ['distributionType=UNIFORM, field='''', magnitude=',num2str(-modelParamStruct.loadingParams.biaxialStress,formatSpec),', name='];   %Note change of sign
        line5 = '''Load-3'', region=mdb.models[''Model-1''].rootAssembly.surfaces[''Surf-4''])';
        lineOut = {line1a;line1b;line2;line3;line4;line5};
    else
        lineOut = '#NOTE: No code written here. modelParamStruct.loadingParams.biaxialFlag = false.';
    end
end

%% Crack tips
%Crack A tip line
if any(strfind(lineIn,'#edgeCoords.A.tip#'))
    line1a = 'mdb.models[''Model-1''].rootAssembly.engineeringFeatures.ContourIntegral(';
    line1b = 'collapsedElementAtTip=SINGLE_NODE, crackFront=Region(';
    line2a = 'edges=mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    line2b = [int_defects_py_line_builder(modelParamStruct.edgeCoords.A.tip, formatSpec),')'];
    
    normVectStr = sprintf('%.0f, ' , modelParamStruct.crackPlaneNormalVectors.CPNVA);
    normVectStr = normVectStr(1:end-2);
    line3 = [', crackNormal=((0.0, 0.0, 0.0), (',normVectStr,')), crackTip=Region('];
    
    line4a = 'edges=mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    line4b = [int_defects_py_line_builder(modelParamStruct.edgeCoords.A.tip, formatSpec),')'];
    line5 = ',extensionDirectionMethod=CRACK_NORMAL,midNodePosition=0.5, name=''CrackA'', symmetric=ON)';
    lineOut = {line1a;line1b;line2a;line2b;line3;line4a;line4b;line5};
end

%Crack B tip line
if any(strfind(lineIn,'#edgeCoords.B.tip#'))
    line1a = 'mdb.models[''Model-1''].rootAssembly.engineeringFeatures.ContourIntegral(';
    line1b = 'collapsedElementAtTip=SINGLE_NODE, crackFront=Region(';
    line2a = 'edges=mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    line2b = [int_defects_py_line_builder(modelParamStruct.edgeCoords.B.tip, formatSpec),')'];
    
    normVectStr = sprintf('%.0f, ' , modelParamStruct.crackPlaneNormalVectors.CPNVB);
    normVectStr = normVectStr(1:end-2);
    line3 = [', crackNormal=((0.0, 0.0, 0.0), (',normVectStr,')), crackTip=Region('];

    line4a = 'edges=mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].edges.findAt(';
    line4b = [int_defects_py_line_builder(modelParamStruct.edgeCoords.B.tip, formatSpec),')'];
    line5 = ',extensionDirectionMethod=CRACK_NORMAL,midNodePosition=0.5, name=''CrackB'', symmetric=ON)';
    lineOut = {line1a;line1b;line2a;line2b;line3;line4a;line4b;line5};
end

%% Cells for assigning mesh controls & element types
if any(strfind(lineIn,'#edgeCoords.meshControls#'))
    %Mesh controls (wedge at crack tips, tets in region surrounding cracks)
    line1a = 'mdb.models[''Model-1''].rootAssembly.setMeshControls(elemShape=WEDGE, regions=';
    line1b = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].cells.findAt(';
    line2 = int_defects_py_line_builder(modelParamStruct.edgeCoords.meshControls.tipZones1, formatSpec);
    line3 = ',)';
    line4a = 'mdb.models[''Model-1''].rootAssembly.setMeshControls(elemShape=TET, sizeGrowthRate=1, regions=';  %Note sizeGrowthRate set to lowest possible value to minimise numerical error
    line4b = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].cells.findAt(';
    line5 = int_defects_py_line_builder(modelParamStruct.edgeCoords.meshControls.aroundCracks, formatSpec);
    line6 = ',)';
    %Element types (quadratic tet in region surrounding cracks)
    line7a = 'mdb.models[''Model-1''].rootAssembly.setElementType(elemTypes=(ElemType(';
    line7b = 'elemCode=C3D20R, elemLibrary=STANDARD), ElemType(elemCode=C3D15,';
    line7c = 'elemLibrary=STANDARD), ElemType(elemCode=C3D10, elemLibrary=STANDARD)),';
    line7d = 'regions=(';
    line7e = 'mdb.models[''Model-1''].rootAssembly.instances[''PlatePartitioned-1''].cells.findAt(';
    line8 = int_defects_py_line_builder(modelParamStruct.edgeCoords.meshControls.aroundCracks, formatSpec);
    line9 = ',))';
    lineOut = {line1a;line1b;line2;line3;line4a;line4b;line5;line6;line7a;line7b;line7c;line7d;line7e;line8;line9};
end

%% Input deck lines for rotating Ellipse B instances
if any(strfind(lineIn,'#rotateEllipseB#'))
    if modelParamStruct.crackPositions.oppositeFlag     %Note that this condition is implied if Crack B is a subsurface defect.
        %NOTE - This part has been surpressed to prevent numerical
        %difficultes which occur with CAE due to the rotation. The rotation
        %is unnecessary now anyway, since I am no longer auto-trimming the
        %shell parts.
        %
        %lineOut = {'mdb.models[''Model-1''].rootAssembly.rotate(angle=180.0, axisDirection=(0.0, 0.0, ';...
        %    '1.0), axisPoint=(0.0, 0.5, 0.0), instanceList=(''TubeShellSmallB-1'', ))';...
        %    'mdb.models[''Model-1''].rootAssembly.rotate(angle=180.0, axisDirection=(0.0, 0.0, ';...
        %    '1.0), axisPoint=(0.0, 0.5, 0.0), instanceList=(''TubeShellLargeB-1'', ))';...
        %    'mdb.models[''Model-1''].rootAssembly.rotate(angle=180.0, axisDirection=(0.0, 0.0, ';...
        %    '1.0), axisPoint=(0.0, 0.5, 0.0), instanceList=(''EdgeShellB-1'', ))'};
        lineOut = '#NOTE: No code written here. Ellipse B does not need to be rotated. See int_defects_py_line_edges.m for more information.';
    else
        lineOut = '#NOTE: No code written here. Ellipse B does not need to be rotated. See int_defects_py_line_edges.m for more information.';
    end
end

%% Input deck lines for creating an additional partition in the crack region if necessary
if any(strfind(lineIn,'#crackZoneDivide#'))
    if modelParamStruct.meshTransitionPositions.divideFlag
        line1 = 'DP1 = mdb.models[''Model-1''].parts[''PlatePartitioned''].DatumPlaneByPrincipalPlane(';
        line2 = 'offset=0.0, principalPlane=YZPLANE)';
        line3 = 'mdb.models[''Model-1''].parts[''PlatePartitioned''].PartitionCellByDatumPlane(';
        line4 = 'cells=';
        line5 = 'mdb.models[''Model-1''].parts[''PlatePartitioned''].cells.findAt(';
        line6 = int_defects_py_line_builder(modelParamStruct.edgeCoords.meshControls.aroundCracks(1,:), formatSpec);
        line7 = ', datumPlane=';
        line8 = 'mdb.models[''Model-1''].parts[''PlatePartitioned''].datums[DP1.id])';
        lineOut = {line1;line2;line3;line4;line5;line6;line7;line8};
    else
        lineOut = '#NOTE - No code has been written here. It was not necessary to divide the crack region.';
    end
end

%% Input deck lines for creating a new datum plane and partitioning subsurface cracks
%Crack A
if any(strfind(lineIn,'#partitionSubsurfaceA#'))
    if modelParamStruct.crackPositions.subsurfaceAFlag
        line1 = 'mdb.models[''Model-1''].parts[''PlatePartitioned''].DatumPlaneByPrincipalPlane(';
        line2 = ['offset=',num2str(modelParamStruct.crackPositions.yA,formatSpec),', principalPlane=XZPLANE)'];
        line3 = 'mdb.models[''Model-1''].parts[''PlatePartitioned''].PartitionCellByDatumPlane(';
        line4 = 'cells=';
        line5 = 'mdb.models[''Model-1''].parts[''PlatePartitioned''].cells.findAt(';
        line6 = int_defects_py_line_builder([modelParamStruct.edgeCoords.ellipseToDivide.A], formatSpec);
        line7 = ', datumPlane= mdb.models[''Model-1''].parts[''PlatePartitioned''].datums[2])';
        lineOut = {line1;line2;line3;line4;line5;line6;line7};
    else
        lineOut = '#NOTE - No code has been written here. It is not necessary to partition Crack B (if it exists), since it is not a subsurface crack.';
    end
end
%Crack B
if any(strfind(lineIn,'#partitionSubsurfaceB#'))
    if modelParamStruct.crackPositions.subsurfaceBFlag
        line1 = 'mdb.models[''Model-1''].parts[''PlatePartitioned''].DatumPlaneByPrincipalPlane(';
        line2 = ['offset=',num2str(modelParamStruct.crackPositions.yB,formatSpec),', principalPlane=XZPLANE)'];
        line3 = 'mdb.models[''Model-1''].parts[''PlatePartitioned''].PartitionCellByDatumPlane(';
        line4 = 'cells=';
        line5 = 'mdb.models[''Model-1''].parts[''PlatePartitioned''].cells.findAt(';
        line6 = int_defects_py_line_builder([modelParamStruct.edgeCoords.ellipseToDivide.B], formatSpec);
        line7 = ', datumPlane= mdb.models[''Model-1''].parts[''PlatePartitioned''].datums[2])';
        lineOut = {line1;line2;line3;line4;line5;line6;line7};
    else
        lineOut = '#NOTE - No code has been written here. It is not necessary to partition Crack B (if it exists), since it is not a subsurface crack.';
    end
end

%% Input deck lines for auto-trimming crack (when they are surface rather than a subsurface defects).
%NOTE - Auto-trimming is currently disabled (see commented-out lines). It
%appears to be an unnecessary step for most (all?) models, and the need to
%specify the geometry object to trim (eg. geometry[3]) is annoying.
if any(strfind(lineIn,'#autoTrimA#'))
%     if ~modelParamStruct.crackPositions.subsurfaceAFlag
%         line1 = 'mdb.models[''Model-1''].sketches[''__sweep__''].autoTrimCurve(curve1=';
%         line2 = 'mdb.models[''Model-1''].sketches[''__sweep__''].geometry[3], point1=(';
%         line3 = [num2str(modelParamStruct.edgeCoords.A.autoTrim(1),formatSpec),', ',...
%             num2str(modelParamStruct.edgeCoords.A.autoTrim(2),formatSpec),...
%             '))'];
%         lineOut = {line1;line2;line3};
%     else
        lineOut = '#NOTE - No code has been written here. It was not necessary to auto-trim the sketch for Crack A sweep.';
%     end
end

if any(strfind(lineIn,'#autoTrimB#'))
%     if ~modelParamStruct.crackPositions.subsurfaceBFlag
%         line1 = 'mdb.models[''Model-1''].sketches[''__sweep__''].autoTrimCurve(curve1=';
%         line2 = 'mdb.models[''Model-1''].sketches[''__sweep__''].geometry[3], point1=(';
%         line3 = [num2str(modelParamStruct.edgeCoords.B.autoTrim(1),formatSpec),', ',...
%             num2str(modelParamStruct.edgeCoords.B.autoTrim(2),formatSpec),...
%             '))'];
%         lineOut = {line1;line2;line3};
%     else
        lineOut = '#NOTE - No code has been written here. It was not necessary to auto-trim the sketch for Crack B sweep.';
%     end
end

