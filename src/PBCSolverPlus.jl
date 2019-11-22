function PBCsetDispPlus(NodeInfor,ElemInfor,OutNodeDict,OutFaceDict)

    NodeA = OutNodeDict["NodeA"];
    NodeB = OutNodeDict["NodeB"];
    NodeC = OutNodeDict["NodeC"];
    NodeD = OutNodeDict["NodeD"];

    NodeBC = OutNodeDict["NodeY1"];
    NodeAD = OutNodeDict["NodeY0"];
    NodeAB = OutNodeDict["NodeX1"];
    NodeDC = OutNodeDict["NodeX0"];


    # if (length(NodeAB)!=length(NodeDC))||(length(NodeAD)!=length(NodeBC))
    #     error("The nodes mismatch!");
    # end

    pointconsnum = 3;
    edgesXnum = length(NodeAB)*2;
    edgesYnum = length(NodeBC)*2;

    edgesXbegin = pointconsnum;
    edgesYbegin = edgesXbegin+edgesXnum;
    constaintotal = edgesYbegin+edgesYnum;

    ConstrainNodeIDArray = Array{Array,1}(undef,constaintotal);
    ConstrainNodePrArray = Array{Array,1}(undef,constaintotal);
    ConstrainNodeDFArray = Array{Array,1}(undef,constaintotal);

    ## Nodes
    # Node A
    ConstrainNodeIDArray[1]=[NodeA,NodeD];  
    ConstrainNodePrArray[1]=[1.0,-1.0];      
    ConstrainNodeDFArray[1]=[2,2];      

    # Node B
    ConstrainNodeIDArray[2]=[NodeB,NodeD,NodeA,NodeC];  
    ConstrainNodePrArray[2]=[1.0,-1.0,-1.0,-1.0];
    ConstrainNodeDFArray[2]=[1,1,1,1];

    ConstrainNodeIDArray[3]=[NodeB,NodeD,NodeC];  
    ConstrainNodePrArray[3]=[1.0,-1.0,-1.0];
    ConstrainNodeDFArray[3]=[2,2,2];

    considnow = 3;

    # Node C    

    # X face
    #=    
        @inbounds @simd for knode = 1:length(NodeAB)
            # #
            considnow = edgesXbegin+2*(knode-1)+1;
            ConstrainNodeIDArray[considnow]=[NodeAB[knode],NodeDC[knode],NodeA];
            ConstrainNodePrArray[considnow]=[1.0,-1.0,-1.0];
            ConstrainNodeDFArray[considnow]=[1,1,1];

            considnow += 1;
            ConstrainNodeIDArray[considnow]=[NodeAB[knode],NodeDC[knode]];
            ConstrainNodePrArray[considnow]=[1.0,-1.0];
            ConstrainNodeDFArray[considnow]=[2,2];

        end

        # Y face
        @inbounds @simd for knode = 1:length(NodeBC)
            # #
            considnow = edgesYbegin+2*(knode-1)+1;
            ConstrainNodeIDArray[considnow]=[NodeBC[knode],NodeAD[knode],NodeC];
            ConstrainNodePrArray[considnow]=[1.0,-1.0,-1.0];
            ConstrainNodeDFArray[considnow]=[1,1,1];

            considnow += 1;
            ConstrainNodeIDArray[considnow]=[NodeBC[knode],NodeAD[knode],NodeC];
            ConstrainNodePrArray[considnow]=[1.0,-1.0,-1.0];
            ConstrainNodeDFArray[considnow]=[2,2,2];

        end
    =#
    
    mapnodeArray,mapweightArray = LineLocMap(NodeInfor,NodeAB,NodeDC,axis=:y);
    @inbounds @simd for knode = 1:length(NodeAB)
        # #
        if length(mapnodeArray[knode])==1
            considnow += 1;
            ConstrainNodeIDArray[considnow]=[NodeAB[knode],mapnodeArray[knode][1],NodeA];
            ConstrainNodePrArray[considnow]=[1.0,-1.0,-1.0];
            ConstrainNodeDFArray[considnow]=[1,1,1];

            considnow += 1;
            ConstrainNodeIDArray[considnow]=[NodeAB[knode],mapnodeArray[knode][1]];
            ConstrainNodePrArray[considnow]=[1.0,-1.0];
            ConstrainNodeDFArray[considnow]=[2,2];
        else
            considnow += 1;
            ConstrainNodeIDArray[considnow]=[NodeAB[knode],mapnodeArray[knode][1],mapnodeArray[knode][2],NodeA];
            ConstrainNodePrArray[considnow]=[1.0,-mapweightArray[knode][1],-mapweightArray[knode][2],-1.0];
            ConstrainNodeDFArray[considnow]=[1,1,1,1];

            considnow += 1;
            ConstrainNodeIDArray[considnow]=[NodeAB[knode],mapnodeArray[knode][1],mapnodeArray[knode][2]];
            ConstrainNodePrArray[considnow]=[1.0,-mapweightArray[knode][1],-mapweightArray[knode][2]];
            ConstrainNodeDFArray[considnow]=[2,2,2];
        end

    end

    mapnodeArray,mapweightArray = LineLocMap(NodeInfor,NodeBC,NodeAD,axis=:x);
    @inbounds @simd for knode = 1:length(NodeBC)
        # #
        if length(mapnodeArray[knode])==1
            considnow += 1;
            ConstrainNodeIDArray[considnow]=[NodeBC[knode],mapnodeArray[knode][1],NodeC];
            ConstrainNodePrArray[considnow]=[1.0,-1.0,-1.0];
            ConstrainNodeDFArray[considnow]=[1,1,1];

            considnow += 1;
            ConstrainNodeIDArray[considnow]=[NodeBC[knode],mapnodeArray[knode][1],NodeC];
            ConstrainNodePrArray[considnow]=[1.0,-1.0,-1.0];
            ConstrainNodeDFArray[considnow]=[2,2,2];
        else
            considnow += 1;
            ConstrainNodeIDArray[considnow]=[NodeBC[knode],mapnodeArray[knode][1],mapnodeArray[knode][2],NodeC];
            ConstrainNodePrArray[considnow]=[1.0,-mapweightArray[knode][1],-mapweightArray[knode][2],-1.0];
            ConstrainNodeDFArray[considnow]=[1,1,1,1];

            considnow += 1;
            ConstrainNodeIDArray[considnow]=[NodeBC[knode],mapnodeArray[knode][1],mapnodeArray[knode][2],NodeC];
            ConstrainNodePrArray[considnow]=[1.0,-mapweightArray[knode][1],-mapweightArray[knode][2]],-1.0;
            ConstrainNodeDFArray[considnow]=[2,2,2,2];
        end

    end

    return ConstrainNodeIDArray,ConstrainNodePrArray,ConstrainNodeDFArray

end

function PBCsetScaleFieldPlus(NodeInfor,ElemInfor,OutNodeDict,OutFaceDict;dof=11)
    `Set the periodic boundary conditions for displacement.`

    NodeA = OutNodeDict["NodeA"];
    NodeB = OutNodeDict["NodeB"];
    NodeC = OutNodeDict["NodeC"];
    NodeD = OutNodeDict["NodeD"];

    NodeBC = OutNodeDict["NodeY1"];
    NodeAD = OutNodeDict["NodeY0"];
    NodeAB = OutNodeDict["NodeX1"];
    NodeDC = OutNodeDict["NodeX0"];

    pointconsnum = 3;
    edgesXnum = length(NodeAB)*1;
    edgesYnum = length(NodeBC)*1;

    edgesXbegin = pointconsnum;
    edgesYbegin = edgesXbegin+edgesXnum;

    constaintotal = edgesYbegin+edgesYnum;

    ConstrainNodeIDArray = Array{Array,1}(undef,constaintotal);
    ConstrainNodePrArray = Array{Array,1}(undef,constaintotal);
    ConstrainNodeDFArray = Array{Array,1}(undef,constaintotal);

    ## Nodes
    # Node A
    ConstrainNodeIDArray[1]=[NodeA,NodeD];  
    ConstrainNodePrArray[1]=[1.0,-1.0];      
    ConstrainNodeDFArray[1]=[dof,dof];      

    # Node B
    ConstrainNodeIDArray[2]=[NodeB,NodeD];  
    ConstrainNodePrArray[2]=[1.0,-1.0];
    ConstrainNodeDFArray[2]=[dof,dof];

    # Node C
    ConstrainNodeIDArray[3]=[NodeC,NodeD];
    ConstrainNodePrArray[3]=[1.0,-1.0];
    ConstrainNodeDFArray[3]=[dof,dof];

    considnow = 3;
    #=
        # X face
            
        @inbounds @simd for knode = 1:length(NodeAB)
            # #
            considnow = edgesXbegin+(knode-1)+1;
            ConstrainNodeIDArray[considnow]=[NodeAB[knode],NodeDC[knode]];
            ConstrainNodePrArray[considnow]=[1.0,-1.0];
            ConstrainNodeDFArray[considnow]=[dof,dof];

        end

        # Y face
        @inbounds @simd for knode = 1:length(NodeBC)
            # #
            considnow = edgesYbegin+(knode-1)+1;
            ConstrainNodeIDArray[considnow]=[NodeBC[knode],NodeAD[knode]];
            ConstrainNodePrArray[considnow]=[1.0,-1.0];
            ConstrainNodeDFArray[considnow]=[dof,dof];
        end
    =#
    mapnodeArray,mapweightArray = LineLocMap(NodeInfor,NodeAB,NodeDC,axis=:y);
    @inbounds @simd for knode = 1:length(NodeAB)
        # #
        if length(mapnodeArray[knode])==1
            considnow += 1;
            ConstrainNodeIDArray[considnow]=[NodeAB[knode],mapnodeArray[knode][1]];
            ConstrainNodePrArray[considnow]=[1.0,-1.0];
            ConstrainNodeDFArray[considnow]=[dof,dof];

        else
            considnow += 1;
            ConstrainNodeIDArray[considnow]=[NodeAB[knode],mapnodeArray[knode][1],mapnodeArray[knode][2]];
            ConstrainNodePrArray[considnow]=[1.0,-mapweightArray[knode][1],-mapweightArray[knode][2]];
            ConstrainNodeDFArray[considnow]=[dof,dof,dof];
        end
    end

    mapnodeArray,mapweightArray = LineLocMap(NodeInfor,NodeBC,NodeAD,axis=:x);
    @inbounds @simd for knode = 1:length(NodeBC)
        # #
        if length(mapnodeArray[knode])==1
            considnow += 1;
            ConstrainNodeIDArray[considnow]=[NodeBC[knode],mapnodeArray[knode][1]];
            ConstrainNodePrArray[considnow]=[1.0,-1.0];
            ConstrainNodeDFArray[considnow]=[dof,dof];
        else
            considnow += 1;
            ConstrainNodeIDArray[considnow]=[NodeBC[knode],mapnodeArray[knode][1],mapnodeArray[knode][2]];
            ConstrainNodePrArray[considnow]=[1.0,-mapweightArray[knode][1],-mapweightArray[knode][2]];
            ConstrainNodeDFArray[considnow]=[dof,dof,dof];
        end
    end



    return ConstrainNodeIDArray,ConstrainNodePrArray,ConstrainNodeDFArray
end

function LineLocMap(NodeInfor,NodeSetPhase,NodeSetBase;axis,tol=1e-6)
#=
    the function LineLocMap return the mapping relationship from NodeSetPhase to the base NodeSetBase.
    mapping - - - - - - - - - - - - - 
    phase:   *  *   *     *    *  *
                |  |   |     |    |  |
    base:    #    #    #     #  # #
    - - - - - - - - - - - - - - - - -
=#
    nphase = length(NodeSetPhase);
    nbase = length(NodeSetBase);

    locphaseArray = zeros(nphase);
    locbaseArray = zeros(nbase);

    if axis == :x
        slicetemp = 1;
    elseif axis == :y
        slicetemp = 2;
    else
        axis = slicetemp;
    end

    @inbounds @simd for kid = 1:nphase
        locphaseArray[kid] = NodeInfor[NodeSetPhase[kid]][slicetemp];
    end
    @inbounds @simd for kid = 1:nbase
        locbaseArray[kid] = NodeInfor[NodeSetBase[kid]][slicetemp];
    end

    mapnodeArray = Array{Array,1}(undef,nphase);
    mapweightArray = Array{Array,1}(undef,nphase);


    # find the two closest nodes in base.
    @inbounds @simd for knode = 1:nphase
        x = locphaseArray[knode];
        locdistArray = abs.(locbaseArray.-x);
        positionArray = sortperm(locdistArray);

        x = locphaseArray[knode];
        x1 = locbaseArray[positionArray[1]];
        x2 = locbaseArray[positionArray[2]];

        baseid1 = NodeSetBase[positionArray[1]];
        baseid2 = NodeSetBase[positionArray[2]];

        w1 = (x2-x)/(x2-x1);
        w2 = (x-x1)/(x2-x1);

        if w1 < tol
            nodeArray = [baseid2];
            weightArray = [1.0];
        elseif w2 < tol
            nodeArray = [baseid1];
            weightArray = [1.0];
        else
            nodeArray = [baseid1;baseid2];
            weightArray = [w1;w2];
        end

        mapnodeArray[knode] = nodeArray;
        mapweightArray[knode] = weightArray;

    end 
    return mapnodeArray,mapweightArray
        
        
end