function PBCsetScaleField(NodeInfor,ElemInfor,OutNodeDict,OutFaceDict;dof=11)
    `Set the periodic boundary conditions for displacement.`

    NodeA = OutNodeDict["NodeA"];
    NodeB = OutNodeDict["NodeB"];
    NodeC = OutNodeDict["NodeC"];
    NodeD = OutNodeDict["NodeD"];

    NodeBC = OutNodeDict["NodeBC"];
    NodeAD = OutNodeDict["NodeAD"];
    NodeAB = OutNodeDict["NodeAB"];
    NodeDC = OutNodeDict["NodeDC"];

    if (length(NodeAB)!=length(NodeDC))||(length(NodeAD)!=length(NodeBC))
        error("The nodes mismatch!");
    end

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


    return ConstrainNodeIDArray,ConstrainNodePrArray,ConstrainNodeDFArray
end

function PBCsetScaleField(NodeInfor,ElemInfor;dof=11)
"""
"""
    OutNodeDict,OutFaceDict = OuterNodesPicker(NodeInfor,ElemInfor);
    ConstrainNodeIDArray,ConstrainNodePrArray,ConstrainNodeDFArray = 
        PBCsetScaleField(NodeInfor,ElemInfor,OutNodeDict,OutFaceDict ;dof=11)
end

function PBCsetDisp(NodeInfor,ElemInfor,OutNodeDict,OutFaceDict)
    `Set the periodic boundary conditions for displacement.`

    NodeA = OutNodeDict["NodeA"];
    NodeB = OutNodeDict["NodeB"];
    NodeC = OutNodeDict["NodeC"];
    NodeD = OutNodeDict["NodeD"];

    NodeBC = OutNodeDict["NodeY1"];
    NodeAD = OutNodeDict["NodeY0"];
    NodeAB = OutNodeDict["NodeX1"];
    NodeDC = OutNodeDict["NodeX0"];

    if (length(NodeAB)!=length(NodeDC))||(length(NodeAD)!=length(NodeBC))
        error("The nodes mismatch!");
    end

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

    # Node C    

    # X face
        
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


    return ConstrainNodeIDArray,ConstrainNodePrArray,ConstrainNodeDFArray
end

function PBCsetDisp(NodeInfor,ElemInfor)
"""
"""
    OutNodeDict,OutFaceDict = OuterNodesPicker(NodeInfor,ElemInfor);
    ConstrainNodeIDArray,ConstrainNodePrArray,ConstrainNodeDFArray = 
        PBCsetDisp(NodeInfor,ElemInfor,OutNodeDict,OutFaceDict);
end
    