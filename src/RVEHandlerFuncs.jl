function AllFaceGet(NodeDict,ElemDict)
"""
"""
    ElemIDArray = collect(keys(ElemDict))
    elemtotal = length(ElemIDArray)

    Face_1_all = zeros(Int64,elemtotal,4);
    Face_1_all_Normal = zeros(Float64,elemtotal,2);
    Face_2_all = zeros(Int64,elemtotal,4);
    Face_2_all_Normal = zeros(Float64,elemtotal,2);
    Face_3_all = zeros(Int64,elemtotal,4);
    Face_3_all_Normal = zeros(Float64,elemtotal,2);
    Face_4_all = zeros(Int64,elemtotal,4);
    Face_4_all_Normal = zeros(Float64,elemtotal,2);

    Threads.@threads for kelem = 1:elemtotal
        #
        elemID = ElemIDArray[kelem];
        #=
           4          3
            C - 3 - B
            |       |
            4       2
            |       |
            D - 1 - A
           1          2
        =#

        node1 = ElemDict[elemID][1];
        node2 = ElemDict[elemID][2];
        node3 = ElemDict[elemID][3];
        node4 = ElemDict[elemID][4];

        Face_1_all[kelem,:] = [elemID 1 node1 node2];
        Face_1_all_Normal[kelem,:] = FaceNormal([node1 node2],NodeDict)[:];
        Face_2_all[kelem,:] = [elemID 2 node2 node3];
        Face_2_all_Normal[kelem,:] = FaceNormal([node2 node3],NodeDict)[:];
        Face_3_all[kelem,:] = [elemID 3 node3 node4];
        Face_3_all_Normal[kelem,:] = FaceNormal([node3 node4],NodeDict)[:];
        Face_4_all[kelem,:] = [elemID 4 node4 node1];
        Face_4_all_Normal[kelem,:] = FaceNormal([node4 node1],NodeDict)[:];

    end

    Face_all = [Face_1_all;Face_2_all;Face_3_all;Face_4_all];
    Face_all_Normal = [Face_1_all_Normal;Face_2_all_Normal;Face_3_all_Normal;
                       Face_4_all_Normal];
    
    return Face_all,Face_all_Normal

end

function FaceNormal(Face::Array,NodeDict::Dict)
`function FaceNormal: return normailized normal vector of a element surface`
#   please input the node in a correct order to make sure that the vector return
# could point out to out.

    VectorF = NodeDict[Face[2]][1:2]-NodeDict[Face[1]][1:2];
    VectorF = VectorF/norm(VectorF);
    VectorNN = [VectorF[2];-VectorF[1]];
    return VectorNN

end

function OuterFaceGet(NodeDict::Dict,ElemDict::Dict)
` Return the outer surface and normal vectors of them.`
    Face_all,Face_all_Normal = AllFaceGet(NodeDict,ElemDict);
    rowsort!(Face_all,[3,4]);

    Facetemp = [Face_all[:,3:4] collect(1:size(Face_all,1))];
    Facetemp = sortslices(Facetemp,dims=1);
    Difftemp = abs.(Facetemp[2:end,1]-Facetemp[1:end-1,1])+abs.(Facetemp[2:end,2]-Facetemp[1:end-1,2]);
    ZeroLoctemp = findall(iszero.(Difftemp));
    RepeatLoc = [ZeroLoctemp.+1;ZeroLoctemp];

    FacetempSingle = Facetemp[setdiff(collect(1:size(Face_all,1)),RepeatLoc),:];
    Face_Out = Face_all[FacetempSingle[:,end],:];
    Face_Out_Normal = Face_all_Normal[FacetempSingle[:,end],:];

    return Face_Out,Face_Out_Normal
end

function rowsort!(InputArray::Array,SortRange::Array = collect(1:size(InputArray,2)))
`The function rowsort! modify serval colume in each row by a specified colume number.`
    ArraytobeSorted = InputArray[:,SortRange]
    ArraybeSorted = sort(ArraytobeSorted,dims=2)
    InputArray[:,SortRange] = ArraybeSorted
end

function OuterNodesPicker(NodeDict::Dict,ElemDict::Dict)
`This function return the charicter geometry elements ,e.g. vertex, edge and surface, nodes. `
    Face_Out,Face_Out_Normal = OuterFaceGet(NodeDict::Dict,ElemDict::Dict);

    FaceOuter_direction_Judger = Int.(round.(Face_Out_Normal[:,1]+2.0*Face_Out_Normal[:,2]));
    FaceXP = Face_Out[findall(FaceOuter_direction_Judger.== 1),:];
    FaceXN = Face_Out[findall(FaceOuter_direction_Judger.==-1),:];
    FaceYP = Face_Out[findall(FaceOuter_direction_Judger.== 2),:];
    FaceYN = Face_Out[findall(FaceOuter_direction_Judger.==-2),:];

    OutFaceDict = Dict();
    OutFaceDict["FaceXP"] = FaceXP;
    OutFaceDict["FaceXN"] = FaceXN;
    OutFaceDict["FaceYP"] = FaceYP;
    OutFaceDict["FaceYN"] = FaceYN;

    NodeXP = unique([FaceXP[:,3];FaceXP[:,4]]);
    NodeXN = unique([FaceXN[:,3];FaceXN[:,4]]);
    NodeYP = unique([FaceYP[:,3];FaceYP[:,4]]);
    NodeYN = unique([FaceYN[:,3];FaceYN[:,4]]);

    # Point
    NodeA=intersect(intersect(NodeXP,NodeYN),NodeZN);
    NodeB=intersect(intersect(NodeXP,NodeYP),NodeZN);
    NodeC=intersect(intersect(NodeXN,NodeYP),NodeZN);
    NodeD=intersect(intersect(NodeXN,NodeYN),NodeZN);
   

    # Edge
    NodeBC=setdiff(intersect(NodeYP,NodeZN),[NodeB;NodeC]);
    NodeAD=setdiff(intersect(NodeYN,NodeZN),[NodeA;NodeD]);
    NodeAB=setdiff(intersect(NodeXP,NodeZN),[NodeA;NodeB]);
    NodeDC=setdiff(intersect(NodeXN,NodeZN),[NodeD;NodeC]);

    # Surface
    NodeX0=setdiff(NodeXN,[NodeC;NodeD;NodeH;NodeG;NodeDC;NodeHG;NodeHD;NodeGC]);
    NodeX1=setdiff(NodeXP,[NodeA;NodeB;NodeE;NodeF;NodeAB;NodeFB;NodeEF;NodeEA]);

    NodeY0=setdiff(NodeYN,[NodeA;NodeE;NodeH;NodeD;NodeEA;NodeHE;NodeHD;NodeAD]);
    NodeY1=setdiff(NodeYP,[NodeB;NodeF;NodeG;NodeC;NodeFB;NodeGF;NodeGC;NodeBC]);


    OutNodeDict = Dict();
    OutNodeDict["NodeA"] = NodeA[1];
    OutNodeDict["NodeB"] = NodeB[1];
    OutNodeDict["NodeC"] = NodeC[1];
    OutNodeDict["NodeD"] = NodeD[1];

    OutNodeDict["NodeBC"] = sortlinenodes(NodeInfor,NodeBC,:x);
    OutNodeDict["NodeAD"] = sortlinenodes(NodeInfor,NodeAD,:x);
    OutNodeDict["NodeAB"] = sortlinenodes(NodeInfor,NodeAB,:y);
    OutNodeDict["NodeDC"] = sortlinenodes(NodeInfor,NodeDC,:y);
    
    OutNodeDict["NodeXP"] = NodeXP;
    OutNodeDict["NodeXN"] = NodeXN;
    OutNodeDict["NodeYP"] = NodeYP;
    OutNodeDict["NodeYN"] = NodeYN;

    OutNodeDict["NodeX0"] = NodeX0;
    OutNodeDict["NodeX1"] = NodeX1;
    OutNodeDict["NodeY0"] = NodeY0;
    OutNodeDict["NodeY1"] = NodeY1;
    
    return OutNodeDict,OutFaceDict
end

function sortlinenodes(nodesDict::Dict,NodeArray::Array,axis::Symbol)
`The function sortlinenodes sort nodes along increase coordinate.`
# input:
# nodesDict: dictionary records nodes information
# NodeArray: the array given for sorting
# axis: sort along which direction
# output:
# NodeArrayNew: the sorted array.

    # if axis in [:x,:X]
    #     CoordArray = [nodesDict[key][1] for key in NodeArray];
    # elseif axis in [:y,:Y]
    #     CoordArray = [nodesDict[key][2] for key in NodeArray];
    # elseif axis in [:z,:Z]
    #     CoordArray = [nodesDict[key][3] for key in NodeArray];
    # end

    # CoordandIDArray = [CoordArray NodeArray];
    if axis in [:x,:X]
        sortnum = 1;
    elseif axis in [:y,:Y]
        sortnum = 2;
    elseif axis in [:z,:Z]
        sortnum = 3;
    end 

    CoordandIDArray = zeros(length(NodeArray),2);
    for knum = 1:length(NodeArray)
        CoordandIDArray[knum,:] = [nodesDict[NodeArray[knum]][sortnum];NodeArray[knum]];
    end

    CoordandIDArraynew = sortslices(CoordandIDArray,dims=1);
    NodeArrayNew = CoordandIDArraynew[:,2];
    NodeArrayNew = Int64.(NodeArrayNew);

    return NodeArrayNew
end