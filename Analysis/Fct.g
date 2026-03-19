# T-space helper for (k,d). Uses polyhedral_common TSPACE_FileFormatConversion.
Read("common.g");
Read("access_points.g");


GetFundamentalInfo:=function(d)
  local res, IsCorrect, eSum, eProd, Dval, eQuot, type_tspace;
  res:=d mod 4;
  IsCorrect:=false;
  eSum:=0;
  eProd:=0;
  if res=0 then
    eSum:=0;
    eProd:=-d/4;
    Dval:=-eProd;
    IsCorrect:=true;
    if IsInt((Dval-1)/4)=true or IsInt(Dval/4)=true then
      IsCorrect:=false;
    fi;
  fi;
  if res=1 then
    eQuot:=(1-d)/4;
    eSum:=1;
    eProd:=eQuot;
    IsCorrect:=true;
  fi;
  if d < 0 then
      type_tspace:="ImagQuad";
  else
      type_tspace:="RealQuad";
  fi;
  return rec(eSum:=eSum, eProd:=eProd, IsCorrect:=IsCorrect, type_tspace:=type_tspace);
end;

FindSplittingCombinatorial:=function(l_gens, EXT)
    local nVert, GRA, iVert, eEXT, eGen, fEXT, jVert, Lconn;
    nVert:=Length(EXT);
    GRA:=NullGraph(Group(()), nVert);
    for iVert in [1..nVert]
    do
        eEXT:=EXT[iVert];
        for eGen in l_gens
        do
            fEXT:=eEXT * eGen;
            jVert:=Position(EXT, fEXT);
            if jVert<>fail and jVert<>iVert then
                AddEdgeOrbit(GRA, [iVert, jVert]);
                AddEdgeOrbit(GRA, [jVert, iVert]);
            fi;
        od;
    od;
    Lconn:=ConnectedComponents(GRA);
#    Print("FindSplittingCombinatorial, Lconn=", Lconn, "\n");
    return List(Lconn, eConn->EXT{eConn});
end;


get_imag_space:=function(n, d)
    local info, ListMat, i, j, eMat, eVal;
    info:=GetFundamentalInfo(d);
    ListMat:=[];
    for i in [1..n]
    do
        for j in [i..n]
        do
            eMat:=NullMat(2*n, 2*n);
            eMat[i][j] := 1;
            eMat[j][i] := 1;
            eVal:= info.eSum / 2;
            eMat[n + i][j] := eVal;
            eMat[n + j][i] := eVal;
            eMat[i][n + j] := eVal;
            eMat[j][n + i] := eVal;
            eMat[n + i][n + j] := info.eProd;
            eMat[n + j][n + i] := info.eProd;
            Add(ListMat, eMat);
        od;
    od;
    for i in [1..n]
    do
        for j in [i+1..n]
        do
            eMat:=NullMat(2*n, 2*n);
            eMat[n + i][j] := 1;
            eMat[j][n + i] := 1;
            eMat[n + j][i] := -1;
            eMat[i][n + j] := -1;
            Add(ListMat, eMat);
        od;
    od;
    return ListMat;
end;

get_perfect_rank:=function(n, d, EXT)
    local k_read, ListMat, n_mat, n_row, TotMat, i_mat, i_row, eEXT, scal;
    k_read:=Length(EXT[1])/2;
    if k_read<>n then
        Error("get_perfect_rank: Incorrect entry");
    fi;
    ListMat:=get_imag_space(n, d);
    n_mat:=Length(ListMat);
    n_row:=Length(EXT);
    TotMat:=NullMat(n_mat, n_row);
    for i_mat in [1..n_mat]
    do
        for i_row in [1..n_row]
        do
            eEXT:=EXT[i_row];
            scal:=eEXT * ListMat[i_mat] * eEXT;
            TotMat[i_mat][i_row]:=scal;
        od;
    od;
    return RankMat(TotMat);
end;




are_space_split:=function(EXT1, EXT2)
    local dim, EXT1_red, EXT2_red, get_int_nsp, EXT1_sat, EXT2_sat, EXT12_sum, NSP1, NSP2, Sum_NSP, rnk, det;
    dim:=Length(EXT1[1]);
    EXT1_red:=RowReduction(EXT1).EXT;
    EXT2_red:=RowReduction(EXT2).EXT;
#    Print("are_space_split |EXT1_red|=", Length(EXT1_red), " |EXT2_red|=", Length(EXT2_red), " dim=", dim, "\n");
    get_int_nsp:=function(EXTin)
        if Length(EXTin)=0 then
            return IdentityMat(dim);
        fi;
        if Length(EXTin)=dim then
            return [];
        fi;
        return NullspaceIntMat(TransposedMat(EXTin));
    end;
    NSP1:=get_int_nsp(EXT1_red);
    NSP2:=get_int_nsp(EXT2_red);
    Sum_NSP:=Concatenation(NSP1, NSP2);
    if Length(Sum_NSP)=0 then
        return false;
    else
        rnk:=RankMat(Sum_NSP);
        if rnk < dim then
            return false;
        fi;
        EXT1_sat:=get_int_nsp(NSP1);
        EXT2_sat:=get_int_nsp(NSP2);
#        Expr1:=List(EXT1, x->SolutionIntMat(EXT1_sat, x));
#        Expr2:=List(EXT2, x->SolutionIntMat(EXT2_sat, x));
#        det1:=AbsInt(DeterminantMat(BaseIntMat(Expr1)));
#        det2:=AbsInt(DeterminantMat(BaseIntMat(Expr2)));
        EXT12_sum:=Concatenation(EXT1_sat, EXT2_sat);
        det:=AbsInt(DeterminantMat(EXT12_sum));
#        Print("det=", det, "\n");
#        Print("det=", det, " det1=", det1, " det2=", det2, "\n");
        if det = 1 then
#            Print("EXT1_sat=\n");
#            PrintArray(EXT1_sat);
#            Print("EXT2_sat=\n");
#            PrintArray(EXT2_sat);
            return true;
        fi;
        return false;
    fi;
end;



get_space_family:=function(EXT, l_spanning_elements)
    local EXTtot, eEXT, e_spann, fEXT;
    EXTtot:=[];
    for eEXT in EXT
    do
        Add(EXTtot, eEXT);
        for e_spann in l_spanning_elements
        do
            fEXT:=eEXT * e_spann;
            Add(EXTtot, fEXT);
        od;
    od;
    return Set(EXTtot);
end;

is_irreducible_ext:=function(rec_tspace, EXT)
    local GRP, l_gens, eElt, EXT_a1, EXT_a2, EXT_b1, EXT_b2, Lpart, Lconn, nConn, ePart, fPart, get_ext_from_part;
    GRP:=Group(rec_tspace.PtStabGens);
    l_gens:=rec_tspace.l_spanning_elements;
    for eElt in GRP
    do
        Add(l_gens, eElt);
    od;
    Lconn:=FindSplittingCombinatorial(l_gens, EXT);
    nConn:=Length(Lconn);
#    Print("nConn=", nConn, "\n");
    get_ext_from_part:=function(eP)
        return Concatenation(Lconn{eP});
    end;
    Lpart:=Filtered(Combinations([1..nConn]), x->Length(x)>0 and Length(x)<nConn);
#    Print("Lpart=", Lpart, "\n");
    for ePart in Lpart
    do
        fPart:=Difference([1..nConn], ePart);
#        Print("ePart=", ePart, " fPart=", fPart, "\n");
        EXT_a1:=get_ext_from_part(ePart);
        EXT_a2:=get_ext_from_part(fPart);
#        Print("EXT_a1=", EXT_a1, "\n");
#        Print("EXT_a2=", EXT_a2, "\n");
        EXT_b1:=get_space_family(EXT_a1, l_gens);
        EXT_b2:=get_space_family(EXT_a2, l_gens);
        if are_space_split(EXT_b1, EXT_b2) then
#            Print("Finding a splitting with\n");
#            Print("EXT_a1=", EXT_a1, "\n");
#            Print("EXT_a2=", EXT_a2, "\n");
            return false;
        fi;
    od;
    return true;
end;


direct_sum:=function(EXT1, EXT2)
    local k1, k2, EXT, eEXT1, e_partA, e_partB, f_partA, f_partB, fEXT, eEXT2;
    k1:=Length(EXT1[1]) / 2;
    k2:=Length(EXT2[1]) / 2;
    EXT:=[];
    for eEXT1 in EXT1
    do
        e_partA:=eEXT1{[1..k1]};
        e_partB:=eEXT1{[k1+1..2*k1]};
        f_partA:=Concatenation(e_partA, ListWithIdenticalEntries(k2, 0));
        f_partB:=Concatenation(e_partB, ListWithIdenticalEntries(k2, 0));
        fEXT:=Concatenation(f_partA, f_partB);
        Add(EXT, fEXT);
    od;
    for eEXT2 in EXT2
    do
        e_partA:=eEXT2{[1..k2]};
        e_partB:=eEXT2{[k2+1..2*k2]};
        f_partA:=Concatenation(ListWithIdenticalEntries(k1, 0), e_partA);
        f_partB:=Concatenation(ListWithIdenticalEntries(k1, 0), e_partB);
        fEXT:=Concatenation(f_partA, f_partB);
        Add(EXT, fEXT);
    od;
    return EXT;
end;


get_cells:=function(k, d, index)
    local desc, only_well_rounded;
    only_well_rounded:=true;
    desc:=GenerateTspaceDescription_imag_quad(k, d, only_well_rounded);
    return PERFCOMP_get_cells(desc, index);
end;

test_equivalent_cells:=function(k, d, EXT1, EXT2)
    local desc, only_well_rounded;
    only_well_rounded:=true;
    desc:=GenerateTspaceDescription_imag_quad(k, d, only_well_rounded);
    return PERFCOMP_test_equivalence(desc, EXT1, EXT2);
end;


get_cell_stabilizer:=function(k, d, EXT)
    local desc, only_well_rounded;
    only_well_rounded:=true;
    desc:=GenerateTspaceDescription_imag_quad(k, d, only_well_rounded);
    return PERFCOMP_stabilizer(desc, EXT);
end;

get_lower_cells:=function(k, d, index)
    local desc, only_well_rounded;
    only_well_rounded:=true;
    desc:=GenerateTspaceDescription_imag_quad(k, d, only_well_rounded);
    return PERFCOMP_get_lower_cells(desc, index);
end;





get_upper_cells:=function(k, d, index)
    local desc, only_well_rounded;
    only_well_rounded:=true;
    desc:=GenerateTspaceDescription_imag_quad(k, d, only_well_rounded);
    return PERFCOMP_get_upper_cells(desc, index);
end;


my_get_rec_tspace:=function(k, d)
    local desc, only_well_rounded;
    only_well_rounded:=true;
    desc:=GenerateTspaceDescription_imag_quad(k, d, only_well_rounded);
    return get_rec_tspace(desc);
end;






# Check whether it is an extension from a lower dimensional cell
is_direct_extension:=function(k, d, EXT)
    local dim_space, perf_rank, EXTone_dim_cell, dim_lower_space, perf_rank_lower, index_lower, RecLower, EXTsum, eEquiv, ListCellLower;
#    Print("is_direct_extension, step 1, k=", k, "\n");
    dim_space:=k*k;
    perf_rank:=get_perfect_rank(k, d, EXT);
#    Print("is_direct_extension, step 2\n");
    EXTone_dim_cell:=get_cells(1, d, 0)[1].EXT;
#    Print("is_direct_extension, step 3\n");
    dim_lower_space:=(k-1)*(k-1);
    perf_rank_lower:=perf_rank - 1;
#    Print("is_direct_extension, step 4\n");
    index_lower:=dim_lower_space - perf_rank_lower;
#    Print("is_direct_extension, step 5, index_lower=", index_lower, "\n");
    ListCellLower:=get_cells(k-1, d, index_lower);
#    Print("is_direct_extension, step 6\n");
    for RecLower in ListCellLower
    do
        EXTsum:=direct_sum(RecLower.EXT, EXTone_dim_cell);
#        Print("is_direct_extension, We have EXTsum\n");
#        Print("EXT=\n");
#        PrintArray(EXT);
#        Print("EXTsum=\n");
#        PrintArray(EXTsum);
#        Print("Set(EXT) = Set(EXTsum)=", Set(EXT) = Set(EXTsum), "\n");
        eEquiv:=test_equivalent_cells(k, d, EXT, EXTsum);
#        Print("is_direct_extension, We have eEquiv\n");
#        Print("eEquiv=", eEquiv, "\n");
        if eEquiv<>fail then
            return true;
        fi;
    od;
    return false;
end;

is_irreducible_both_method:=function(k, d, EXT)
    local rec_tspace, test1, test2, test2_not;
#    Print("---------------------------------------------------------------------------------------\n");
#    Print("is_irreducible_both_method, step 1 |EXT|=", Length(EXT), "\n");
    rec_tspace:=my_get_rec_tspace(k, d);
#    Print("is_irreducible_both_method, step 2\n");
    test1:=is_irreducible_ext(rec_tspace, EXT);
#    Print("is_irreducible_both_method, step 3, test1=", test1, "\n");
    test2:=is_direct_extension(k, d, EXT);
#    Print("is_irreducible_both_method, step 4, test2=", test2, "\n");
    test2_not:=not test2;
    if test1<>test2_not then
        Print("EXT=\n");
        PrintArray(EXT);
        Error("Both method return different results");
    fi;
#    Print("is_irreducible_both_method, step 5\n");
    return test1;
end;


get_cells_with_irreducibility:=function(k, d, index)
    local ListCells, n_cell, i_cell, is_irreducible;
    Print("get_cells_with_irreducibility, step 1\n");
    ListCells:=get_cells(k, d, index);
    Print("get_cells_with_irreducibility, step 2\n");
    n_cell:=Length(ListCells);
    for i_cell in [1..n_cell]
    do
        is_irreducible:=is_irreducible_both_method(k, d, ListCells[i_cell].EXT);
        ListCells[i_cell].is_irreducible:=is_irreducible;
    od;
    Print("get_cells_with_irreducibility, step 3\n");
    return ListCells;
end;

get_graph:=function(ListEXT1, ListEXT2)
    local n_vert1, n_vert2, GRA, i1, i2;
    n_vert1:=Length(ListEXT1);
    n_vert2:=Length(ListEXT2);
    GRA:=NullGraph(Group(()), n_vert1 + n_vert2);
    for i1 in [1..n_vert1]
    do
        for i2 in [1..n_vert2]
        do
            if IsSubset(ListEXT2[i2], ListEXT1[i1]) then
                AddEdgeOrbit(GRA, [i1, n_vert1+i2]);
                AddEdgeOrbit(GRA, [n_vert1+i2, i1]);
            fi;
        od;
    od;
    return GRA;
end;


get_upper_graphs:=function(k, d, index)
    local list_upp0, list_upp1, list_cell1, ListGraphInfo, ListEXT1, ListEXT2, ListEXT1_iOrb, ListEXT2_iOrb, pos, eMap, fMap, i, len, EXT, EXT2, GRA, graph_info, i_ent;
    list_upp0:=get_upper_cells(k, d, index);
    list_upp1:=get_upper_cells(k, d, index-1);
    list_cell1:=get_cells(k, d, index-1);
    ListGraphInfo:=[];
    for i in [1..Length(list_upp0)]
    do
        ListEXT1:=List(list_upp0[i].ListEXT, Set);
        ListEXT2:=Set([]);
        for eMap in list_upp0[i].ListMap
        do
            len:=Length(list_upp1[eMap.jOrb].ListEXT);
            for i_ent in [1..len]
            do
                EXT:=list_upp1[eMap.jOrb].ListEXT[i_ent];
                EXT2:=Set(EXT * eMap.M);
                AddSet(ListEXT2, EXT2);
            od;
        od;
        ListEXT1_iOrb:=[];
        ListEXT2_iOrb:=ListWithIdenticalEntries(Length(ListEXT2), -1);
        for eMap in list_upp0[i].ListMap
        do
            Add(ListEXT1_iOrb, eMap.jOrb);
            len:=Length(list_upp1[eMap.jOrb].ListEXT);
            for i_ent in [1..len]
            do
                EXT:=list_upp1[eMap.jOrb].ListEXT[i_ent];
                fMap:=list_upp1[eMap.jOrb].ListMap[i_ent];
                EXT2:=Set(EXT * eMap.M);
                pos:=Position(ListEXT2, EXT2);
                if pos<>fail then
                    ListEXT2_iOrb[pos]:=fMap.jOrb;
                fi;
            od;
        od;
        GRA:=get_graph(ListEXT1, ListEXT2);
        graph_info:=rec(Graph:=GRA, ListEXT1:=ListEXT1, ListEXT2:=ListEXT2, ListEXT1_iOrb:=ListEXT1_iOrb, ListEXT2_iOrb:=ListEXT2_iOrb);
        Add(ListGraphInfo, graph_info);
    od;
    return ListGraphInfo;
end;

get_lower_graphs:=function(k, d, index)
    local list_low0, list_low1, list_cell1, ListGraphInfo, ListEXT1, ListEXT2, ListEXT1_iOrb, ListEXT2_iOrb, pos, eBnd, fBnd, i, len, EXT, EXT1, GRA, graph_info, i_ent;
    list_low0:=get_lower_cells(k, d, index);
    list_low1:=get_lower_cells(k, d, index+1);
    list_cell1:=get_cells(k, d, index-1);
    ListGraphInfo:=[];
    for i in [1..Length(list_low0)]
    do
        ListEXT2:=List(list_low0[i].ListBnd, x->Set(x.EXT));
        ListEXT1:=Set([]);
        for eBnd in list_low0[i].ListBnd
        do
            len:=Length(list_low1[eBnd.jOrb].ListBnd);
            for i_ent in [1..len]
            do
                EXT:=list_low1[eBnd.jOrb].ListBnd[i_ent].EXT;
                EXT1:=Set(EXT * eBnd.M);
                AddSet(ListEXT1, EXT1);
            od;
        od;
        ListEXT1_iOrb:=ListWithIdenticalEntries(Length(ListEXT2), -1);
        ListEXT2_iOrb:=[];
        for eBnd in list_low0[i].ListBnd
        do
            Add(ListEXT2_iOrb, eBnd.jOrb);
            len:=Length(list_low1[eBnd.jOrb].ListBnd);
            for i_ent in [1..len]
            do
                fBnd:=list_low1[eBnd.jOrb].ListBnd[i_ent];
                EXT:=fBnd.EXT;
                EXT1:=Set(fBnd.EXT * eBnd.M);
                pos:=Position(ListEXT1, EXT1);
                ListEXT1_iOrb[pos]:=fBnd.jOrb;
            od;
        od;
        GRA:=get_graph(ListEXT1, ListEXT2);
        graph_info:=rec(Graph:=GRA, ListEXT1:=ListEXT1, ListEXT2:=ListEXT2, ListEXT1_iOrb:=ListEXT1_iOrb, ListEXT2_iOrb:=ListEXT2_iOrb);
        Add(ListGraphInfo, graph_info);
    od;
    return ListGraphInfo;
end;

