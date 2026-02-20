# T-space helper for (k,d). Uses polyhedral_common TSPACE_FileFormatConversion.

# Local copies from polyhedral_common/CI_tests/common.g
RemoveFileIfExist:=function(FileName)
    if IsExistingFile(FileName) then
        RemoveFile(FileName);
    fi;
end;

GetBinaryFilename:=function(FileName)
    local TmpFile, list_lines;
    TmpFile:=Filename(DirectoryTemporary(), "Test.in");
    Exec("which ", FileName, " > ", TmpFile);
    list_lines:=ReadTextFile(TmpFile);
    if Length(list_lines)=0 then
        return fail;
    fi;
    return list_lines[1];
end;


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

get_rec_tspace:=function(k, d)
    local info, FileNml, FileOut, binary, output, cmd, tspace;

    info:=GetFundamentalInfo(d);
    if info.IsCorrect=false then
        Print("Discriminant d=", d, " is not valid, skipping\n");
        return fail;
    fi;

    FileNml:=Concatenation("Tspace_", String(k), "_", String(d), ".nml");
    FileOut:=Concatenation("Tspace_", String(k), "_", String(d), ".gap");

    RemoveFileIfExist(FileNml);
    RemoveFileIfExist(FileOut);

    output:=OutputTextFile(FileNml, true);
    AppendTo(output, "&TSPACE\n");
    AppendTo(output, " TypeTspace = \"", info.type_tspace, "\"\n");
    AppendTo(output, " FileLinSpa = \"unset.linspa\"\n");
    AppendTo(output, " SuperMatMethod = \"NotNeeded\"\n");
    AppendTo(output, " ListComm = \"Use_realimag\"\n");
#    AppendTo(output, " PtGroupMethod = \"Trivial\"\n");
    AppendTo(output, " PtGroupMethod = \"Compute\"\n");
    AppendTo(output, " FileListSubspaces = \"unset\"\n");
    AppendTo(output, " RealImagDim = ", k, "\n");
    AppendTo(output, " RealImagSum = ", info.eSum, "\n");
    AppendTo(output, " RealImagProd = ", info.eProd, "\n");
    AppendTo(output, "/\n");
    CloseStream(output);

    binary:=GetBinaryFilename("TSPACE_FileFormatConversion");
    cmd:=Concatenation(binary, " ", FileNml, " GAP ", FileOut);
    Exec(cmd);

    if IsExistingFile(FileOut)=false then
        Error("The output file is not existing. That qualifies as a fail");
    fi;

    tspace:=ReadAsFunction(FileOut)();

    RemoveFileIfExist(FileNml);
    RemoveFileIfExist(FileOut);

    return tspace;
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
    local dim, NSP1, NSP2, Sum_NSP, rnk;
    dim:=Length(EXT1);
    NSP1:=NullspaceMat(TransposedMat(EXT1));
    NSP2:=NullspaceMat(TransposedMat(EXT2));
    Sum_NSP:=Concatenation(NSP1, NSP2);
    if Length(Sum_NSP)=0 then
        return false;
    else
        rnk:=RankMat(Sum_NSP);
        if rnk < dim then
            return false;
        fi;
        return true;
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
    return EXTtot;
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
    get_ext_from_part:=function(eP)
        return Concatenation(Lconn{eP});
    end;
    Lpart:=Filtered(Combinations([1..nConn]), x->Length(x)>0 and Length(x)<nConn);
    for ePart in Lpart
    do
        fPart:=Difference([1..nConn], ePart);
        EXT_a1:=get_ext_from_part(ePart);
        EXT_a2:=get_ext_from_part(fPart);
        EXT_b1:=get_space_family(EXT_a1, rec_tspace.l_spanning_elements);
        EXT_b2:=get_space_family(EXT_a2, rec_tspace.l_spanning_elements);
        if are_space_split(EXT_b1, EXT_b2) then
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
        f_partA:=Concatenation(ListWithIdenticalEntries(k2, 0), e_partA);
        f_partB:=Concatenation(ListWithIdenticalEntries(k2, 0), e_partB);
        fEXT:=Concatenation(f_partA, f_partB);
        Add(EXT, fEXT);
    od;
    return EXT;
end;

append_initial_tspace:=function(output, k, d)
    local info, CacheFile;
    info:=GetFundamentalInfo(d);
    CacheFile:=Concatenation("Cache_", String(k), "_", String(d));
    AppendTo(output, "&TSPACE\n");
    AppendTo(output, " TypeTspace = \"", info.type_tspace, "\"\n");
    AppendTo(output, " FileLinSpa = \"unset.linspa\"\n");
    AppendTo(output, " SuperMatMethod = \"NotNeeded\"\n");
    AppendTo(output, " ListComm = \"Use_realimag\"\n");
#    AppendTo(output, " PtGroupMethod = \"Trivial\"\n");
    AppendTo(output, " PtGroupMethod = \"Compute\"\n");
    AppendTo(output, " FileListSubspaces = \"unset\"\n");
    AppendTo(output, " RealImagDim = ", k, "\n");
    AppendTo(output, " RealImagSum = ", info.eSum, "\n");
    AppendTo(output, " RealImagProd = ", info.eProd, "\n");
    AppendTo(output, "/\n");
    AppendTo(output, "\n");
    AppendTo(output, "&DATA\n");
    AppendTo(output, " arithmetic_T = \"gmp_rational\"\n");
    AppendTo(output, " arithmetic_Tint = \"gmp_integer\"\n");
    AppendTo(output, " OnlyWellRounded = T\n");
    AppendTo(output, " ComputeBoundary = T\n");
    AppendTo(output, " ComputeContractingHomotopy = T\n");
    AppendTo(output, " CacheFile = \"", CacheFile, "\"\n");
    AppendTo(output, "/\n");
end;




get_cells:=function(k, d, idim)
    local TmpDir, FileNml, FileCells, output, binary, cmd, ListCells;
    TmpDir:=DirectoryTemporary();

    FileNml:=Filename(TmpDir, "Input.nml");
    FileCells:=Filename(TmpDir, "Input.gap");

    output:=OutputTextFile(FileNml, true);
    append_initial_tspace(output, k, d);
    AppendTo(output, "\n");
    AppendTo(output, "&QUERIES\n");
    AppendTo(output, " FileCells = \"", FileCells, "\"\n");
    AppendTo(output, " DimCell = ", idim, "\"\n");
    AppendTo(output, "/\n");
    CloseStream(output);

    binary:=GetBinaryFilename("PERF_SerialPerfectComputation");
    cmd:=Concatenation(binary, " ", FileNml);
    Exec(cmd);

    if IsExistingFile(FileCells)=false then
        Error("The output file is not existing. That qualifies as a fail");
    fi;

    ListCells:=ReadAsFunction(FileCells)();
    RemoveFileIfExist(FileCells);
    RemoveFileIfExist(FileNml);
    return ListCells;
end;

test_equivalent_cells:=function(k, d, EXT1, EXT2)
    local TmpDir, FileNml, FileEquiv, FileEquivOutput, output, binary, cmd, ListCells, ListEquivOutput;
    TmpDir:=DirectoryTemporary();

    FileNml:=Filename(TmpDir, "Input.nml");
    FileEquiv:=Filename(TmpDir, "Computation");
    FileEquivOutput:=Filename(TmpDir, "Computation.output");

    output:=OutputTextFile(FileEquiv, true);
    AppendTo(output, "2\n");
    WriteMatrix(output, EXT1);
    WriteMatrix(output, EXT2);
    CloseStream(output);

    output:=OutputTextFile(FileNml, true);
    append_initial_tspace(output, k, d);
    AppendTo(output, "\n");
    AppendTo(output, "&QUERIES\n");
    AppendTo(output, " FileEquivalenceQueries = \"", FileEquiv, "\"\n");
    AppendTo(output, "/\n");
    CloseStream(output);

    binary:=GetBinaryFilename("PERF_SerialPerfectComputation");
    cmd:=Concatenation(binary, " ", FileNml);
    Exec(cmd);

    if IsExistingFile(FileEquivOutput)=false then
        Error("The output file is not existing. That qualifies as a fail");
    fi;

    ListEquivOutput:=ReadAsFunction(FileEquivOutput)();
    RemoveFileIfExist(FileNml);
    RemoveFileIfExist(FileEquiv);
    RemoveFileIfExist(FileEquivOutput);
    return ListEquivOutput[1];
end;

# Check whether it is an extension from a lower dimensional cell
is_direct_extension:=function(k, d, EXT)
    local dim_space, perf_rank, EXTone_dim_cell, dim_lower_space, perf_rank_lower, idim_lower, RecLower, EXTsum, eEquiv, ListCellLower;
    Print("is_direct_extension, step 1, k=", k, "\n");
    dim_space:=k*k;
    perf_rank:=get_perfect_rank(k, d, EXT);
    Print("is_direct_extension, step 2\n");
    EXTone_dim_cell:=get_cells(1, d, 0)[1].EXT;
    Print("is_direct_extension, step 3\n");
    dim_lower_space:=(k-1)*(k-1);
    perf_rank_lower:=perf_rank - 1;
    Print("is_direct_extension, step 4\n");
    idim_lower:=dim_lower_space - perf_rank_lower;
    Print("is_direct_extension, step 5, idim_lower=", idim_lower, "\n");
    ListCellLower:=get_cells(k-1, d, idim_lower);
    Print("is_direct_extension, step 6\n");
    for RecLower in ListCellLower
    do
        EXTsum:=direct_sum(RecLower.EXT, EXTone_dim_cell);
        Print("is_direct_extension, We have EXTsum\n");
        eEquiv:=test_equivalent_cells(k, d, EXT, EXTsum);
        Print("is_direct_extension, We have eEquiv\n");
        if eEquiv<>fail then
            return true;
        fi;
    od;
    return false;
end;

is_irreducible_both_method:=function(k, d, EXT)
    local rec_tspace, test1, test2, test2_not;
    Print("is_irreducible_both_method, step 1\n");
    rec_tspace:=get_rec_tspace(k, d);
    Print("is_irreducible_both_method, step 2\n");
    test1:=is_irreducible_ext(rec_tspace, EXT);
    Print("is_irreducible_both_method, step 3\n");
    test2:=is_direct_extension(k, d, EXT);
    Print("is_irreducible_both_method, step 4\n");
    test2_not:=not test2;
    if test1<>test2_not then
        Error("Both method return different results");
    fi;
    Print("is_irreducible_both_method, step 5\n");
    return test1;
end;


get_cells_with_irreducibility:=function(k, d, idim)
    local ListCells, n_cell, i_cell, is_irreducible;
    Print("get_cells_with_irreducibility, step 1\n");
    ListCells:=get_cells(k, d, idim);
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

