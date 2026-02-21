# T-space helper for (k,d). Uses polyhedral_common TSPACE_FileFormatConversion.

# Local copies from polyhedral_common/CI_tests/common.g
RemoveFileIfExist:=function(FileName)
    if IsExistingFile(FileName) then
        RemoveFile(FileName);
    fi;
end;

WriteMatrixStream:=function(output, EXT)
    local eEXT, eVal;
    AppendTo(output, Length(EXT), " ", Length(EXT[1]), "\n");
    for eEXT in EXT
    do
        for eVal in eEXT
        do
            AppendTo(output, " ", eVal);
        od;
        AppendTo(output, "\n");
    od;
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
    Print("FindSplittingCombinatorial, Lconn=", Lconn, "\n");
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
    local dim, EXT1_red, EXT2_red, NSP1, NSP2, Sum_NSP, rnk;
    dim:=Length(EXT1[1]);
    EXT1_red:=RowReduction(EXT1).EXT;
    EXT2_red:=RowReduction(EXT2).EXT;
    Print("are_space_split |EXT1_red|=", Length(EXT1_red), " |EXT2_red|=", Length(EXT2_red), " dim=", dim, "\n");
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
        Print("det=", det, "\n");
#        Print("det=", det, " det1=", det1, " det2=", det2, "\n");
        if det = 1 then
            Print("EXT1_sat=\n");
            PrintArray(EXT1_sat);
            Print("EXT2_sat=\n");
            PrintArray(EXT2_sat);
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
    Print("nConn=", nConn, "\n");
    get_ext_from_part:=function(eP)
        return Concatenation(Lconn{eP});
    end;
    Lpart:=Filtered(Combinations([1..nConn]), x->Length(x)>0 and Length(x)<nConn);
    Print("Lpart=", Lpart, "\n");
    for ePart in Lpart
    do
        fPart:=Difference([1..nConn], ePart);
        Print("ePart=", ePart, " fPart=", fPart, "\n");
        EXT_a1:=get_ext_from_part(ePart);
        EXT_a2:=get_ext_from_part(fPart);
#        Print("EXT_a1=", EXT_a1, "\n");
#        Print("EXT_a2=", EXT_a2, "\n");
        EXT_b1:=get_space_family(EXT_a1, l_gens);
        EXT_b2:=get_space_family(EXT_a2, l_gens);
        if are_space_split(EXT_b1, EXT_b2) then
            Print("Finding a splitting with\n");
            Print("EXT_a1=", EXT_a1, "\n");
            Print("EXT_a2=", EXT_a2, "\n");
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
    AppendTo(output, "\n");
end;




get_cells:=function(k, d, index)
    local TmpDir, FileNml, FileCells, output, binary, cmd, ListCells;
    TmpDir:=DirectoryTemporary();

    FileNml:=Filename(TmpDir, "Input.nml");
    FileCells:=Filename(TmpDir, "Input.gap");

    output:=OutputTextFile(FileNml, true);
    append_initial_tspace(output, k, d);
    AppendTo(output, "&QUERIES\n");
    AppendTo(output, " FileCells = \"", FileCells, "\"\n");
    AppendTo(output, " IndexCell = ", index, "\"\n");
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
    local TmpDir, FileNml, FileEquiv, FileEquivOutput, output, binary, cmd, ListEquivOutput, eEquiv, EXT1_img, EXT2_img;
    TmpDir:=DirectoryTemporary();

    FileNml:=Filename(TmpDir, "Input.nml");
    FileEquiv:=Filename(TmpDir, "Computation");
    FileEquivOutput:=Filename(TmpDir, "Computation.output");

    output:=OutputTextFile(FileEquiv, true);
    AppendTo(output, "2\n");
    WriteMatrixStream(output, EXT1);
    WriteMatrixStream(output, EXT2);
    CloseStream(output);

    output:=OutputTextFile(FileNml, true);
    append_initial_tspace(output, k, d);
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
    eEquiv:=ListEquivOutput[1];
    if eEquiv<>fail then
        Print("eEquiv=", eEquiv, "\n");
        EXT1_img:=Set(EXT1 * eEquiv);
        EXT2_img:=Set(EXT2);
        Print("EXT1_img=", EXT1_img, "\n");
        Print("EXT2_img=", EXT2_img, "\n");
        if EXT1_img<>EXT2_img then
            Error("The equivalence is not actually one");
        fi;
    fi;
    RemoveFileIfExist(FileNml);
    RemoveFileIfExist(FileEquiv);
    RemoveFileIfExist(FileEquivOutput);
    return eEquiv;
end;


get_cell_stabilizer:=function(k, d, EXT)
    local TmpDir, FileNml, FileStab, FileStabOutput, output, binary, cmd, ListStabOutput, eStab, eGen, SetEXT, EXT_img_set;
    TmpDir:=DirectoryTemporary();

    FileNml:=Filename(TmpDir, "Input.nml");
    FileStab:=Filename(TmpDir, "Computation");
    FileStabOutput:=Filename(TmpDir, "Computation.output");

    output:=OutputTextFile(FileStab, true);
    AppendTo(output, "1\n");
    WriteMatrixStream(output, EXT);
    CloseStream(output);

    output:=OutputTextFile(FileNml, true);
    append_initial_tspace(output, k, d);
    AppendTo(output, "&QUERIES\n");
    AppendTo(output, " FileStabilizerQueries = \"", FileStab, "\"\n");
    AppendTo(output, "/\n");
    CloseStream(output);

    binary:=GetBinaryFilename("PERF_SerialPerfectComputation");
    cmd:=Concatenation(binary, " ", FileNml);
    Exec(cmd);

    if IsExistingFile(FileStabOutput)=false then
        Error("The output file is not existing. That qualifies as a fail");
    fi;

    ListStabOutput:=ReadAsFunction(FileStabOutput)();
    eStab:=ListStabOutput[1];
    SetEXT:=Set(EXT);
    for eGen in GeneratorsOfGroup(eStab)
    do
        EXT_img_set:=Set(EXT * eGen);
        if EXT_img_set<>SetEXT then
            Error("The matrix does not preserve EXT");
        fi;
    od;
    RemoveFileIfExist(FileNml);
    RemoveFileIfExist(FileStab);
    RemoveFileIfExist(FileStabOutput);
    return eStab;
end;

get_lower_cells:=function(k, d, index)
    local TmpDir, FileNml, FileListLower, output, binary, cmd, ListLower;
    TmpDir:=DirectoryTemporary();

    FileNml:=Filename(TmpDir, "Input.nml");
    FileListLower:=Filename(TmpDir, "Result");

    output:=OutputTextFile(FileNml, true);
    append_initial_tspace(output, k, d);
    AppendTo(output, "&QUERIES\n");
    AppendTo(output, " FileListLowerBoundary = \"", FileListLower, "\"\n");
    AppendTo(output, " IndexLowerBoundary = ", index, "\n");
    AppendTo(output, "/\n");
    CloseStream(output);

    binary:=GetBinaryFilename("PERF_SerialPerfectComputation");
    cmd:=Concatenation(binary, " ", FileNml);
    Exec(cmd);

    if IsExistingFile(FileListLower)=false then
        Error("The output file is not existing. That qualifies as a fail");
    fi;

    ListLower:=ReadAsFunction(FileListLower)();
    RemoveFileIfExist(FileNml);
    RemoveFileIfExist(FileListLower);
    return ListLower;
end;

get_upper_cells:=function(k, d, index)
    local TmpDir, FileNml, FileListUpper, output, binary, cmd, ListUpper;
    TmpDir:=DirectoryTemporary();

    FileNml:=Filename(TmpDir, "Input.nml");
    FileListUpper:=Filename(TmpDir, "Result");

    output:=OutputTextFile(FileNml, true);
    append_initial_tspace(output, k, d);
    AppendTo(output, "&QUERIES\n");
    AppendTo(output, " FileListUpperBoundary = \"", FileListUpper, "\"\n");
    AppendTo(output, " IndexUpperBoundary = ", index, "\n");
    AppendTo(output, "/\n");
    CloseStream(output);

    binary:=GetBinaryFilename("PERF_SerialPerfectComputation");
    cmd:=Concatenation(binary, " ", FileNml);
    Exec(cmd);

    if IsExistingFile(FileListUpper)=false then
        Error("The output file is not existing. That qualifies as a fail");
    fi;

    ListUpper:=ReadAsFunction(FileListUpper)();
    RemoveFileIfExist(FileNml);
    RemoveFileIfExist(FileListUpper);
    return ListUpper;
end;






# Check whether it is an extension from a lower dimensional cell
is_direct_extension:=function(k, d, EXT)
    local dim_space, perf_rank, EXTone_dim_cell, dim_lower_space, perf_rank_lower, index_lower, RecLower, EXTsum, eEquiv, ListCellLower;
    Print("is_direct_extension, step 1, k=", k, "\n");
    dim_space:=k*k;
    perf_rank:=get_perfect_rank(k, d, EXT);
    Print("is_direct_extension, step 2\n");
    EXTone_dim_cell:=get_cells(1, d, 0)[1].EXT;
    Print("is_direct_extension, step 3\n");
    dim_lower_space:=(k-1)*(k-1);
    perf_rank_lower:=perf_rank - 1;
    Print("is_direct_extension, step 4\n");
    index_lower:=dim_lower_space - perf_rank_lower;
    Print("is_direct_extension, step 5, index_lower=", index_lower, "\n");
    ListCellLower:=get_cells(k-1, d, index_lower);
    Print("is_direct_extension, step 6\n");
    for RecLower in ListCellLower
    do
        EXTsum:=direct_sum(RecLower.EXT, EXTone_dim_cell);
        Print("is_direct_extension, We have EXTsum\n");
#        Print("EXT=\n");
#        PrintArray(EXT);
#        Print("EXTsum=\n");
#        PrintArray(EXTsum);
#        Print("Set(EXT) = Set(EXTsum)=", Set(EXT) = Set(EXTsum), "\n");
        eEquiv:=test_equivalent_cells(k, d, EXT, EXTsum);
        Print("is_direct_extension, We have eEquiv\n");
        Print("eEquiv=", eEquiv, "\n");
        if eEquiv<>fail then
            return true;
        fi;
    od;
    return false;
end;

is_irreducible_both_method:=function(k, d, EXT)
    local rec_tspace, test1, test2, test2_not;
    Print("---------------------------------------------------------------------------------------\n");
    Print("is_irreducible_both_method, step 1 |EXT|=", Length(EXT), "\n");
    rec_tspace:=get_rec_tspace(k, d);
    Print("is_irreducible_both_method, step 2\n");
    test1:=is_irreducible_ext(rec_tspace, EXT);
    Print("is_irreducible_both_method, step 3, test1=", test1, "\n");
    test2:=is_direct_extension(k, d, EXT);
    Print("is_irreducible_both_method, step 4, test2=", test2, "\n");
    test2_not:=not test2;
    if test1<>test2_not then
        Print("EXT=\n");
        PrintArray(EXT);
        Error("Both method return different results");
    fi;
    Print("is_irreducible_both_method, step 5\n");
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

