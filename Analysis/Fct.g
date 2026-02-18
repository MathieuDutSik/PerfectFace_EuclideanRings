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

GetTspace:=function(k, d)
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
    local nVert, GRA, iVert, eEXT, eGen, fEXT, jVert;
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
    LConn:=ConnectedComponents(GRA);
    return List(Lconn, eConn->EXT{eConn});
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

IsIrreducibleFromEXT:=function(rec_tspace, EXT)
    local GRP, l_gens, eElt;
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


get_cells:=function(k, d, idim)
    local TmpDir, FileNml, FileCells, CacheFile, output, binary, cmd;
    TmpDir:=DirectoryTemporary();

    FileNml:=Filebame(TmpDir, "Conf_", String(k), "_", String(d), ".nml");
    FileCells:=Filename(TmpDir, "Conf_", String(k), "_", String(d), ".gap");
    CacheFile:=Concatenation("Cache_", String(k), "_", String(d));

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
    AppendTo(output, "\n");
    AppendTo(output, "&DATA\n");
    AppendTo(output, " arithmetic_T = \"gmp_rational\"\n");
    AppendTo(output, " arithmetic_Tint = \"gmp_integer\"\n");
    AppendTo(output, " OnlyWellRounded = T\n");
    AppendTo(output, " ComputeBoundary = T\n");
    AppendTo(output, " ComputeContractingHomotopy = T\n");
    AppendTo(output, "/\n");
    AppendTo(output, "\n");
    AppendTo(output, "&QUERIES\n");
    AppendTo(output, " CacheFile = \"", CacheFile, "\"\n");
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
