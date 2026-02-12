# T-space helper for (k,d). Uses polyhedral_common TSPACE_FileFormatConversion.

# Local copies from polyhedral_common/CI_tests/common.g
RemoveFileIfExist:=function(FileName)
    if IsExistingFile(FileName) then
        RemoveFile(FileName);
    fi;
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
    local info, FileNml, FileOut, output, cmd, tspace;

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
    AppendTo(output, " PtGroupMethod = \"Trivial\"\n");
    AppendTo(output, " FileListSubspaces = \"unset\"\n");
    AppendTo(output, " RealImagDim = ", k, "\n");
    AppendTo(output, " RealImagSum = ", info.eSum, "\n");
    AppendTo(output, " RealImagProd = ", info.eProd, "\n");
    AppendTo(output, "/\n");
    CloseStream(output);

    cmd:=Concatenation(
        "/Users/mathieudutoursikiric/GITall/GITmathieu/polyhedral_common/src_latt/TSPACE_FileFormatConversion ",
        FileNml, " GAP ", FileOut
    );
    Exec(cmd);

    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return fail;
    fi;

    tspace:=ReadAsFunction(FileOut)();

    RemoveFileIfExist(FileNml);
    RemoveFileIfExist(FileOut);

    return tspace;
end;
