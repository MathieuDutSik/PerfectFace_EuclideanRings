Read("Fct.g");




get_grp:=function()
    local ListValue, ListEntries, ListGRP, value, entry, RecTspace, GRP;

    ListValue:=[-3, -4, -7, -8, -11];
    #ListValue:=[-3];
    ListEntries:=[];
    ListGRP:=[];

    for value in ListValue
    do
        RecTspace:=get_rec_tspace(3, value);
        entry:=rec(value:=value, rec_tspace:=RecTspace);
        GRP:=Group(entry.rec_tspace.PtStabGens);
        Add(ListGRP, GRP);
        Add(ListEntries, entry);
    od;
end;


process_example:=function(d)
    local ListCells2, ListCells3, ListCells4, ListLower3, ListLower4, ListUpper2, ListUpper3;
    ListCells2:=get_cells_with_irreducibility(3, d, 6);
    ListCells3:=get_cells_with_irreducibility(3, d, 5);
    ListCells4:=get_cells_with_irreducibility(3, d, 4);

    ListLower3:=get_lower_cells(3, d, 5);
    ListLower4:=get_lower_cells(3, d, 4);

    ListUpper2:=get_lower_cells(3, d, 6);
    ListUpper3:=get_lower_cells(3, d, 5);

    RecSave:=rec(ListCells2:=ListCells2, ListCells3:=ListCells3, ListCells4:=ListCells4,
                 ListLower3:=ListLower3, ListLower4:=ListLower4,
                 ListUpper2:=ListUpper2, ListUpper3:=ListUpper3);
    FileSave:=Concatenation("Enumeration_3_", String(d), ".gap");
    SaveDataToFile(FileSave, RecSave);
end;

process_example(-3);
#process_example(-4);
#process_example(-7);
#process_example(-8);
#process_example(-11);

