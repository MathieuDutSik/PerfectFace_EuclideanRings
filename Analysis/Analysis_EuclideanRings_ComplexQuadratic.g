Read("Fct.g");




get_grp:=function()
    local ListValue, ListEntries, ListGRP, only_well_rounded, value, entry, RecTspace, GRP, desc;

    ListValue:=[-3, -4, -7, -8, -11];
    #ListValue:=[-3];
    ListEntries:=[];
    ListGRP:=[];
    only_well_rounded:=true;
    for value in ListValue
    do
        desc:=GenerateTspaceDescription_imag_quad(3, value, only_well_rounded);
        RecTspace:=get_rec_tspace(desc);
        entry:=rec(value:=value, rec_tspace:=RecTspace);
        GRP:=Group(entry.rec_tspace.PtStabGens);
        Add(ListGRP, GRP);
        Add(ListEntries, entry);
    od;
end;


process_example:=function(d)
    local ListCells2, ListCells3, ListCells4, ListLower3, ListLower4, ListUpper2, ListUpper3, ListGraph2, RecSave, FileSave;
    ListCells2:=get_cells_with_irreducibility(3, d, 6);
    ListCells3:=get_cells_with_irreducibility(3, d, 5);
    ListCells4:=get_cells_with_irreducibility(3, d, 4);

    ListLower3:=get_lower_cells(3, d, 5);
    ListLower4:=get_lower_cells(3, d, 4);

    ListUpper2:=get_upper_cells(3, d, 6);
    ListUpper3:=get_upper_cells(3, d, 5);

    ListGraph2:=get_upper_graphs(3, d, 6);

    RecSave:=rec(ListCells2:=ListCells2, ListCells3:=ListCells3, ListCells4:=ListCells4,
                 ListLower3:=ListLower3, ListLower4:=ListLower4,
                 ListUpper2:=ListUpper2, ListUpper3:=ListUpper3,
                 ListGraph2:=ListGraph2);
    FileSave:=Concatenation("Enumeration_3_", String(d), ".gap");
    SaveDataToFile(FileSave, RecSave);
end;

process_example(-3);
process_example(-4);
process_example(-7);
process_example(-8);
process_example(-11);

