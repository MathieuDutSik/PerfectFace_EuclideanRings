Read("Fct.g");




get_grp:=function()
    local ListValue, ListEntries, ListGRP, value, entry, RecTspace, GRP;

    ListValue:=[-3, -4, -7, -8, -11];
    #ListValue:=[-3];
    ListEntries:=[];
    ListGRP:=[];

    for value in ListValue
    do
        RecTspace:=GetTspace(3, value);
        entry:=rec(value:=value, rec_tspace:=RecTspace);
        GRP:=Group(entry.rec_tspace.PtStabGens);
        Add(ListGRP, GRP);
        Add(ListEntries, entry);
    od;
end;

ListEXT:=get_cells(3, -3, 4);

