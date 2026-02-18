Read("Fct.g");


ListValue:=[-3, -4, -7, -8, -11];
ListEntries:=[];

for value in ListValue
do
    RecTspace:=GetTspace(7, value);
    entry:=rec(value:=value, rec_tspace:=RecTspace);
    Add(ListEntries, entry);
od;
