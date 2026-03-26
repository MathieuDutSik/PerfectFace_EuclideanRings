Read("Fct.g");




get_grp:=function()
    local ListValue, ListEntries, ListGRP, only_well_rounded, value, entry, RecTspace, GRP, desc;

    ListValue:=[-3, -4, -7, -8, -11];
    #ListValue:=[-3];
    ListEntries:=[];
    ListGRP:=[];
    for value in ListValue
    do
        RecTspace:=my_get_rec_tspace(3, value);
        entry:=rec(value:=value, rec_tspace:=RecTspace);
        GRP:=Group(entry.rec_tspace.PtStabGens);
        Add(ListGRP, GRP);
        Add(ListEntries, entry);
    od;
end;



get_rec_save_kernel:=function(d)
    local ListCells2, ListCells3, ListCells4, ListLower3, ListLower4, ListUpper2, ListUpper3, ListUpperGraph2, ListLowerGraph4, RecSave, FileSave, homology_data3_f, homology_data3_t;
    Print("Startin of process_example for d=", d, "\n");
    ListCells2:=get_cells_with_irreducibility(3, d, 6);
    ListCells3:=get_cells_with_irreducibility(3, d, 5);
    ListCells4:=get_cells_with_irreducibility(3, d, 4);

    ListLower3:=get_lower_cells(3, d, 5);
    ListLower4:=get_lower_cells(3, d, 4);

    ListUpper2:=get_upper_cells(3, d, 6);
    ListUpper3:=get_upper_cells(3, d, 5);

    ListUpperGraph2:=get_upper_graphs(3, d, 6);
    ListLowerGraph4:=get_lower_graphs(3, d, 4);

    return rec(ListCells2:=ListCells2, ListCells3:=ListCells3, ListCells4:=ListCells4,
               ListLower3:=ListLower3, ListLower4:=ListLower4,
               ListUpper2:=ListUpper2, ListUpper3:=ListUpper3,
               ListUpperGraph2:=ListUpperGraph2,
               ListLowerGraph4:=ListLowerGraph4);
end;


get_rec_save:=function(d)
    local FileSave, RecSave;
    FileSave:=Concatenation("Enumeration_3_", String(d), ".gap");
    if IsExistingFile(FileSave) then
        return ReadAsFunction(FileSave)();
    else
        RecSave:=get_rec_save_kernel(d);
        SaveDataToFile(FileSave, RecSave);
        return RecSave;
    fi;
end;

get_initial_state:=function(RecSave)
    local l_cell2, l_cell3, l_cell4;
    l_cell2:=Filtered([1..Length(RecSave.ListCells2)], x->x.is_irreducible=false);
    l_cell3:=Filtered([1..Length(RecSave.ListCells3)], x->x.is_irreducible=false);
    l_cell4:=Filtered([1..Length(RecSave.ListCells4)], x->x.is_irreducible=false);
    return rec(l_cell2:=l_cell2, l_cell3:=l_cell3, l_cell4:=l_cell4);
end;

get_missing_cells:=function(RecSave)
    local l_choice, i;
    l_choice:=[];
    for i in [1..Length(RecSave.ListCells2)]
    do
        if RecSave.ListCells2[i].is_irreducible=true then
            Add(l_choice, rec(dim:=2, index:=i));
        fi;
    od;
    for i in [1..Length(RecSave.ListCells3)]
    do
        if RecSave.ListCells3[i].is_irreducible=true then
            Add(l_choice, rec(dim:=3, index:=i));
        fi;
    od;
    for i in [1..Length(RecSave.ListCells4)]
    do
        if RecSave.ListCells4[i].is_irreducible=true then
            Add(l_choice, rec(dim:=4, index:=i));
        fi;
    od;
    return l_choice;
end;





get_subgraph:=function(g, l_status1, l_status2)
    local n_vert1, n_vert2, LVert, i;
    n_vert1:=Length(l_status1);
    n_vert2:=Length(l_status2);
    LVert:=[];
    for i in [1..n_vert1]
    do
        if l_status1[i] then
            Add(LVert, i);
        fi;
    od;
    for i in [1..n_vert2]
    do
        if l_status2[i] then
            Add(LVert, n_vert1 + i);
        fi;
    od;
    return InducedSubgraph(g, LVert);
end;


number_edges:=function(g)
    local n_entry, n_vert, i, LAdj, n_edges;
    n_entry:=0;
    n_vert:=OrderGraph(g);
    for i in [1..n_vert]
    do
        LAdj:=Adjacency(g, i);
        n_entry:=n_entry + Length(LAdj);
    od;
    n_edges:=n_entry / 2;
    return n_edges;
end;



IsTree:=function(g)
    local n_vert, n_edge;
    if IsConnectedGraph(g)=false then
        return false;
    fi;
    n_vert:=OrderGraph(g);
    n_edge:=number_edges(g);
    if n_edge<>n_vert-1 then
        return false;
    fi;
    return true;
end;


is_allowed_extension:=function(RecSave, rec_cells, choice)
    local l_cell2, l_cell3, l_cell4, RecGRA, l_status1, l_status2, h, l_lower_jOrb3, l_upper_jOrb3, l_y3, l_Xi;
    l_cell2:=rec_cells.l_cell2;
    l_cell3:=rec_cells.l_cell3;
    l_cell4:=rec_cells.l_cell4;
    if choice.dim=2 then
        # Building X_i^{> \sigma}
        RecGRA:=RecSave.ListUpperGraph2[choice.index];
        l_status1:=List(RecGRA.ListEXT1_iOrb, x->Position(x, l_cell3)<>fail);
        l_status2:=List(RecGRA.ListEXT2_iOrb, x->Position(x, l_cell4)<>fail);
        h:=get_subgraph(RecGRA.Graph, l_status1, l_status2);
        # 0-connectedness is equivalent to classic graph connectedness
        return IsConnectedGraph(h);
    fi;
    if choice.dim=3 then
        l_lower_jOrb3:=List(RecSave.ListLower3[choice.index].ListBnd, x->x.jOrb);
        l_upper_jOrb3:=List(RecSave.ListUpper3[choice.index], x->x.jOrb);
        # Y = Vor \diagdown X_i^{< \sigma}
        l_y3:=Filtered(l_lower_jOrb3, x->Position(x, l_cell2)=fail);
        # X_i^{> \sigma}
        l_Xi:=Filtered(l_upper_jOrb3, x->Position(x, l_cell4)<>fail);
        if Length(l_y3)=1 then
            return true;
        fi;
        if Length(l_Xi)=1 then
            return true;
        fi;
        if Length(l_Xi)>0 and Length(l_y3)=0 then
            return true;
        fi;
        return false;
    fi;
    if choice.dim=4 then
        # Building Y = Vor \diagdown X_i^{< \sigma}
        RecGRA:=RecSave.ListLowerGraph4[choice.index];
        l_status1:=List(RecGRA.ListEXT1_iOrb, x->Position(x, l_cell2)=fail);
        l_status2:=List(RecGRA.ListEXT1_iOrb, x->Position(x, l_cell3)=fail);
        h:=get_subgraph(RecGRA.Graph, l_status1, l_status2);
        if IsTree(h) then
            # Condition 1: Then Y is contractible
            return true;
        fi;
        if Order(h)=0 then
            # Condition 2: Y is empty
            return true;
        fi;
        return false;
    fi;
    Error("The choice is not correct");
end;



get_one_ordering:=function(RecSave, rec_cells, l_miss)
    local ret_ordering, n_miss, rec_cells_work, l_idx_done, get_one_index, append_index, i, index;
    ret_ordering:=[];
    n_miss:=Length(l_miss);
    rec_cells_work:=StructuralCopy(rec_cells);
    l_idx_done:=ListWithIdenticalEntries(n_miss, false);
    get_one_index:=function()
        local i, choice;
        for i in [1..n_miss]
        do
            if l_idx_done[i]=false then
                choice:=l_miss[i];
                if is_allowed_extension(RecSave, rec_cells, choice) then
                    return i;
                fi;
            fi;
        od;
        return fail;
    end;
    append_index:=function(idx)
        local dim, u;
        dim:=l_miss[idx].dim;
        u:=l_miss[idx].index;
        if dim=2 then
            Add(rec_cells_work.l_cell2, u);
        fi;
        if dim=3 then
            Add(rec_cells_work.l_cell3, u);
        fi;
        if dim=4 then
            Add(rec_cells_work.l_cell4, u);
        fi;
        Add(ret_ordering, l_miss[idx]);
        l_idx_done[idx]:=true;
    end;
    for i in [1..n_miss]
    do
        index:=get_one_index();
        if index=fail then
            return fail;
        fi;
        append_index(index);
    od;
    return ret_ordering;
end;








search_cell_ordering:=function(d)
    local RecSave, rec_cells, l_miss, n_miss, GRP, iter, eGen, l_miss_work, result;
    RecSave:=get_rec_save(d);
    rec_cells:=get_initial_state(RecSave);
    l_miss:=get_missing_cells(RecSave);
    n_miss:=Length(l_miss);
    GRP:=SymmetricGroup(n_miss);
    iter:=0;
    while(true)
    do
        iter:=iter+1;
        Print("Before get_one_ordering, iter=", iter, "\n");
        eGen:=Random(GRP);
        l_miss_work:=Permuted(l_miss, eGen);
        result:=get_one_ordering(RecSave, rec_cells, l_miss_work);
        if result<>fail then
            return result;
        fi;
    od;
end;







#get_rec_save(-3);
#get_rec_save(-4);
get_rec_save(-7);
get_rec_save(-8);
get_rec_save(-11);

