import os
import h5py
import uuid
import hashlib
import numpy as np
import pandas as pd
from tardis.atomic import AtomData


old_path = os.path.abspath("kurucz_cd23_chianti_H_He.h5")
#old_path = os.path.expanduser("h_alpha_p.h5")


def single_level_selection(locnr,lines,levels):
    line = lines.loc[[locnr]].copy()
    level_idx = np.concatenate((line.level_number_lower.values , line.level_number_upper.values))
    lvls = levels.loc[ levels.level_number.isin(level_idx) ].copy()
    return line,lvls

def level_selection(locs,lines,levels,**kwargs):
    lns,lvls = single_level_selection(locs[0],lines,levels)
    for loc in locs[1:]:
        line,level = single_level_selection(loc,lines,levels)
        lns = lns.append(line)
        lvls = lvls.append(level)

    lvls.loc[:, "metastable"]  = lvls.loc[:, "metastable"].values.astype(bool)
    dtypes = zip(lvls.dtypes.index, lvls.dtypes)
    fully_ionized = np.array( [(1, 1, 0, 0.0, 1, True)], dtype=dtypes)
    lvls = lvls.drop_duplicates()
    lvls = lvls.append(pd.DataFrame(fully_ionized))
    return lns,lvls

def macro_atom_selection(lines, levels, synpp_refs, ionization_data, macro_atom_data, macro_atom_references,**kwargs):
    sel = sorted(levels.level_number.values)[1:]
    mcr_dat = macro_atom_data.loc[ macro_atom_data.source_level_number.isin(sel) 
                                   & macro_atom_data.destination_level_number.isin(sel) ].copy()
    mcr_ref = macro_atom_references.loc[ macro_atom_references.source_level_number.isin(sel) ].copy()
    
    for lev in sel:
        all,up = (mcr_dat.source_level_number == lev).sum(), ((mcr_dat.source_level_number == lev) & (mcr_dat.destination_level_number > lev)).sum() 
        mcr_ref.count_up.loc[mcr_ref.source_level_number == lev]    = up 
        mcr_ref.count_down.loc[mcr_ref.source_level_number == lev]  = all - up
        mcr_ref.count_total.loc[mcr_ref.source_level_number == lev] = all

    dtypes = zip(macro_atom_references.dtypes.index, macro_atom_references.dtypes)
    fully_ionized = np.array( [(1, 1, 0, 0, 0, 0)], dtype=dtypes)
    mcr_ref =  mcr_ref .append(pd.DataFrame(fully_ionized))
    return mcr_dat, mcr_ref


def species_selection(atom,ion,lines, levels, synpp_refs, zeta_data, atom_data, ionization_data, macro_atom_data, macro_atom_references, t_rad, data_sources, **kwargs):
    zeta_data_selected = zeta_data[atom-1,None]
    atom_data_selected = atom_data[atom-1,None]
    lines_selected  = lines.loc[(lines["atomic_number"]==atom) & (lines["ion_number"]==ion)].copy()
    lines_selected.reset_index(drop=True, inplace=True)
    levels_selected = levels.loc[(levels["atomic_number"]==atom) & (levels["ion_number"]==ion)].copy()
    levels_selected.reset_index(drop=True, inplace=True)
    synpp_refs_selected = synpp_refs.loc[(synpp_refs["atomic_number"]==atom) & (synpp_refs["ion_number"]==ion)].copy()
    ionization_data_selected = ionization_data.loc[(ionization_data["atomic_number"]==atom) & (ionization_data["ion_number"]==(ion+1))].copy()
    macro_atom_data_selected = macro_atom_data.loc[(macro_atom_data["atomic_number"]==atom) & (macro_atom_data["ion_number"]==ion)].copy()
    macro_atom_references_selected = macro_atom_references.loc[(macro_atom_references["atomic_number"]==atom) & (macro_atom_references["ion_number"]==ion)].copy()    
    return {"t_rad" : t_rad,
            "zeta_data":zeta_data_selected,
            "atom_data":atom_data_selected,
            "data_sources" : data_sources,
            "lines": lines_selected,
            "levels":levels_selected,
            "synpp_refs":synpp_refs_selected,
            "ionization_data": ionization_data_selected,
            "macro_atom_data": macro_atom_data_selected,
            "macro_atom_references": macro_atom_references_selected}

def load_from_HDFstore(old_path):
    with h5py.File(old_path) as f:
        zeta_data = f["zeta_data"][:]
        atom_data = f["basic_atom_data"][:]
        data_sources = f.attrs["data_sources"]
        t_rad = f["zeta_data"].attrs["t_rad"]

    with pd.HDFStore(old_path) as store:
        return {"t_rad" : t_rad,
                "zeta_data" : zeta_data,
                "atom_data" : atom_data,
                "data_sources" : data_sources,
                "lines" : store["lines_data"],
                "levels" : store["levels_data"],
                "synpp_refs" : store["synpp_refs"],
                "ionization_data" : store["ionization_data"],
                "macro_atom_data" : store["macro_atom_data"],
                "macro_atom_references" : store["macro_atom_references"]}

def save_to_H5(filename, zeta_data, atom_data, data_sources, lines, levels, synpp_refs, ionization_data, macro_atom_data, macro_atom_references, **kwargs):
    fpath = os.path.abspath(filename)
    with h5py.File(fpath, "w") as f:
        f["zeta_data"] = zeta_data
        f['zeta_data'].attrs['t_rad'] = np.arange(2000, 42000, 2000)
        f['zeta_data'].attrs['source'] = 'Used with kind permission from Knox Long'
        f["basic_atom_data"] = atom_data
        f["lines_data"] = lines.to_records(index=False)
        f["levels_data"] = levels.to_records(index=False)
        f["synpp_refs"] = synpp_refs.to_records(index=False)
        f["ionization_data"] = ionization_data.to_records(index=False)
        f["macro_atom_data"] = macro_atom_data.to_records(index=False)
        f["macro_atom_references"] = macro_atom_references.to_records(index=False)
        f.attrs['data_sources'] = data_sources
        f.attrs['database_version'] = 'v0.9'


        md5_hash = hashlib.md5()
        for dataset in f.values():
            md5_hash.update(dataset.value.data)
        uuid1 = uuid.uuid1().hex
        f.attrs['md5'] = md5_hash.hexdigest()
        f.attrs['uuid1'] = uuid1
    

frames = load_from_HDFstore(old_path)
frm = species_selection(1,0,**frames)
H = frm.copy()
frm["lines"],frm["levels"] = level_selection([51],**frm)
frm["macro_atom_data"],frm["macro_atom_references"] = macro_atom_selection(**frm)

print "Creating database with following lines"
print frm["lines"]

save_to_H5("test_H_15t8.h5",**frm)

