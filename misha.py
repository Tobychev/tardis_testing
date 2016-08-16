import os
import h5py
import uuid
import hashlib
import numpy as np
import pandas as pd
from tardis.atomic import AtomData

with pd.HDFStore(old_path) as store:
    ionization_data = store["ionization_data"]
    lines = store["lines_data"]
    levels = store["levels_data"]
    macro_atom_data = store["macro_atom_data"]
    macro_atom_references = store["macro_atom_references"]
    synpp_refs = store["synpp_refs"]

lines_h0 = lines.loc[(lines["atomic_number"]==1)&(lines["ion_number"]==0)].copy()
levels_h0 = levels.loc[(levels["atomic_number"]==1)&(levels["ion_number"]==0)].copy()

macro_atom_data_h0 = macro_atom_data.loc[ (macro_atom_data["atomic_number"]==1)&(macro_atom_data["ion_number"]==0) ].copy()
macro_atom_references_h0 = macro_atom_references.loc[ (macro_atom_references["atomic_number"]==1)&(macro_atom_references["ion_number"]==0) ].copy()
ionization_data_h0 = ionization_data.loc[ (ionization_data["atomic_number"]==1)&(ionization_data["ion_number"]==1) ].copy()  
synpp_refs_h0 = synpp_refs.loc[(synpp_refs["atomic_number"]==14)&(synpp_refs["ion_number"]==1)].copy()

lines_h0 = lines_h0.reset_index(drop=True)
levels_h0 = levels_h0.reset_index(drop=True)

lines_cut = lines_h0.loc[[27]].copy()
levels_cut = levels_h0.loc[levels_h0.level_number.isin([3,8])].copy()
levels_cut.loc[:, "metastable"] = np.array([True, False], dtype=bool)
fully_ionized = [(1, 1, 0, 0.0, 1, True)]
dtypes = zip(levels.dtypes.index, levels_cut.dtypes)
fully_ionized = np.array(fully_ionized, dtype=dtypes)
levels_cut = levels_cut.append(pd.DataFrame(fully_ionized))
print lines_cut
print levels_cut
print levels_cut.dtypes
macro_atom_data_cut = macro_atom_data_h0.loc[ macro_atom_data_h0.source_level_number.isin([3,8]) 
                                            & macro_atom_data_h0.destination_level_number.isin([3,8]) ].copy()
macro_atom_references_cut = [(1, 0, 3, 0, 1, 1), (1, 0, 8, 1, 0, 2), (1, 1, 0, 0, 0, 0)]
dtypes = zip(macro_atom_references.dtypes.index, macro_atom_references.dtypes)
macro_atom_references_cut = np.array(macro_atom_references_cut, dtype=dtypes)
macro_atom_references_cut = pd.DataFrame(macro_atom_references_cut)

with h5py.File(old_path) as f:
    zeta_data = f["zeta_data"][:]
    atom_data = f["basic_atom_data"][:]
    data_sources = f.attrs["data_sources"]
    t_rad = f["zeta_data"].attrs["t_rad"]
zeta_data_h0 = zeta_data[0,None]
atom_data_h0 = atom_data[0, None]

fpath = os.path.expanduser("test_alpha.h5")
with h5py.File(fpath, "w") as f:
    f["basic_atom_data"] = atom_data_h0
    f["ionization_data"] = ionization_data_h0.to_records(index=False)
    f["lines_data"] = lines_cut.to_records(index=False)
    f["levels_data"] = levels_cut.to_records(index=False)
    f["macro_atom_data"] = macro_atom_data_cut.to_records(index=False)
    f["macro_atom_references"] = macro_atom_references_cut.to_records(index=False)
    f["synpp_refs"] = synpp_refs_h0.to_records(index=False)
    f["zeta_data"] = zeta_data_h0
    f['zeta_data'].attrs['t_rad'] = np.arange(2000, 42000, 2000)
    f['zeta_data'].attrs['source'] = 'Used with kind permission from Knox Long'
    
    f.attrs['data_sources'] = data_sources
    f.attrs['database_version'] = 'v0.9'

    md5_hash = hashlib.md5()
    for dataset in f.values():
        md5_hash.update(dataset.value.data)
    uuid1 = uuid.uuid1().hex
    f.attrs['md5'] = md5_hash.hexdigest()
    f.attrs['uuid1'] = uuid1

ad = AtomData.from_hdf5(fpath)
print ad.atom_data
print ad.ionization_data
print ad.levels
print ad.lines
print ad.macro_atom_data_all
print ad.macro_atom_references_all
print ad.synpp_refs
