from typing import List,Dict

def slice_frames(nb_frames:int, nb_process:int) -> List[List[int]]:
    nb_frames_by_pool = int(nb_frames / nb_process)
    nb_chunks = nb_frames if nb_frames_by_pool == 0 else nb_process
    #nb_frames_by_pool = 1 if nb_frames_by_pool == 0 else nb_frames_by_pool

    intervals = []
    rest = nb_frames - nb_frames_by_pool * nb_process

    i = 0
    for j in range(rest):
        intervals.append([i, i + nb_frames_by_pool + 1])
        i = i + nb_frames_by_pool + 1
    
    for j in range(nb_chunks - rest):
        intervals.append([i, i + nb_frames_by_pool])
        i = i + nb_frames_by_pool

    return intervals


def copy_universe(frames_chunks: List[List[int]], system, mol_to_store:str) -> Dict:
    dict_store = {}
    for i,frame_int in enumerate(frames_chunks):
        dict_store[i] = {}
        copied_universe = system.copy()
        atomgroup = copied_universe.select_atoms(f"resname {mol_to_store}")
        if not atomgroup:
            raise error.MoleculeNotFound(f"Molecule {mol_to_store} doesn't exist in system")

        dict_store[i]["universe"] = copied_universe
        dict_store[i]["atomgroup"] = atomgroup
        dict_store[i]["i"] = frame_int[0]
        dict_store[i]["j"] = frame_int[1]

    return dict_store