from itertools import product, permutations, combinations
from collections import defaultdict
from rdkit import Chem

__author__ = 'pavel'


def get_submol(mol, atom_ids):
    bond_ids = []
    for pair in combinations(atom_ids, 2):
        b = mol.GetBondBetweenAtoms(*pair)
        if b:
            bond_ids.append(b.GetIdx())
    return Chem.PathToSubmol(mol, bond_ids)


def __get_mmp_context_env(mol, radius):
    # mol is context consisting of one or more groups with single attachment point
    bond_ids = set()
    for a in mol.GetAtoms():
        if a.GetSymbol() == "*":
            i = radius
            b = Chem.FindAtomEnvironmentOfRadiusN(mol, i, a.GetIdx())
            while not b and i > 0:
                i -= 1
                b = Chem.FindAtomEnvironmentOfRadiusN(mol, i, a.GetIdx())
            bond_ids.update(b)
    return Chem.PathToSubmol(mol, list(bond_ids))


def __replace_att(mol, repl_dict):
    for a in mol.GetAtoms():
        map_num = a.GetAtomMapNum()
        if map_num in repl_dict:
            a.SetAtomMapNum(repl_dict[map_num])


def __get_maps_and_ranks(context):
    """Return list of attachment point map numbers and list of ranks (canonical SMILES without mapped attachment points)"""
    tmp_mol = Chem.Mol(context)
    maps = []
    ranks = []
    for comp in Chem.GetMolFrags(tmp_mol, asMols=True, sanitizeFrags=False):
        for a in comp.GetAtoms():
            atom_num = a.GetAtomMapNum()
            if atom_num:
                maps.append(atom_num)
                a.SetAtomMapNum(0)
                break
        ranks.append(Chem.MolToSmiles(comp))
    return maps, ranks


def __standardize_att_by_context(context, core):
    """
    Set attachment point numbers in core and context according to canonical ranks of attachment points in context
    Ties are broken
    Makes changes in place
    """
    maps, ranks = __get_maps_and_ranks(context)
    new_att = {m: i+1 for i, (r, m) in enumerate(sorted(zip(ranks, maps)))}
    __replace_att(core, new_att)
    __replace_att(context, new_att)


def __get_att_permutations(context):
    """
    Return possible permutations of attachment point map numbers as a tuple of dicts,
    where each dict: key - old number, value - new number
    """
    maps, ranks = __get_maps_and_ranks(context)

    d = defaultdict(list)
    for rank, att in zip(ranks, maps):
        d[rank].append(att)

    c = []
    for v in d.values():
        c.append([dict(zip(v, x)) for x in permutations(v, len(v))])

    return tuple(__merge_dicts(*item) for item in product(*c))


def __permute_att(mol, d):
    new_mol = Chem.Mol(mol)
    for a in new_mol.GetAtoms():
        i = a.GetAtomMapNum()
        if i in d:
            a.SetAtomMapNum(d[i])
    return new_mol


def __merge_dicts(*dicts):
    res = dicts[0].copy()
    for item in dicts[1:]:
        res.update(item)
    return res


def get_std_context_core_permutations(context, core, keep_stereo, radius, return_mols=False):
    # input are SMILES or Mol objects

    # returns standardized environment as SMILES (context with specified radius) and
    # the list of cores with permuted att. point numbers according to environment

    # if creation of core/context Mol object failed returns None (this is possible if SMILES were supplied as input)

    if isinstance(context, str):
        context = Chem.MolFromSmiles(context)
    if isinstance(core, str):
        core = Chem.MolFromSmiles(core)

    if core and context:

        att_num = len(Chem.GetMolFrags(context))

        if not keep_stereo:
            Chem.RemoveStereochemistry(context)
            Chem.RemoveStereochemistry(core)

        env = __get_mmp_context_env(context, radius)

        if att_num == 1:
            if not return_mols:
                return Chem.MolToSmiles(env, isomericSmiles=keep_stereo), \
                       [Chem.MolToSmiles(core, isomericSmiles=keep_stereo)]
            else:
                return env, [core]
        else:
            res = []
            __standardize_att_by_context(env, core)
            p = __get_att_permutations(env)
            # permute attachment point numbering only in core, since permutations in env will give the same canonical smiles
            if len(p) > 1:
                for d in p:
                    c = __permute_att(core, d)
                    res.append(c)
            else:
                res.append(core)
            if not return_mols:
                return Chem.MolToSmiles(env, isomericSmiles=keep_stereo), \
                       list(set(Chem.MolToSmiles(item, isomericSmiles=keep_stereo) for item in res))
            else:
                return env, res

    else:

        return None


def get_canon_context_core(context, core, keep_stereo, radius):
    # context and core are Mol
    # returns SMILES by default
    res = get_std_context_core_permutations(context, core, keep_stereo, radius)
    if res:
        env, cores = res
        return env, sorted(cores)[0]
    else:
        return None
