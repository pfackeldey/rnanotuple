from __future__ import annotations

try:
    import ROOT
except ImportError as err:
    raise ImportError(
        """install the 'ROOT' package with:

conda install -c conda-forge root
"""
    ) from err


from tqdm import tqdm

RNTupleModel = ROOT.Experimental.RNTupleModel
RNTupleWriter = ROOT.Experimental.RNTupleWriter


def convert(filename, max_entries=None, progress=True):
    with ROOT.TFile(filename, "read") as f:
        tree = f["Events"]
        branch_names = [b.GetName() for b in tree.GetListOfBranches()]

        jaggedness_branches = [n for n in branch_names if n.startswith("n")]

        attribute_branches = [[] for _ in range(len(jaggedness_branches))]
        for b in branch_names:
            for i, j in enumerate(jaggedness_branches):
                if b.startswith(j[1:] + "_"):
                    attribute_branches[i].append(b)
                    break
        all_attribute_branches = []
        for a in attribute_branches:
            all_attribute_branches.extend(a)

        group_names = list(
            {
                b.split("_")[0]
                for b in branch_names
                if "_" in b and b not in all_attribute_branches
            }
        )

        groups = [[] for _ in range(len(group_names))]
        for b in branch_names:
            if "_" not in b or b in all_attribute_branches:
                continue
            for i, g in enumerate(group_names):
                if b.startswith(g + "_"):
                    groups[i].append(b)
                    break
        all_group_branches = []
        for g in groups:
            all_group_branches.extend(g)

        independent_branches = [
            b
            for b in branch_names
            if b not in jaggedness_branches
            and b not in all_attribute_branches
            and b not in all_group_branches
        ]

        model = RNTupleModel.Create()

        # Independent branches
        for b in tree.GetListOfBranches():
            if b.GetName() not in independent_branches:
                continue
            t = b.GetLeaf(b.GetName()).GetTypeName()
            model.MakeField[t](b.GetName())

        # Grouped branches
        for i, b in enumerate(group_names):
            struct_def = f"struct {b} {{ \n"
            for a in groups[i]:
                br = tree.GetListOfBranches()[branch_names.index(a)]
                t = br.GetLeaf(br.GetName()).GetTypeName()
                struct_def += f"  {t} "
                struct_def += f"  {a[len(b)+1:]};\n"
            struct_def += "};"
            ROOT.gInterpreter.Declare(struct_def)
            model.MakeField[b](b)

        # Jagged branches
        for i, b in enumerate(jaggedness_branches):
            struct_def = f"struct {b[1:]} {{ \n"
            for a in attribute_branches[i]:
                br = tree.GetListOfBranches()[branch_names.index(a)]
                t = br.GetLeaf(br.GetName()).GetTypeName()
                struct_def += f"  {t} "
                struct_def += f"  {a[len(b):]};\n"
            struct_def += "};"
            ROOT.gInterpreter.Declare(struct_def)
            model.MakeField["std::vector<" + b[1:] + ">"](f"{b[1:]}")

        with RNTupleWriter.Recreate(
            model, "Events", filename[:-5] + "_rntuple.root"
        ) as writer:
            n_entries = (
                tree.GetEntries()
                if max_entries is None
                else min(tree.GetEntries(), max_entries)
            )
            iterator = tqdm(range(n_entries)) if progress else range(n_entries)
            for n in iterator:
                tree.GetEntry(n)
                entry = writer.CreateEntry()

                # Independent branches
                for b in independent_branches:
                    entry[b] = getattr(tree, b)

                # Grouped branches
                for i, b in enumerate(group_names):
                    obj = entry.GetPtr(b)
                    for a in groups[i]:
                        setattr(obj, a[len(b) + 1 :], getattr(tree, a))

                # Jagged branches
                for i, b in enumerate(jaggedness_branches):
                    n_objects = getattr(tree, b)
                    if n_objects == 0:
                        continue

                    for o in range(n_objects):
                        entry.GetPtr(b[1:]).emplace_back()
                        obj = entry.GetPtr(b[1:])[o]
                        for a in attribute_branches[i]:
                            assert n_objects == len(getattr(tree, a))
                            setattr(obj, a[len(b) :], getattr(tree, a)[o])

                writer.Fill(entry)
