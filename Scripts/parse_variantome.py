import pandas as pd

variantome = pd.read_csv("/mnt/Disc4T/Projects/PhD_Project/Data/variantome.csv",
                         na_values="#VALUE!")

gene_dict = {}
with open("/mnt/Disc4T/Projects/PhD_Project/Data/gene_aliases.csv", "r+") as infile:
    for line in infile:
        linelist = line.strip().split("\t")
        newid = linelist[0]
        for oldname in linelist[1:]:
            gene_dict[oldname.upper()] = newid

newids = []
for oldname in variantome["Old_id"]:
    oldid = oldname.upper()
    if oldid in gene_dict.keys():
        newids.append(gene_dict[oldid])
    else:
        newids.append(oldid+"_oldname")


variantome["Gene_id"] = newids
variantome[variantome['Gene_id'].str.contains('_oldname')]

# Finally we need to agregate because some old ids correspond to the same newid.
final = variantome.groupby(variantome["Gene_id"]).aggregate("mean").reset_index()

final.to_csv("/mnt/Disc4T/Projects/PhD_Project/Data/variantome_parsed.csv", index=None)
