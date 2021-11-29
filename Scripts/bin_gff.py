import os
import pybedtools as py
import numpy as np

def bin_region(start, stop, nbins=10):

    cuts = np.linspace(start, stop, num=(nbins+1), dtype = "int")
    bins = []

    for i in range(0,len(cuts)-1):
        bins.append((cuts[i], cuts[i+1]))

    return(bins)

def genome_to_dict(genome_file):

    chr_sizes = {}

    with open(genome_file) as infile:
        for line in infile:
            vals = line.strip().split()
            chr_sizes[vals[0]] = int(vals[1])

    return(chr_sizes)

def getGeneId(gene_bed):
    gid = gene_bed.fields[8].split(";")[0].replace("ID=", "")
    return(gid)

def elongate_and_bin_GFF(gff, genome, nbins=5):

    # Ensure output is overwritten
    prefix = "/mnt/Disc4T/Projects/PhD_Project/Data/"
    name = gff.rsplit( ".", 1 )[0].rsplit("/", 1)[1]
    sufix = "_bin"+str(nbins)+"_2prevGenes"
    outname = "".join([prefix, name, sufix])+".bed"

    print("Ouput will be written in:\n %s" % outname)

    fl = open(outname, "w+")
    fl.close()

    ## Load genome
    chr_sizes = genome_to_dict(genome)

    ## Load GFF for annotation
    ref = py.BedTool(gff)
    ref = ref.filter(lambda x: x[2] == "gene")
    ref = ref.sort()

    ## Create a variable for number of gene-bis
    gbins = nbins*3

    for gene in ref:

        ## Get gene ID

        #gene = ref[1]

        gid = getGeneId(gene)
        #print(gid)
        chrom = gene.chrom
        start = int(gene.start)+1 #Compensate a base (maybe it come from gff/bed)
        stop = int(gene.stop)

        ## Create a mini-bed for the gene and look for
        ## closest 2 genes before and after it (but not overlapping)
        linebed = "\t".join([gene.chrom, str(gene.start), str(gene.stop)])
        gene_bed = py.BedTool(linebed, from_string=True)

        pregenes = gene_bed.closest(ref, D = "ref", id = True, io = True, k = 2)
        postgenes = gene_bed.closest(ref, D = "ref", iu = True, io = True, k = 2)

        ## Resolve cases in which a gene has none or just one gene after/before it.
        ## In case there is no gene before/after:
        ## Set pre region to gene-start -2000 or start of chrom.
        ## Set post region to gene-stop + 2000 or end of chrom.

        if len(pregenes) == 1:

            pregene = pregenes[0]

            if pregene.fields[3] == ".": #No gene before
                pre_stop = max([start-2000, 0])
                pre_start = pre_stop
                prepre_start, prepre_stop = 0, 0

            else: #Just 1 gene before
                pre_start, pre_stop = pregene.fields[6], pregene.fields[7]
                prepre_start, prepre_stop = 0, 0

        else:

            pre_start, prepre_start = pregenes[0].fields[6], pregenes[1].fields[6]
            pre_stop, prepre_stop = pregenes[0].fields[7], pregenes[1].fields[7]

        if len(postgenes) == 1:

            postgene = postgenes[0]

            if postgene.fields[3] == ".": #No gene after
                post_start = min([stop+2000, chr_sizes[chrom]])
                post_stop = post_start
                postpost_start, postpost_stop = 0, 0

            else: #Just 1 gene after
                post_start, post_stop = postgene.fields[6], postgene.fields[7]
                postpost_start, postpost_stop = 0, 0

        else:

            post_stop, postpost_stop = postgenes[0].fields[7], postgenes[1].fields[7]
            post_start, postpost_start = postgenes[0].fields[6], postgenes[1].fields[6]


        ## Create bins

        prepre_gen_cov = bin_region(int(prepre_start), int(prepre_stop), nbins=2)
        pre_gen_cov = bin_region(int(pre_start), int(pre_stop), nbins=2)

        postpost_gen_cov = bin_region(int(postpost_start), int(postpost_stop), nbins=2)
        post_gen_cov = bin_region(int(post_start), int(post_stop), nbins=2)

        pre = bin_region(int(pre_stop), start, nbins=nbins)
        body = bin_region(start, stop, nbins=nbins)
        post = bin_region(stop, int(post_start), nbins=nbins)

        ## Print output taking into accound strandness
        ## We cannot sort here (reverse if strand is "-")
        ## because we will tabix afterwards.
        ## We create a column for sorting afterwards.

        output = []

        regions = [prepre_gen_cov,
                   pre_gen_cov,
                   pre,
                   body,
                   post,
                   postpost_gen_cov,
                   post_gen_cov]

        if gene.strand == "-":
            order = range(gbins+8, 0, -1)
        else:
            order = range(1, gbins+9)

        i = 0
        for reg in regions:
            for interval in reg:
                output.append([chrom,
                               str(interval[0]),
                               str(interval[1]),
                               gid, str(order[i])])
                i += 1

        for line in output:
            with open(outname, "a+") as outfile:
                outfile.write("\t".join(line)+"\n")


#gff = "/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-45_Pfalciparum3D7.gff"
gff = "/mnt/Disc4T/Projects/PhD_Project/Data/PlasmoDB-46_Pfalciparum3D7_withGDV1_ncRNAs.gff"
gen = "/home/lucas/ISGlobal/Gen_Referencies/Pf3D7.genome"

elongate_and_bin_GFF(gff, gen)
