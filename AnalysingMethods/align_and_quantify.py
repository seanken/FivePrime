

def add_gene_intervals_to_rsem_output(genes_file, rsem_raw, rsem_output_file,
                                      annotation="ucsc",
                                      working_dir="."):

    if annotation != "ucsc":
        id = "gene_id"
        id_idx = 0
    else:
        id = "transcript_id"
        id_idx = 1

    with open(genes_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        known_genes = {}
        row_index = 0
        for row in reader:
            gene = {}
            gene["chr"] = row[0]
            gene["start"] = row[1]
            gene["end"] = row[2]
            gene["strand"] = row[4]
            gene["name"] = row[3]
            # map the gene name field to the gene:
            known_genes[gene["name"]] = gene

            # and also see if there are more transcript ids in another field. map them
            # all to the same gene as well:
            try:
                if ";" in row[6]:
                    v = row[6].replace(",", "") # delete the commas that bedtools merge added.
                    vals = v.split(";")
                    for a in vals:
                        if len(a) > 0:
                            key, value = a.split()
                            if key == id:
                                gene_id = str(value.replace('\"', ''))
                                if row_index < 10:
                                    print("Storing id: " + gene_id)
                                known_genes[gene_id] = gene
                else:
                    vals = row[6].split(",")
                    for a in vals:
                        known_genes[a] = gene
            except KeyError:
                continue  # nothing, this row has no attrs field
            except ValueError:
                continue
            row_index += 1

    tmp = working_dir + "/rsem.unsorted.bed"
    with open(tmp, 'wb') as out:
        writer = csv.writer(out, delimiter="\t")
        notfound = []
        with open(rsem_raw + ".genes.results", 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            row_index = 0
            for row in reader:
                gene_id = row[id_idx]
                if not gene_id.startswith(id):  # ignore the header line.
                    if "," in gene_id:
                        gene_id = gene_id.split(",")[0]#[1:]
                    try:
                        gene_model = known_genes[gene_id]
                        # bed format:
                        # chr start end name score strand
                        gene_info = [gene_model["chr"], gene_model["start"], gene_model["end"],
                                     gene_model["name"], ".", gene_model["strand"] ]
                        #print(gene_info)
                        writer.writerow(gene_info + row)
                    except KeyError:
                        if row_index < 10:
                            print("Could not find transcript id [" + gene_id + "] in ucsc genes table!")
                        notfound.append(gene_id)
                row_index += 1

    if len(notfound) > 1:
        missing = "notfound_transcripts.txt"
        with open(missing, 'w+') as f:
            for item in notfound:
                f.write("%s\n" % item)

import sys
import csv
if __name__=="__main__":
	genes_file=sys.argv[1]
	rsem_raw=sys.argv[2]
	rsem_output_file=sys.argv[2]
	direct=sys.argv[3]
	add_gene_intervals_to_rsem_output(genes_file, rsem_raw, rsem_raw,annotation="ucsc",working_dir=direct)
	
