library(scales)
library(bayesboot)
library(Rmisc)

hdi.95range.bayesboot <- function(x)
{
	y = bayesboot(x, weighted.mean, use.weights=T, R = 500)
	z = data.frame(summary(y))
	r = range(z[grep("hdi", z$measure),]$value)
	r[2] - r[1]
}

ci.95.range <- function(x)
{
	r = range(CI(x))
	r[2] - r[1]
}


# Store the results of Fisher's exact test 
Fisher.Result <- setClass("Fisher.Result", slots = c(
 											p.values="numeric", 
 											min.pvalue="numeric",
 											min.pvalue.between="character"))
setMethod("show", "Fisher.Result",
          function(object) {
            cat("An object of class ", class(object), "\n", sep = "")
            cat(sprintf("Minimum p-value: %s between %s \n", signif(object@min.pvalue,3), object@min.pvalue.between))
            invisible(NULL)
          }
)

# fishers exact test for significant differential TSS use in a gene.
# expects scaled, rounded percentage of gene reads data. eg a TSS that
# got 90.14% of reads for a given gene will have an entry of 90 in the table.
# 	source("~/dev/adam/fiveprime/r/diff_dexseq.R")
fishers_exact <- function(g, use.chi=F, monte.carlo=F, verbose=F)
{
	#print(gene_df)
	comps = combn(1:ncol(g), 2)
	comp.names = combn(colnames(g), 2)
	n_compare = ncol(comps)
	if(verbose){
		if(use.chi){
		print(sprintf("Running %s pairwise ChiSquared tests", n_compare))
		}else{
			print(sprintf("Running %s pairwise exact tests", n_compare))
		}
	}
	
	rval = vector()
	for(i in 1:n_compare)
	{
		pair = unlist(comps[,i])
		conting_table = g[, pair]
		compare.name =  paste(unlist(comp.names[,i]), collapse=" vs ")
		if(verbose){print(compare.name); print(conting_table)}
		if(sum(conting_table > 0) > 0)
		{
			rval[i]=tryCatch({
			  	if(!use.chi){
			  		fisher.test(conting_table, workspace = 2e8, simulate.p.value=monte.carlo)$p.value
		  		}else{
		  			chisq.test(conting_table)$p.value
		  		}
			  	
			}, printing = function(w) {
			  	print(w)
			  	NA
			}, error = function(e) {
				print(e)
			  	NA
			})
			p = rval[i]
			if(verbose){print(p)}
		}else
		{
			print("All elements of table are 0!")
			rval[i] = NA
		}
		names(rval)[i] <- compare.name
	}
	
	return (new("Fisher.Result",
				p.values=rval,
				min.pvalue=min(rval),
				min.pvalue.between=names(rval)[which.min(rval)]))
}

# extract the scaled, rounded percentage values from 
# the main dataframe. 
# if filter > 0, drop samples with less total reads in the gene, than
# that from the returned frame

# # assume that columns are labeled in the way implied by the gene.reads.id,
# # and peak.reads.id
# get.gene <- function(diff, gene_symbol, scaled=T, filter=100, filter.on.peaks=T)
# {
# 	gene.reads.id = "_gene_reads"
# 	peak.reads.id = "_peak_reads"

# 	g.all = diff[diff$Gene==toupper(gene_symbol),]
# 	#print(g.all)
# 	rownames(g.all) = paste(g.all$Gene, g.all$PeakID, sep="_")
	
# 	if(filter.on.peaks)
# 	{
# 		peak_reads = colSums(g.all[, grep(peak.reads.id, colnames(g.all))])
# 		keep.samples = which(peak_reads > filter)
# 	}else{
# 		gene_reads = colSums(g.all[, grep(gene.reads.id, colnames(g.all))])
# 		keep.samples = which(gene_reads > filter)
# 	}	
	

# 	if(scaled)
# 	{
# 		g = g.all[, grep("scaled", colnames(g.all))]
# 		g = round(100*g)
# 		if(filter > 0)
# 		{
# 			# print(length(keep.samples))
# 			if(any(!keep.samples %in% 1:ncol(g)))
# 			{
# 				screwed.up = keep.samples[!keep.samples %in% 1:ncol(g)]
# 				error("Looking in this gene:")
# 				print(g)
# 				stop(sprintf("Could not find column %s\n", screwed.up))
# 			}
# 			if(length(keep.samples)<2) # no point returning the gene if less than 2 samples are expressed above the filter
# 			{
# 				g = NULL
# 			}else{
# 				g = g[, keep.samples]
# 			}
# 		}
# 	}else{
# 		g = g.all[, grep(peak.reads.id, colnames(g.all))]
# 		if(filter > 0)
# 		{
# 			if(length(keep.samples)<2) # no point returning the gene if less than 2 samples are expressed above the filter
# 			{
# 				g = NULL
# 			}else{
# 				g = g[, keep.samples]
# 			}
# 		}

# 	}
# 	return(g)
# }


# adapted the above function to Sean's new table. dec13, 2016
# assume that columns are labeled in the way implied by the gene.reads.id,
# and peak.reads.id
get.gene <- function(diff, gene_symbol, scaled=T, filter=0, round=T)
{
	peak.reads.id = "PR_"
	g.all = diff[toupper(diff$Gene)==toupper(gene_symbol),]
	if(nrow(g.all)==0){print(sprintf("%s not found!", toupper(gene_symbol))); return (NULL)}
	#print(g.all)
	rownames(g.all) = make.names(g.all$Name, unique=T)

	if(filter>0){
		peak_reads = colSums(g.all[, grep(peak.reads.id, colnames(g.all))])
		keep.samples = which(peak_reads > filter)
	}
	
	g = g.all[, grep(peak.reads.id, colnames(g.all))]
	if(filter > 0)
	{
		if(length(keep.samples)<2) # no point returning the gene if less than 2 samples are expressed above the filter
		{
			g = NULL
		}else{
			g = g[, keep.samples]
		}
	}
	if(scaled & !is.null(g)){

		g = 100 * sweep(g, 2, colSums(g), "/")
		if(round){g = round(g)}
	}
	return(g)
}

# removed most of the filtering
filter_diff <- function(tss_counts_matrix, keep.scaled=F)
{
	# d = tss_counts_matrix[,c("Gene", "PeakID", "Chromosome", "Start", "End", "Strand", 
	# 	"Average.Score", "Mean.Expression", grep("peak_reads|scaled", colnames(d), value=T))]
	print("Filtering TSS matrix")
	d = tss_counts_matrix
	# colnames(d) = gsub("X", "", colnames(d))

	# before = nrow(d)
	# # print(head(d))
	# d = d[d$Gene != "Intergenic", ]
	# print(sprintf("Dropped %s intergenic peak regions, %s remain [%s genes]", before-nrow(d), nrow(d), length(unique(d$Gene))))

	# #drop infs and nas
	# before = nrow(d)
	# d = na.omit(d)
	# print(sprintf("Dropped %s peak regions containing Nas or Infs, %s remain [%s genes]", before-nrow(d), nrow(d), length(unique(d$Gene))))

	#drop genes where the peak is ambiguous (i.e the genes overlap on the same strand, so its not clear which peak the gene is in)
	before = nrow(d)
	d = d[grep(",", d$Gene, invert=T),]
	print(sprintf("Dropped %s peaks in genes where the peak is ambiguous, %s remain [%s genes]", 
		before-nrow(d), nrow(d), length(unique(d$Gene))))
	
	# # filter any peak regions that get a low fraction of the gene's reads in all samples.
	# # use the scaled reads/gene columns to calculate this, but then drop those columns.
	# # print("Relying on _scaled suffix for scaled columns that was manually added. will fail if the diffTSS table has been regenerated")
	# scaled.cols = grep("scaled", colnames(d), value=T)
	# min_frac_of_gene_reads = 0.02
	# before = nrow(d)
	# d = d[ rowSums(d[,scaled.cols] > min_frac_of_gene_reads) > 0,]
	# print(sprintf("Dropped %s peak regions below %s of gene reads in all samples, %s remain [%s genes]", 
	# 	before-nrow(d), percent(min_frac_of_gene_reads), nrow(d), length(unique(d$Gene))))
	# if(!keep.scaled){
	# 	d = d[, grep("scaled", colnames(d), invert=T)]
	# }
	

	#keep only genes for which there is more than one peak.
	before = nrow(d)
	duplicated.genes = d$Gene[which(duplicated(d$Gene))];
	d = d[which(d$Gene %in% duplicated.genes),]
	print(sprintf("Dropped %s peak regions for which there is only one per gene, %s remain [%s genes]", before-nrow(d), nrow(d), length(unique(d$Gene))))

	print("Filtering complete!")
	cat("\n\n")
	return(d)
}

# significantly DE TSS
calculate_all <- function(d, filter=100, p.thresh=0.0001,NotSave=T)
{
	
	numSamps=length(grep("^PR_",colnames(d)))
	# d = filter_diff(d)
	cat("=======================================================================================\n")
	print("Searching for differential TSS use (fishers exact test)")
	cat("=======================================================================================\n\n")

	print("Filtering diffTSS matrix")
	d = filter_diff(d, keep.scaled=T)

	print(sprintf("Filter on peaks (Min reads= %s)", filter))
	genes = as.character(unique(d$Gene))

	# print("Removing NAs")
	# d[is.na(d)] <- 0
	# d = d[, grep("scaled", colnames(d))]
	# d = round(100*d)
	ngene = length(genes)
	pvals = vector()
	means = vector()
	compars = vector()
	stds = vector()
	print(sprintf("Genes: %s, p.threshold: %s", ngene, p.thresh))
	for(i in 1:ngene)
	{
		gene = unlist(genes[i])

		if(i %% 100 == 0)
		{
			print(sprintf("Testing %s [%s percent complete]", gene, 100*round(i / ngene, 3)))
		}
		# extract the scaled peaks matrix to test for differential TSS
		#print(gene)
		data = get.gene(d, gene, filter=filter)
		
		if(is.null(data))
		{
			#print(sprintf("No TSS pass filter for %s", gene))
			p = NA
		}else{
			if(nrow(data) ==0)
			{
				#print(sprintf("No TSS pass filter for %s", gene))
				p = NA
				comp = NA
			}else
			{
				res = fishers_exact(data)
				p = res@min.pvalue
				comp = res@min.pvalue.between
			}
		}
		pvals[i] = as.double(p)
		compars[i] = comp

		g = get.gene(d, gene, scaled=F, filter=filter)
		#print(g)
		if(!is.null(g))
		{
			peak_reads = colSums(g)
		}else
		{
			peak_reads = 0
		}
		means[i] = as.double(mean(peak_reads))
		stds[i] = as.double(sd(peak_reads))
		
		#print(sprintf("Samples above filter: %s, mean: %s, std.dev: %s", length(gene_reads), means[i], stds[i]))
	}
	rval = data.frame(genes, pvals, compars, means, stds)
	colnames(rval) = c("Gene", "p", "Comparison", "Mean.Expression", "Std.Dev")
	
	# realised the multiple testing was done wrong so added this but haven't tested it or 
	# looked at how it changes the results
	print("Correcting for multiple comparisons, FDR")
	numTests=numSamps*(numSamps-1)/2
	pval=as.numeric(rval[,"p"])*numTests
	pval[pval>1]=1
	print(pval[1:10])
	rval$p.adj = p.adjust(pval, method = "fdr")


	rval$CV = rval$Std.Dev/rval$Mean.Expression
	rval = rval[order(rval$p.adj), ]
	rval$LogP = -log10(rval$p.adj)
	print(head(rval))
	if(NotSave)
	{
	return(rval)
	}
	# optional filtering steps?
	# rval = rval[rval$Std.Dev < 1000,]
	#rval = rval[rval$Mean.Expression > 100,]

	g=ggplot(rval, aes(x=CV, y=LogP, label=Gene)) + 
		geom_point(size=1) + stat_density2d() + theme_bw() + 
		geom_text(data=subset(rval, LogP>15 & CV<0.8), size=1, vjust=2)
	ggsave(g, filename="diff_cv_p_plot.pdf")
	output.table = "ranked_fisher_exact_all.txt"
	output.table.sig = "ranked_fisher_exact_significant.txt"
	write.table(rval, file=output.table, sep="\t", row.names=F, quote=F)
	sig = rval[rval$p.adj<p.thresh,]
	print(sprintf("Writing %s genes significant below p=%s to %s", nrow(sig), p.thresh, output.table.sig))
	write.table(sig, file=output.table.sig, sep="\t", row.names=F, quote=F)
	
	ft = gen.full.table(d, rval)
	#ft = ft[grep(",", ft$Gene, invert=T, value=T),]
	draw.variable.tss(ft)
	ft = unsupervised(ft) # cluster the TSS and add clustering to the full table.

	write.table(ft, file="ranked_fisher_exact_full_table_clustered.txt", 
		sep="\t", quote=F, row.names=F)

	g = list("InVitro-Org"=c(rep(NA, 6), "Organoid", "In-Vitro"), 
		"FrontalLobe-InVitro"=c(rep(NA, 4), "FrontalLobe", "FrontalLobe", "In-Vitro", "In-Vitro"), 
		"AllNeurons-InVitro"=c(rep(NA, 2), rep("AllNeurons", 4), "In-Vitro", "In-Vitro"))
	supervised(ft, groups=g)

	return(rval)

}


draw.variable.tss <- function(ft, cutoffs=c(0.01, 0.05, 0.1, 100, 500), wd="variable_TSS_Figure_5")
{
	# extract variable TSS -- figure 5
	# =======================================
	init.dir = getwd()
	dir.create(wd)
	setwd(wd)
	for(mvf in cutoffs)
	{
		y = get.counts(ft, remove.ambig=T)
		ft= y$full
		s = y$scaled
		cts = y$counts
		tss.names = y$names
		
		ft$var = apply(100*s, 1, var)
		ft = ft[order(ft$var, decreasing = T),]
		print(sprintf("Variable TSS cutoff --> %s", mvf))
		if(mvf < 1)
		{
			n.var = round(nrow(ft)*mvf)
			var.tss = tss.names[1:n.var]
			id = sprintf("%spct", 100*mvf)
			print(sprintf("Extracted the %s%% most variable TSS", 100*mvf))
		}else{
			var.tss = tss.names[1:mvf]
			id = mvf
			print(sprintf("Extracted the %s most variable TSS", mvf))
		}

		maxvar.cts = cts[var.tss,]
		maxvar.s = s[var.tss,]

		colnames(maxvar.cts) = c("HBVSMC-1105","HBMEC-1005","Neuron-1525","Astrocyte-1805","Fetal frontal lobe","Adult frontal lobe","Organoids","H66 In vitro Neuron")
		colnames(maxvar.s) = colnames(maxvar.cts)

		library(NMF)
		nmf.options(grid.patch=TRUE) #avoid an extra blank page in the pdfs
	
		hmap.file = sprintf("top_%s_var_expression_TSS.pdf",id)
		hmap.file.s = sprintf("top_%s_var_relative_TSS.pdf",id)
		hmap.file.un = sprintf("top_%s_var_expression_TSS_unscaled.pdf",id)
		hmap.file.s.un = sprintf("top_%s_var_relative_TSS_unscaled.pdf",id)

		print(sprintf("Drawing them to %s", hmap.file))
		scaled.cols = rev(brewer.pal(9, "RdBu"))
		library(colorspace)
		unscaled.cols = matplotlib.viridis #ad.cubehelix(50) #"-RdYlBu2:100" #brewer.pal(9, "Blues") #colorRamps::matlab.like(10)

		cexRow=0 #and also add  to the aheatmap call
		#unscaled.cols = rev(heat_hcl(20, c = c(80, 30), l = c(30, 90), power = c(1/5, 2)))
		pdf(hmap.file, width=6, height=12); aheatmap(cexRow=cexRow,color = scaled.cols, t(scale(t(log2(maxvar.cts+1))))); dev.off()
		pdf(hmap.file.s, width=6, height=12); aheatmap(cexRow=cexRow,color = scaled.cols, t(scale(t(maxvar.s)))); dev.off()
		pdf(hmap.file.un, width=6, height=12); aheatmap( cexRow=cexRow,log2(maxvar.cts+1)); dev.off()
		pdf(hmap.file.s.un, width=6, height=12); aheatmap( cexRow=cexRow,maxvar.s); dev.off()
	}
	setwd(init.dir)
}

gen.full.table <- function(d, sig_fisher)
{
	ft = d[d$Gene %in% sig_fisher$Gene, ]
	genes = unlist(ft$Gene)
	pvals = vector()
	for(i in 1:nrow(ft))
	{
		g = as.character(genes[i])
		#print(g)
		pval = sig_fisher[sig_fisher$Gene==g, "p.adj"]
		pvals[i] = pval
	}
	ft$p.adj = pvals
	ft = ft[order(ft$p.adj), ]
	print("Calculating variability of TSS")
	print("Bayesian bootstrap high density - 95")
	ft$hdi95.range=apply(100*ft[,grep("scaled", colnames(ft))], 1, hdi.95range.bayesboot)
	print("CI95")
	ft$ci95.range= apply(100*ft[,grep("scaled", colnames(ft))], 1, ci.95.range)
	print("Variance")
	ft$var= apply(100*ft[,grep("scaled", colnames(ft))], 1, var)

	print("Renaming columns (readability)")

	ft = clean.colnames(ft, "_scaled", "S")
	ft = clean.colnames(ft, "_gene_reads", "G")
	ft = clean.colnames(ft, "_gene_inpeak_reads", "GP")
	ft = clean.colnames(ft, "_peak_reads", "P")

	inf = data.frame(c("S", "G", "GP", "P"), c("Scaled reads. Fraction of the gene's reads in that peak", 
		"Total reads in the gene", "Reads in the gene that are in peaks", "Reads in the peak"))
	colnames(inf) = c("Prefix", "Explanation")
	write.table(inf, "print.txt", sep="\t", quote=F, row.names=F)
	return(ft)
}

# utility to remove the long suffix from colnames and replace with 
# shorter prefix that I think makes the table easier to read in e.g
# excel.
clean.colnames <- function(ft, str, short.str){
	cols = grep(str, colnames(ft))
	colnames(ft)[cols] <- paste(short.str, colnames(ft)[cols], sep="_")
	colnames(ft)[cols] <- gsub(str, "", colnames(ft)[cols])
	ft
}

# provide a named list of group vectors to compare. calculate the
# TSS with the largest difference	
supervised <- function(ft, groups.list)
{
	init.dir = getwd()
	wd = "supervised"
	dir.create(wd)
	setwd(wd)
	for(i in 1:length(groups.list))
	{
		groups = groups.list[[i]]
		name = names(groups.list)[i]
		print(sprintf("Running supervised analysis %s", name))
		tab = do.supervised(ft, groups, name=name)
	}
	setwd(init.dir)
}

# run a supervised comparison for a given groups vector
do.supervised <- function(ft, groups, name=NULL) #, grouped.test="DEXSeq")
{
	init.dir = getwd()
	if(is.null(name))
	{
		name="unnamed_comparison"
	}
	library(NMF)
	nmf.options(grid.patch=TRUE) #avoid an extra blank page in the heatmap pdfs

	dir.create(name)
	setwd(name)
	print(sprintf("Running supervised analysis [%s]", name))
	#tab = ft[!ft$Hclust %in% c("Ambiguous", "NotVariable"),]
	tab = ft[grep(",", ft$Gene, invert=T),]
	genes = tab$Gene
	peaks = tab$PeakID
	rownames(tab) = paste(genes, peaks, sep="_")
	tab = tab[,grep("S_", colnames(tab))]

	colnames(tab) = c("HBVSMC-1105","HBMEC-1005","Neuron-1525","Astrocyte-1805","Fetal frontal lobe","Adult frontal lobe","Organoids","H66 In vitro Neuron")
	print(sprintf("Scaled TSS table size: %s", paste(dim(tab), collapse=",")))
	print("Comparing groups:")
	groups = factor(groups)
	print(table(groups))
	min.group.size=min(table(groups))
	if(min.group.size==1)
	{
		print("Groups have single members. Not using a test, just ranking by difference")
		diff = vector()
		log2fc = vector()
		compare = which(!is.na(groups))
		if(length(compare) != 2){
			print("For single member groups, there should be exactly 2 samples to compare")
			print(compare)
			return(FALSE)
		}
		print(sprintf("Comparing %s to %s", groups[compare[1]], groups[compare[2]]))
		for(i in 1:nrow(tab))
		{
			x=unlist(tab[i,])

			diff[i] = x[compare[1]] - x[compare[2]]
			log2fc[i] = log2(x[compare[1]]/x[compare[2]])
		}
		res = data.frame(diff, log2fc, genes, peaks)
		rownames(res) <- paste(genes, peaks, sep="_")
		res = res[order(res$diff),]
		threshold = 0.75
		sig = rownames(res[abs(res$diff) > threshold,])
		sig.up = rownames(res[res$diff > threshold,])
		sig.down = rownames(res[res$diff < threshold,])
		pdf(sprintf("%s_hmp.pdf", name)); aheatmap(tab[sig,]); dev.off()
		pdf(sprintf("%s_comparedonly_hmp.pdf", name)); aheatmap(tab[sig,compare]); dev.off()
		print(sprintf("%s TSS are differentially used (diff<%s)", sum(res$Diff), threshold))
	}else{
		p = vector()
		log2fc = vector()
		for(i in 1:nrow(tab))
		{
			z = tryCatch({
			  	t.test(unlist(tab[i,])~groups)
			}, printing = function(w) {
			  	print(w)
			  	NA
			}, error = function(e) {
				NA
			})
			if(length(z) > 1)
			{
				p[i]= z$p.value
			  	log2fc[i]=  log2(z$estimate[1] /z$estimate[2])
		  	}else{
		  		p[i]= NA
			  	log2fc[i]=  NA
		  	}
		}
		res = data.frame(p, log2fc, genes, peaks)
		rownames(res) <- paste(genes, peaks, sep="_")
		res = res[order(res$p),]
		
		y= res[!is.na(res$p),]
		sig = rownames(y[y$p<0.05,])
		sig.up = rownames(y[(y$p<0.05) & (y$log2fc>0),])
		sig.down = rownames(y[(y$p<0.05) & (y$log2fc<0),])

		draw =tab[sig,]
		pdf(sprintf("%s_hmp.pdf", name)); aheatmap(draw); dev.off()
		pdf(sprintf("%s_comparedonly_hmp.pdf", name)); aheatmap(draw[,which(!is.na(groups))]); dev.off()
		print(sprintf("%s TSS are differentially used (p<0.05)", sum(y$p<0.05)))
	}
	print("Differentially used TSS:")
	print(sig)
	write.table(res, file=sprintf("%s_table.txt", name), sep="\t", quote=F, row.names=F)
	source("~/dev/adam/rna_seq/r/go_analysis.r")

	print("Getting background list (all genes that were detected in the 5' neurons experiment")
	ct= read.delim(paste(path.expand("~"), 
		"/Desktop/projects/Fiveprime/ahaber/analysis/Differential_TSS/rerun_resequenced/correct_scaling_march'16/allPeaks/diffTSS.txt", sep="/"))
	background = as.character(unique(ct$Gene))
	background = grep("Intergenic", background, invert=T, value=T)
	print(sprintf("Got list of %s genes to use as background", length(background)))
	first.field <- function(x){strsplit(x, "_", fixed=T)[[1]][1]}
	go.analysis(unlist(lapply(sig.up, first.field)), background, sprintf("%s_Upreg", name), annotation="hg19")
	go.analysis(unlist(lapply(sig.down, first.field)), background,sprintf("%s_Downreg", name), annotation="hg19")
	setwd(init.dir)
	return(res)
}


get.counts <- function(ft, remove.ambig=F)
{
	#extract scaled (s) and counts (cts) tables from the full table (ft)
	if(remove.ambig)
	{
		before = nrow(ft)
		ft = ft[grep(",", ft$Gene, invert=T),]
		print(sprintf("Removed %s ambiguous TSS", before-nrow(ft)))
	}
	cts = ft[, grep("P_", colnames(ft))]
	cts = cts[, 1:8]
	s = ft[, grep("S_", colnames(ft))]

	tss.names = paste(ft$Gene, ft$PeakID, sep="_")
	rownames(cts) = tss.names
	rownames(s)   = tss.names
	rownames(ft)  = tss.names
	return(list("scaled"=s, "counts"=cts, "names"=tss.names, "full"=ft))
}

# run unsupervised clustering and analysis of TSS that are in genes where
# there is a significant p-value for differential TSS use. TSS with a 
# variance of 100 have an average change of 10% across the samples.
unsupervised <- function(ft=NULL, most.var.fraction = 0.1, min.var=100)
{
	if(is.null(ft))
	{
		print("Loading full table from")
		ft = read.delim("ranked_fisher_exact_full_table.txt", check.names = F)
		print("Running unsupervised analysis")
	}else
	{
		print("Running unsupervised analysis on given table")
	}

	y = get.counts(ft)
	s = y$scaled
	cts = y$counts
	tss.names = y$names
	ft$var = apply(100*s, 1, var)
	ft = ft[order(ft$var, decreasing = T),]

	# filter out TSS that don't vary before clustering.
	x = data.frame(ft$var)
	colnames(x) = "Variance"
	ambig.tss = grep(",", tss.names)

	print("Keeping only variable TSS")
	keep.tss = which(ft$var > min.var)
	keep.tss = keep.tss[!keep.tss %in% ambig.tss]
	
	print(sprintf("There are %s ambiguous TSS:", length(ambig.tss)))
	print(tss.names[ambig.tss])

	filtered.cts = cts[keep.tss,]
	filtered.s = s[keep.tss,]
	n.removed = nrow(ft) - length(ambig.tss) - length(keep.tss)
	ggsave(ggplot(x, aes(x=Variance)) + geom_histogram() + 
		theme_bw() + geom_vline(xintercept=min.var) + 
		ggtitle(sprintf("Distribution of TSS variance \n Removed %s less than %s variance, clustered the rest", 
			n.removed, min.var)), filename="variance_filtering_TSS.pdf", width=10, height=8)
 	# filtered = cts[1:2095,] # remove the 20% least variable TSS
	print("Clustering the %s remaining TSS", length(keep.tss))
	
	cluster.methods = c("Mclust", "Pam", "Hclust") #,"clusGap")
	#cluster.methods="Hclust"
	cls = list()


	print("Running Mclust")
	library(mclust)
	y = Mclust(filtered.s)
	cls[[1]] <-  y$classification

	print("Running pam")
	library(fpc)
	p = pamk(filtered.s)
	cls[[2]] <- p$pamobject$clustering
	
	print("Running hclust")
	dm   = as.dist(1-cor(t(filtered.s)))
	dend = hclust(dm)
	library(dynamicTreeCut)
	cls[[3]] <- dynamicTreeCut::cutreeDynamic(dend, 
				distM=as.matrix(dm),
				deepSplit=2,
				minClusterSize=100)

	# print("Running pam w gap statistic")
	# pam1 <- function(x,k) list(cluster = pam(x,k, cluster.only=TRUE))
	# gap = clusGap(filtered.s, FUN = pam1, K.max = 10, B = 60)
	# g = maxSE(gap$Tab[,"gap"],gap$Tab[,"SE.sim"])
	# cls[[4]] <- pam(as.dist(1-cor(t(filtered.s))), k=g)$clustering


	script("heatmap")
	print("Ambig TSS indices:")
	print(ambig.tss)
	for(i in 1:length(cluster.methods))
	{
		cl.meth = cluster.methods[i]
		cl = cls[[i]]
		print(sprintf("%s found %s clusters", cl.meth, length(unique(cl))))
		print(table(cl))
		
		fn = sprintf("%s.pdf", cl.meth)

		z = average.heatmap(t(filtered.s),  cl, NULL, cols = hmap.cols(0, 4, 9), 
			pdf.output = T, pdf.name = paste0("Relative_aggregated", fn), width=12, height=9, cexCol=0.6,  cexRow=0.6, fontsize=16)
		z = average.heatmap(t(filtered.cts),  cl, NULL, cols = hmap.cols(0, 4, 9), 
			pdf.output = T, pdf.name = paste0("Expression_aggregated", fn), width=12, height=9, cexCol=0.6, cexRow=0.6, fontsize=16)

		print("Drawing correlation matrix:")
		# pdf(paste0("TSS_filtered_correlation", fn), width=12, height=12); aheatmap(color = scaled.cols, Rowv=dend, cexCol=0, annRow=factor(as.character(cl)), cor(t(filtered.s))); dev.off()
		# pdf(paste0("TSS_filtered_all", fn), width=7, height=12); aheatmap(color = scaled.cols, Rowv=dend, annRow=factor(as.character(cl)), filtered.s); dev.off()
		ft[, cl.meth] <- "NotVariable"
		ft[keep.tss, cl.meth] <- as.character(cl)
		ft[ambig.tss, cl.meth] <- "Ambiguous"
	}
	return (ft)
}


