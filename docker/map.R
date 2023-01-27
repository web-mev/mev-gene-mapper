suppressMessages(suppressWarnings(library("optparse")))

# args from command line:
args<-commandArgs(TRUE)

option_list <- list(
    make_option(
        c('-f','--input_file'),
        help='Path to the input to remap.'
    ),
    make_option(
        c('-m','--map_file'),
        help='The file which has the mapping between identifiers.'
    ),
    make_option(
        c('-i', '--initial_id'),
        help = 'The initial identifier type'
    ),
    make_option(
        c('-t', '--target_id'),
        help = 'The target identifier type'
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

# the input file is a special-format string which includes the resource type
split_str = strsplit(opt$input_file, ':::')[[1]]
input_file_path = split_str[1]
resource_type = split_str[2]

# change the working directory to co-locate with the counts file:
working_dir <- dirname(input_file_path)
setwd(working_dir)

# the original file we wish to remap
df = read.table(input_file_path, header=T, check.names=F, row.names=1, sep='\t')

# the organism-specific mapping file. Looks like:
# "ALIAS"	"ENSEMBL"	"SYMBOL"	"REFSEQ"
# "1"	"A1B"	"ENSG00000121410"	"A1BG"	"NM_130786"
# "2"	"A1B"	"ENSG00000121410"	"A1BG"	"NP_570602"
# "3"	"A1B"	"ENSG00000172164"	"SNTB1"	"NM_021021"
# "4"	"A1B"	"ENSG00000172164"	"SNTB1"	"NP_066301"
# "5"	"A1B"	"ENSG00000172164"	"SNTB1"	"XM_011517239"
mapping = read.table(opt$map_file, header=T)

initial_id = toupper(opt$initial_id)
target_id = toupper(opt$target_id)

if (initial_id == target_id) {
    message('The initial and target gene identifiers were the same. No re-mapping will be performed.')
    quit(status=1)
}

# we want to preserve the original sample names. However, we also 
# want to use `make.names`, so we don't run into any issues with
# the dataframe ops below.
orig_names = colnames(df)
new_names = make.names(orig_names)
colnames(df) = new_names

# Subset to only include unique rows for the two systems
# we are mapping between
mapping = unique(mapping[,c(initial_id, target_id)])
mapping = na.omit(mapping)

# Perform an inner join on the rownames
merged_df = merge(df, mapping, by.x=0, by.y=initial_id)
if (dim(merged_df)[1] == 0) {
    message('The mapping database had zero identifiers in common with your data. Did you specify an incorrect gene identifier or the incorrect organism?')
    quit(status=1)
}
# merged_df looks like:
#         Row.names oncogene  pval  SYMBOL
# 1 ENSG00000121410        Y 0.010    A1BG
# 2 ENSG00000197953        Y 0.850 AADACL2
# 3 ENSG00000261846        Y 0.001 AADACL2
# 4 ENSG00000277796        N 0.050  CCL3L1
# 5 ENSG00000277796        N 0.050  CCL3L3

# Especially for mapping from ENSG to symbol, we encounter situations where 
# multiple ENSG Ids (e.g. ENSG1, ENSG2) map to a single symbol (gX). In such a case,
# the re-mapped table with rowname gX can have multiple "sources". There's no way
# around this ambiguity however. 
# Below, we do some work to find all these ambiguous mappings, determine the 
# median absolute deviation, and use the max.
# This is done under the premise that in this situation where both ENSG Ids
# map to CLN3, the top one would have more "potential to be interesting"
# since it's more varied. This is no better or worse than choosing the 
# first one we see, and users will have to use their judgement. This is why
# we provide the full gene mapping as an output.
# ENSG00000188603	219	266	207	178	203	258	233	361	274	191	343	229	272
# ENSG00000261832	1	1	3	0	1	3	3	1	0	1	0	1	2

# a boolean array which sets FALSE for the first occurence of each target_id.
# Later ones have been seen prior, so they get marked as TRUE
duplicated_targets_filter = duplicated(merged_df[target_id])

# A list of the unique target_ids (e.g. symbol) that are seen more than once
duplicated_targets = unique(merged_df[duplicated_targets_filter, target_id])

# Split the matrix to those that have repeats and those that don't
has_been_duplicated = merged_df[,target_id] %in% duplicated_targets
mm_subset = merged_df[has_been_duplicated,]

# if we have a resource type that is all numeric, we can perform a MAD calculation 
if (resource_type %in% c('MTX','I_MTX','EXP_MTX', 'RNASEQ_COUNT_MTX')){

    # to calculate the MAD, need to cast this dataframe as a matrix
    # retaining `new_names` removes the 'extra' columns related to the identifier systems we're mapping between
    mtx = as.matrix(mm_subset[new_names]) 
    mad_vals = apply(mtx,1,mad)
    mm_subset['mad'] = mad_vals

    # We next do a sort first by the ID (which effectively "groups")
    # the genes. Then the second sort on MAD sorts within each group.
    # To do the sort, we need the location of the columns
    cols = colnames(mm_subset)
    mad_col_idx = match('mad', cols)
    target_col_idx = match(target_id, cols)

    # perform the sort
    mm_subset = mm_subset[
    order(
        mm_subset[,target_col_idx], -mm_subset[,mad_col_idx]
    ), 
    ]

    # Now, use the `duplicated` method to create a boolean array
    # which effectively marks the row (within each ambiguous gene)
    # that has the maximum MAD
    is_duped = duplicated(mm_subset[target_id])
    mm_subset=mm_subset[!is_duped,]

    # Now drop the 'mad' column since we don't need it anymore
    mm_subset = subset(mm_subset, select=-c(mad))
} else {

    # if it was not a "matrix" resource type (which guarantees all numerical entries)
    # then we can't use the MAD as a proxy for the 'most interesting' gene in a multi-map
    # situation. In this case, we simply grab the first occurrence.

    # sort by the target symbol (which effectively groups everything)
    mm_subset = mm_subset[order(mm_subset[,target_id]), ]

    # Now, use the `duplicated` method to create a boolean array
    # which effectively marks the secondary occurrences, which is the records
    # we will drop
    is_duped = duplicated(mm_subset[target_id])
    mm_subset=mm_subset[!is_duped,]

}

# Now, concatenate the dataframe which contains the un-duplicated genes
# and the MAD-selected duplicated/ambiguous genes.
X = rbind(merged_df[!has_been_duplicated,], mm_subset)

# Now that the `target_id` column has only unique entries, we can
# set the rownames
rownames(X) = X[, target_id]

# Due to the failure of mapping to be 1:1/onto, we provide the mapping file as an additional output
# HOWEVER, there can be issues with validating files without unique rownames. Since the X dataframe
# above has unique rownames, we then use that to subset the original `mapping` dataframe which can
# cause issues (in WebMeV's file validation schemes) if it's directly written to file. Beyond that, 
# users would only care about where their target gene came from.
final_mapping <- X[,c('Row.names', target_id)]
colnames(final_mapping) <- c(initial_id, target_id)

# subset to keep only the columns from the original matrix (e.g. no extra mapping cols)
# and set the sample names to the originals
X = X[,new_names]
colnames(X) = orig_names

# Write the output files.
remapped_output_filename = paste(target_id, 'remapped_results.tsv', sep='_')
remapped_output_filename = paste(working_dir, remapped_output_filename, sep='/')
mapping_filename = paste(initial_id, 'to', target_id, 'mapping.tsv', sep='_')
mapping_filename = paste(working_dir, mapping_filename, sep='/')
write.table(X, remapped_output_filename, sep='\t', quote=F)
write.table(final_mapping[,c(target_id, initial_id)], mapping_filename, sep='\t', quote=F, row.names=F)

json_str = paste0(
       '{"remapped_file":{ "path": "', remapped_output_filename, 
       '","resource_type": "', resource_type, '"},',
       '"mapping":"', mapping_filename, '"}'
)
output_json <- paste(working_dir, 'outputs.json', sep='/')
write(json_str, output_json)
