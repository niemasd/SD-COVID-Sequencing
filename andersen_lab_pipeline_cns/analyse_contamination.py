import sys
import numpy as np
from contaminant_utils import *
from path import Path



if __name__ == "__main__":
    # grab user args
    out_pth = Path(snakemake.output)
    expt_pth = out_pth.parent
    input_pths = [Path(p) for p in snakemake.input]
    # print(input_pths)
    num_samples = len(input_pths)
    # grab data from each sample
    ans = load_all_data(input_pths)
    # generate paired reads
    ans['paired_read'] = ans.apply(lambda x: x['forward_barcode'] + '-' + x['reverse_barcode'], axis=1)
    # compute log of read counts 
    ans['log_count'] = ans['paired_read_count'].apply(lambda x: np.log(x+1))
    # generate heatmap matrix of (logged: optional) read counts per sample per paired read
    hmap, data, x, y = get_heatmap_data(ans)
    general_hmap = generate_heatmap(data, x, y)
    # generate heatmap matrix of only samples suspected of contamination
    cont_data, cont_x, cont_y = get_contaminated_data(data, hmap)
    cont_hmap = generate_heatmap(cont_data, cont_x, cont_y)
    # table of barcode read counts for contaminated samples only
    cont_table = generate_table(ans, cont_x)
    # number of contaminated samples
    num_conts = len(cont_x)
    # generate html string
    html_output = generate_html(general_hmap, cont_hmap, cont_table,
                                num_samples, num_conts, out_pth)
    # save report to file
    save_html(html_output, out_pth)
#     # grab data from each sample
#     ans = load_all_data(input_pths)
#     # generate paired reads
#     ans['paired_read'] = ans.apply(lambda x: x['forward_barcode'] + '-' + x['reverse_barcode'], axis=1)
#     # compute log of read counts 
#     ans['log_count'] = ans['paired_read_count'].apply(lambda x: np.log(x+1))
#     # generate heatmap matrix of (logged) read counts per sample per paired read
#     hmap = (ans.pivot_table(index=["sample"], columns=["paired_read"], values="paired_read_count")
#                .replace([np.inf, -np.inf], np.nan)
#                .fillna(0))
#     # drop unknown-unknown reads only
#     hmap = hmap.drop(columns='unknown-unknown')
#     # grab read counts 
#     counts = hmap.values
#     # normalize counts per sample (to address visualization issue) - IGNORE IF USING LOG
#     summed_counts = counts.sum(axis=1)[:, np.newaxis]
#     # normalize counts [OPTIONAL]
#     counts = counts / np.where(summed_counts > 0, summed_counts, 1)
#     # prepare data for identifying potential contaminants
#     flag = (ans.groupby('sample')
#                .agg(uniq_forward_bcodes = ('forward_barcode', get_unique_barcodes),
#                     uniq_reverse_bcodes = ('reverse_barcode', get_unique_barcodes)))

#     # create boolean column that identifies potential contamination
#     flag['contamination'] = flag.apply(is_contaminant, axis=1)
#     # getting only samples with potential contamination, not used later but useful to look at
#     # contaminants = flag[flag['contamination']==True]
#     # merge with original data to include the contamination flags
#     contaminants_flag = (hmap.join(flag, how='inner')['contamination']
#                              .apply(lambda x: np.where(x==True, counts.max(), counts.min())))
#     # add contaminant flag column to the read counts
#     data = np.hstack((counts, contaminants_flag[:, np.newaxis]))
#     # list of all sample IDs
#     x = hmap.index.values
#     # list of all paired reads and an extra flag column for contamination
#     y = hmap.columns.tolist() + ['contamination']
#     # heatmap of barcode read counts for all samples
#     general_hmap = generate_heatmap(data, x, y)
#     # choosing only samples with potential contaminants
#     contaminated_data = data[data[:, -1] !=0][:, :-1]
#     # list of samples with potential contamination only
#     cont_y = hmap.columns.values[~np.all(contaminated_data == 0, axis=0)].tolist()
#     # data (read counts) for contaminated samples only
#     contaminated_data = contaminated_data[:, ~np.all(contaminated_data == 0, axis=0)]
#     # list of paired reads with potential contamination only
#     cont_x = hmap.index.values[data[:, -1] !=0].tolist()
#     # heatmap of barcode read counts for contaminated samples only
#     cont_hmap = generate_heatmap(contaminated_data, cont_x, cont_y)
#     # table of barcode read counts for contaminated samples only
#     cont_table = generate_table(ans, cont_x)
#     # number of contaminated samples
#     num_conts = len(cont_x)
#     # generate html string
#     html_output = generate_html(general_hmap, cont_hmap, cont_table.to_html(),
#                                 num_samples, num_conts, expt_pth)
#     # save report to file
#     save_html(html_output, out_pth)