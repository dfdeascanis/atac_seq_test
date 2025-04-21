import argparse
import anndata as ad
import pybedtools
from pybedtools import featurefuncs
import scanpy as sc
import csv
import os
import sys
import numpy as np

def make_bed (data):
    """
    Takes as input anndata object and converts var_names in the format 'chr:start-end' to bed file format
    """
    bed_entries = []
    if isinstance(data, ad.AnnData):
        for line in data.var_names.tolist():
            chrom = line.split(':')[0] #extract the chromosome
            start, stop = line.split(':')[1].split('-') #extract the start and stop positions
            bed6 = [chrom, start, stop, '.', '.', '.'] #convert to bed6 format
            bed6 = '\t'.join(map(str, bed6)) #tab delimmited
            bed_entries.append(bed6)
        
        # Join all entries with newline characters to create a single string and store as a bedfile
        bed_string = '\n'.join(bed_entries)

        return bed_string
    
    elif isinstance(data, pybedtools.BedTool):
        for line in data:
            chrom = line.chrom
            start = line.start
            stop = line.stop
            strand = line.strand
            gene = [attr.split('"')[1] for attr in line[8].split(';') if 'gene_name' in attr][0]
            bed6 = [chrom, start, stop, gene, '.', strand]
            bed6 = '\t'.join(map(str, bed6))
            bed_entries.append(bed6)
    
        bed_string = '\n'.join(bed_entries)

        return bed_string

    else:
        raise ValueError('Input must be an anndata object or a pybedtools object')

def gtf_to_pybedtools(gtf, feature, chromosomes=None, saveas=True, ondisk=False, filename = None):
    """
    Takes as input a GTF file and converts it to a pybedtools object
    """
    if gtf is None:
        raise ValueError("Must supply GTF file")

    gtf_bedtool = pybedtools.BedTool(gtf)

    transcripts = gtf_bedtool.filter(
        lambda x: (
            x.fields[2] == feature and
            'gene_type "protein_coding"' in [field.strip() for field in x.fields[-1].split(';')] and
            (chromosomes is None or x.chrom in chromosomes)
        )
    ).sort()

    if saveas:
        if ondisk:
            if filename is None:
                raise ValueError("If saving to disk, must provide a filename")
            return transcripts.saveas(filename)
        return transcripts.saveas()
    return transcripts

def annotate_peaks_to_genes(filename,
                            gtf_file,
                            upstream=1000,
                            downstream=100,
                            distal_dist=2e4,
                            feature_type='transcript',
                            outfile='atac_peak_annotation.tsv'):
    """
    1. Annotate peaks/windows/features to overlapping Promoter regions (+ 1kb upstream and - 100bp downstream) according to 10X Genomics pipeline.
    It is possible to extend the search as needed.
    2. Annotate peaks as 'distal' regions if within a 20kb window (default) of the TSS. Upstream distances are negative and downstream distances are positive. Strandedness is also considered.

    Parameters:
    - filname (str): Path to filtered peak .h5 file from CellRanger.
    - gtf_file (str): Path to the GTF file used in CellRanger pipeline.
    - upstream (int): Number of bases upstream of the TSS start to consider for mapping. Default is 1kb.
    - downstream (int): Number of bases downstream of TSS to consider for mapping. Default is 100bp.
    - distal_dist (int): Number of bases to consider for distal peaks. Default is 20kb.
    - feature_type (str): Type of feature to search for in the GTF file (e.g., 'gene', 'transcript'). Default is 'transcript'.
    - outfile (str): Path to save the TSV file with peak annotations. Default is 'atac_peak_annotation.tsv'.

    Returns:
    - outfile (str) The function creates a peak annotation tsv file.
    """

    print(f'reading {filename} into anndata object...')

    #using scanpy to read in the 10x atac-seq data
    adata = sc.read_10x_h5(filename, gex_only=False)

    #Assuming adata.var_names.tolist() returns a list of strings like 'chr1:100-200'
    #atac seq data should be formatted as contig such as chromosme:start-end

    peaks = pybedtools.BedTool(make_bed(adata), from_string=True).filter(lambda x: x.chrom.startswith('chr')).sort()

    chromosomes = {interval.chrom for interval in peaks} #extract chromosomes from peaks

    #checks if gtf file is provided and performs filtering for genes annotated as protein coding
    #can modify in the future to filter all gene types

    print(f'reading {gtf_file} and subsetting by protein-coding {feature_type}s...')

    tmp = gtf_to_pybedtools(gtf=gtf_file, 
    feature=feature_type, 
    chromosomes=chromosomes,
    saveas=True, 
    ondisk=True, 
    filename='transcripts.gtf').sort()

    promoters = tmp.each(featurefuncs.five_prime, upstream=upstream, downstream=downstream).sort() #create a bed file of the promoters and sort
    tss = tmp.each(featurefuncs.five_prime, upstream=0, downstream=0).sort() #just the TSS

    peaks_on_promoters = pybedtools.BedTool(make_bed(promoters), from_string=True).intersect(peaks, wa=True, wb=True) #intersect the peaks with the promoters to find overlaps and write the peaks that overlap to a new bed file

    print('annotating peaks to promoters...')

    results = {}
    #rule 1 : if a peak overlaps with the promoter of any TSS, the peak is assigned a promoter label for that gene
    for i in peaks_on_promoters:
        # Extract coordinates and gene name
        feat_chrom, feat_start, feat_stop = i[-6], int(i[-5]), int(i[-4]) #coordinates of the peak
        var_idx = f'{feat_chrom}:{feat_start}-{feat_stop}' #index of the peak
        distance = [0] #create a distance of zero
        gene_name = [i[3]] #extract the gene name
        peak_type = ['promoter'] #label the peak as a promoter
        
        # Initialize the dictionary entry if not present
        if var_idx not in results:
            results[var_idx] = {'gene_names': gene_name, 'distance': distance, 'peak_type': peak_type} #creating a nested dictionary with the gene name, distance and peak type
        elif gene_name[0] not in results[var_idx]['gene_names']: #if the same peak overlaps with multiple genes, make sure to append the gene name to the list
            results[var_idx]['gene_names'].append(gene_name[0])
            results[var_idx]['distance'].append(distance[0])
            results[var_idx]['peak_type'].append(peak_type[0])

    nearby = pybedtools.BedTool(make_bed(tss), from_string=True).window(peaks, w=distal_dist, sw=True) #20kb window around the TSS

    print('annotating peaks to distal regions...')
    #rule 2 : If a peak is within 20kb of the closest TSS, and if it is not a promoter peak of the gene of the closest TSS, it will be annotated as a distal peak of that gene

    for i in nearby:

        feat_chrom, feat_start, feat_stop = i[-6], int(i[-5]), int(i[-4]) #pull the feature chrom, start, and stop
        tss_chrom, tss_start, tss_stop = i.chrom, int(i.start), int(i.stop) #pull the TSS chrom, start, and stop
        var_idx = f'{feat_chrom}:{feat_start}-{feat_stop}' #create a string index for the peak
        gene_name = [i[3]] #extract the gene name

        # Calculate the distance between the peak and the TSS
        #we need to take into account strandedness

        if i.strand == '-': #process for negative strand first
            if feat_start > tss_stop: #if the peak start is greater than the promoter stop, calculate the upstream distance
                distance = [tss_stop - feat_start]
            elif feat_stop < tss_start: #if the peak stop is less than the promoter start, calculate the downstream distance
                distance = [tss_start - feat_stop]
            
        elif i.strand == '+': 
            if feat_stop < tss_start: #if the peak stop is less than the promoter start, calculate the upstream distance
                distance = [feat_stop - tss_start]
            elif feat_start > tss_stop: #if the peak start is greater than the promoter stop, calculate the downstream distance
                distance = [feat_start - tss_stop]

        if var_idx not in results:
            results[var_idx] = {'gene_names': gene_name, 'distance': distance, 'peak_type': ['distal']} #for unseen peaks, create a new dictionary entry

        if gene_name[0] not in results[var_idx]['gene_names']: #if the feature is within the window of the promoter, but not the promoter itself, label it as distal
            results[var_idx]['gene_names'].append(gene_name[0])
            results[var_idx]['distance'].append(distance[0])
            results[var_idx]['peak_type'].append('distal')
        
        else:
            idx = results[var_idx]['gene_names'].index(gene_name[0])
            if results[var_idx]['peak_type'][idx] != 'promoter': #if the feature is already labeled as a promoter, do not append
                cdist = results[var_idx]['distance'][idx]
                if abs(distance[0]) < abs(cdist):
                    results[var_idx]['distance'][idx] = distance[0]

    gene_body = pybedtools.BedTool(make_bed(tmp), from_string=True).intersect(peaks, wa=True, wb=True) #intersect the peaks with the gene bodies

    print('Identifying peaks that overlap with gene bodies...')
    #rule 3 : If a peak overlaps with the gene body, and if it is not a promoter peak of the gene, it will be annotated as a distal peak of that gene with distance equal to 0
    for line in gene_body:
        feat_chrom, feat_start, feat_stop = line[-6], int(line[-5]), int(line[-4]) #pull the feature chrom, start, and stop
        var_idx = f'{feat_chrom}:{feat_start}-{feat_stop}' #create a string index for the peak
        gene_name =  [line[3]] #extract the gene name of the overlapping feature

        if var_idx in results:
            if gene_name[0] in results[var_idx]['gene_names']:
                idx = results[var_idx]['gene_names'].index(gene_name[0])
                if results[var_idx]['peak_type'][idx] == 'distal' and results[var_idx]['distance'][idx] != 0:
                    results[var_idx]['distance'][idx] = 0
        else:
            results[var_idx] = {'gene_names': gene_name, 'distance': [0], 'peak_type': ['distal']}

    print('annotating intergenic peaks...')
    #rule 4 : If a peak does not overlap with any gene, it will be annotated as intergenic
    for i in peaks:
        chr, start, stop = i.chrom, i.start, i.stop
        var_idx = f'{chr}:{start}-{stop}'
        if var_idx not in results:
            results[var_idx] = {'gene_names': [], 'distance' : [], 'peak_type' : ['intergenic']}
    
    #formatting the results to write to a tsv file that downstream tools (i.e., muon) can read
    for peak, data in results.items():
        gene_names = data['gene_names']
        distances = [str(d) for d in data['distance']]  # Convert each distance to a string
        peak_types = data['peak_type']

        if len(gene_names) > 1:
            results[peak]['gene_names'] = ';'.join(gene_names)
            results[peak]['distance'] = ';'.join(distances)
            results[peak]['peak_type'] = ';'.join(peak_types)
        elif peak_types[0] != 'intergenic':
            # For consistency, ensure single entries are strings and not lists
            results[peak]['gene_names'] = gene_names[0] if gene_names else ""
            results[peak]['distance'] = distances[0] if gene_names else ""
            results[peak]['peak_type'] = peak_types[0] if gene_names else ""
        else:
            results[peak]['gene_names'] = ""
            results[peak]['distance'] = ""
            results[peak]['peak_type'] = 'intergenic'

    # Open the TSV file for writing
    with open(outfile, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
            
        # Write the header row
        writer.writerow(['peak', 'gene', 'distance', 'peak_type'])
            
        # Write the data rows
        for peak, data in results.items():
            writer.writerow([
                str(peak), 
                str(data['gene_names']), 
                str(data['distance']), 
                str(data['peak_type'])
            ])

def main():

            """
            Main function to parse command-line arguments and execute the mapping of genes to regions.
            Usage in command line:
                python script_name.py --input_file path/to/adata.h5 --gtf_file path/to/gtf.gtf --upstream 1000\
                      --downstream 100 --distal_dist 2e4 --feature_type transcript --outfile atac_peak_annotation.tsv
            """
            parser = argparse.ArgumentParser(description="Map peaks to overlapping TSS regions.")
            parser.add_argument('--input_file', type=str, help="Path to the .h5 file.")
            parser.add_argument('--gtf_file', type=str, help="Path to the GTF file.")
            parser.add_argument('--upstream', type=int, default=1000, help="Number of bases upstream to consider.")
            parser.add_argument('--downstream', type=int, default=100, help="Number of bases downstream to consider.")
            parser.add_argument('--distal_dist', type=int, default=2e4, help="Number of bases to consider for distal peaks.")
            parser.add_argument('--feature_type', type=str, default='transcript', help="Feature type to search in the GTF file.")
            parser.add_argument('--outfile', type=str, default='atac_peak_annotation.tsv', help="Path to save the tsv file")

            args = parser.parse_args()

            # Check if required arguments are provided
            if not args.input_file or not args.gtf_file:
                parser.print_help()
                print(f"Error: The input_file and gtf_file are required.")
                sys.exit(1)

            # Check if files exist
            for file_path in [args.input_file, args.gtf_file]:
                if not os.path.exists(file_path):
                    print(f"Error: The file {file_path} does not exist.")
                    sys.exit(1)

            # Map genes to regions
            annotate_peaks_to_genes(
                filename=args.input_file,
                gtf_file=args.gtf_file,
                upstream=args.upstream,
                downstream=args.downstream,
                feature_type=args.feature_type,
                outfile=args.outfile)
    
if __name__ == '__main__':
    main()