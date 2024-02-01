#!/usr/bin/env python3
import sys, os, gzip, argparse, datetime

# (c) Angel G. Rivera-Colon

PROG = sys.argv[0].split('/')[-1]
DESC = '''Subset an annotation GTF file to extract annotations with a specified
distance from a set of target outlier SNPs.'''

class GtfElement:
    def __init__(self, gene_id: str, chromosome: str, start_bp: int, end_bp: str,
                 direction: str, annotation_type: str, gene_name: str ='NA'):
        assert direction in ['+','-']
        self.id    = gene_id
        self.chr   = chromosome
        self.start = start_bp
        self.end   = end_bp
        self.dir   = direction
        self.name  = gene_name
        self.type  = annotation_type
        self.snps  = list()
    def add_snp(self, snp_id):
        self.snps.append(snp_id)
    def print_row(self, delimeter=','):
        snps_str = ';'.join(self.snps)
        row_list = [self.chr, self.type, str(self.start), str(self.end), self.dir, self.id, str(self.name), snps_str]
        row = delimeter.join(row_list)
        return f'{row}\n'
    def __str__(self):
        # Chr   Type   StartBP   EndBP   Dir   GeneID   GeneName   SNPS
        row = f'{self.chr}\t{self.type}\t{self.start}\t{self.end}\t{self.dir}\t{self.id}\t{self.name}\t{self.snps}'
        return row

def parse_args():
    p = argparse.ArgumentParser(prog=PROG, description=DESC)
    p.add_argument('-g', '--gtf', required=True, help='(str) Path to annotation GTF')
    p.add_argument('-l', '--outlier-list', required=True, help='(str) Path to outlier list CSV')
    p.add_argument('-o', '--outdir', required=False, default='.', help='(str) Path to output directory')
    p.add_argument('-b', '--buffer', required=False, type=int, default=10_000, help='(int) Distance in bp plus/minus target SNP to extract a gene [default=10000]')
    args = p.parse_args()
    if not os.path.exists(args.gtf):
        sys.exit(f'Error: `{args.gtf}`: GTF does not exist.')
    if not os.path.exists(args.outlier_list):
        sys.exit(f'Error: `{args.outlier_list}`: outlier list does not exist.')
    if not os.path.exists(args.outdir):
        sys.exit(f'Error: `{args.outdir}`: output directory does not exist.')
    if args.buffer < 0:
        sys.exit(f'Error: buffer distance ({args.buffer:,}) must be >= 0.')
    return args

def now():
    '''Print the current date and time.'''
    return f'{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'

def write_output_csv(gtf_elements: list, per_gene_snps: dict, outdir: str='.'):
    '''Write the output CSV containing the retained GTF elements and corresponding SNP IDs'''
    print('Writing output...')
    # Tallies
    n_records = 0
    n_genes = 0
    n_exons = 0
    # Open file handle
    fh = open(f'{outdir}/candidate_annotations.csv', 'w')
    fh.write('#Chrom,Type,StartBP,EndBP,Dir,GeneID,GeneName,CandidateSNP\n')
    # Loop over the retained elements
    for element in gtf_elements:
        assert isinstance(element, GtfElement)
        n_records += 1
        if element.type == 'gene':
            n_genes += 1
        elif element.type == 'exon':
            n_exons += 1
        # Add SNPs to the annotation element
        element.snps = per_gene_snps[element.id]
        # Save to file
        fh.write(element.print_row())
    fh.close()
    # Print tallies
    print(f'    After filtering, kept {n_records:,} records, containing {n_genes:,} genes and {n_exons:,} exons.\n')


def parse_gtf(gtf_f: str, outliers: dict,  outdir: str, buffer: int=10_000):
    '''Parse an input GTF file and extract the annotations close to the target SNPs.
    Save the annotations to a file.'''
    print('Parsing GTF...')
    assert os.path.exists(gtf_f)
    # Temp data keeping gene ID and gene name pairs
    gene_names = dict()
    # Temp data storing the SNPs in each gene ID
    gene_snps = dict()
    # Temp data storing the GtfElements kept
    gtf_elements = list()
    # Line tallies
    n_records = 0
    n_genes = 0
    n_exons = 0
    # Check input
    infh = None
    if gtf_f.endswith('gz'):
        infh = gzip.open(gtf_f, 'rt')
    else:
        infh = open(gtf_f, 'r')
    with infh:
        for i, line in enumerate(infh):
            gtf_element = None
            if line.startswith('#'):
                continue
            fields = line.strip('\n').split('\t')
            # Check format
            assert len(fields) == 9, 'Error: GTF line must contain 9 fields.'
            n_records += 1
            # Skip unwanted annotation types
            annotation_type = fields[2]
            if annotation_type not in ['gene', 'exon']:
                continue
            chrom = fields[0]
            start_bp = int(fields[3])
            end_bp = int(fields[4])
            direction = fields[6]
            gene_id = None
            gene = None
            gene_name = None
            description = None
            attributes = fields[8]
            # Process Gene annotations
            if annotation_type == 'gene':
                n_genes += 1
                # Filter the annotations based on the outlier SNPs
                if chrom not in outliers:
                    # Skip if the chromosome has no outliers
                    continue
                # Check possible SNPs in the range of the gene
                snps_in_gene = list()
                for snp in sorted(outliers[chrom]):
                    snp_range = range((snp-buffer), (snp+buffer+1))
                    gene_range = range(start_bp, (end_bp+1))
                    # Skip
                    # Find an overlap
                    overlap = set(snp_range).intersection(gene_range)
                    if len(overlap) > 0:
                        snps_in_gene.append(f'{chrom}:{snp}')
                # Skip if gene has no nearby SNPs
                if len(snps_in_gene) == 0:
                    continue
                # Process the attributes line of the GTF
                attributes = attributes.split(';')
                for attribute in attributes:
                    attribute = attribute.lstrip(' ')
                    if len(attribute) == 0:
                        continue
                    if attribute.startswith('gene_id'):
                        gene_id = attribute.split(' ')[1]
                        gene_id = gene_id.strip('"')
                    elif attribute.startswith('gene_name'):
                        gene_name = attribute[10:]
                        gene_name = gene_name.strip('"')
                        gene_name = gene_name.strip('\'')
                    elif attribute.startswith('gene '):
                        gene_name = attribute[5:]
                        gene_name = gene_name.strip('"')
                        gene_name = gene_name.strip('\'')
                    elif attribute.startswith('description'):
                        description = attribute[12:].strip('"')
                # Match the gene ID and gene names
                if gene_name is None:
                    if gene is not None:
                        gene_name = gene
                    elif description is not None:
                        gene_name = description
                    else:
                        gene_name = gene_id
                gene_names[gene_id] = gene_name
                # Match the gene ID with the SNPs
                gene_snps[gene_id] = snps_in_gene
                # Store the annotation
                gtf_element = GtfElement(gene_id, chrom, start_bp, end_bp, direction, annotation_type, gene_name)
                gtf_elements.append(gtf_element)
            # Process Exon annotations
            elif annotation_type == 'exon':
                n_exons += 1
                # Filter the annotations based on the outlier SNPs
                if chrom not in outliers:
                    # Skip if the chromosome has no outliers
                    continue
                # Process the attributes line of the GTF
                attributes = attributes.split(';')
                for attribute in attributes:
                    attribute = attribute.lstrip(' ')
                    if len(attribute) == 0:
                        continue
                    if attribute.startswith('gene_id'):
                        gene_id = attribute.split(' ')[1]
                        gene_id = gene_id.strip('"')
                # Skip if this exon doesn't have a matching gene
                gene_name = gene_names.get(gene_id, None)
                if gene_name is None:
                    continue
                # Store the annotation
                gtf_element = GtfElement(gene_id, chrom, start_bp, end_bp, direction, annotation_type, gene_name)
                gtf_elements.append(gtf_element)
            # if i > 1_000:
            #     break
            # print(gtf_element)
    # Report GTF stats
    print(f'    Parsed {n_records:,} records from the GTF, containing {n_genes:,} genes and {n_exons:,} exons.\n')
    # Generate output
    write_output_csv(gtf_elements, gene_snps, outdir)


def parse_outliers(outlier_list_f: str) -> dict:
    '''Parse the outlier list file'''
    print('Parsing outlier list...')
    n_snps = 0
    # Main output:
    # { chr_id : [ bp_1, bp_2, ..., bp_n ] }
    outliers = dict()
    with open(outlier_list_f) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.strip('\n').strip('\r').split(',') # Also removes the \r line endings
            assert len(fields) == 2, 'Error: outlier list must be a CSV with two columns, chromosome ID and basepair'
            n_snps += 1
            chromosome = fields[0]
            basepair = int(fields[1])
            # Populate dictionary
            outliers.setdefault(chromosome, [])
            outliers[chromosome].append(basepair)
    print(f'    Loaded {n_snps:,} outlier SNPs.\n')
    return outliers


def main():
    # Report start
    print(f'{PROG} started on {now()}\n')
    args = parse_args()
    outliers = parse_outliers(args.outlier_list)
    parse_gtf(args.gtf, outliers, args.outdir, args.buffer)
    # End report
    print(f'{PROG} finished on {now()}')

# Run Code
if __name__ == '__main__':
    main()
