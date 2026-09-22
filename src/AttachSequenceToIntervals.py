from Bio.Seq import Seq
import click
from functools import partial
import polars as pl
from pyfaidx import Fasta


def GetSequencesForDF(genome, row):
    row_seq = (str(genome[row['chrom']][row['chromStart'] - 1:row['chromEnd']]).upper())
    if row['strand'] == '+':
        return row_seq
    else:
        return str(Seq(row_seq).reverse_complement())


@click.command()
@click.option("--input", required=True)
@click.option("--fasta", required=True)
@click.option("--output", required=True)
def main(input, fasta, output):
    df1 = pl.read_csv(input, separator='\t')
    genome = Fasta(fasta)
    GetSequencesForDF2 = partial(GetSequencesForDF, genome)
    df2 = df1.with_columns(
        pl.Series("sequence", map(GetSequencesForDF2, df1.iter_rows(named=True)))
    )
    df2.write_csv(output, separator='\t')

if __name__ == '__main__':
    main()