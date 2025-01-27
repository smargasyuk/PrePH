import click
import polars as pl


@click.command()
@click.option("--input", required=True)
@click.option("--input-genes", required=True)
@click.option("--output", required=True)
def main(input, input_genes, output):
    df_ci = pl.read_csv(input, separator="\t", has_header=False, new_columns=['chrom', 'chromStart', 'chromEnd', 'gene'])\
        .with_columns(pl.col("gene").str.split(";"))\
        .explode(["gene"])
    
    df_genes = pl.read_csv(input_genes, separator="\t", has_header=False, new_columns=['chrom', 'start_gene', 'end_gene', 'gene', 'score', 'strand'])

    dfm = df_ci.join(
        df_genes.select(['start_gene', 'end_gene', 'gene', 'strand']),
        on='gene'
    ).unique(
        subset = ['chrom', 'chromStart', 'chromEnd', 'strand', 'start_gene', 'end_gene'], keep = 'first'
    )\
    .with_row_index(offset=0, name="name")\
    .with_columns(
        pl.lit(1).alias("score")
    )\
    .sort(by =['chrom', 'chromStart', 'chromEnd'])\
    
    dfm.select(
        ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'start_gene', 'end_gene']
    ).write_csv(output, separator="\t")

    print('selected intervals')


if __name__ == "__main__":
    main()
