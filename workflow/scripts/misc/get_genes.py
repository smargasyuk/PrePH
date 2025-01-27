import click
import polars as pl


def get_genes(in_df, output, gene_type):
    df_genes = in_df\
        .filter(pl.col("feature") == "gene")\
        .with_columns(
            (pl.col("gene_id") + "_" + pl.col("gene_name")).alias("gene_id_name")
        )
    
    if gene_type == "coding":
        print('Selected coding genes')
        df1 = df_genes.filter(
               pl.col("gene_type") == "protein_coding"
        )
    elif gene_type == "noncoding":
        print('Selected non-coding genes')
        df1 = df_genes.filter(
               pl.col("gene_type") != "protein_coding"
        )
    else:
        print('Selected all genes')
        df1 = df_genes

    df1\
        .filter(
            ~pl.col("gene_id_name").is_null() & ~pl.col("gene_id_name").str.contains(" ")
        )\
        .select(
            ["seqname", "start", "end", "gene_id_name", "score", "strand"]
        )\
        .write_csv(output, separator="\t", include_header=False)


def get_exons(in_df, output, gene_type):

    if gene_type == "coding":
        df_exons = in_df.filter(
               pl.col("feature") == "CDS"
        )
    else:
        df_exons = in_df.filter(
               pl.col("feature") == "exon"
        )

    df_exons = df_exons\
        .with_columns(
            (pl.col("gene_id") + "_" + pl.col("gene_name")).alias("gene_id_name")
        )

    df_exons\
        .select(
            ["seqname", "start", "end", "gene_id_name", "score", "strand"]
        )\
        .write_csv(output, separator="\t", include_header=False)


@click.command()
@click.option("--input", required=True)
@click.option("--output-dir", required=True)
@click.option("--gene-type", default="coding")
def main(input, output_dir, gene_type):
    df_full = pl.read_parquet(input)
    get_genes(df_full, f"{output_dir}/genes.bed", gene_type)
    get_exons(df_full, f"{output_dir}/CDS.bed", gene_type)


if __name__ == "__main__":
    main()
