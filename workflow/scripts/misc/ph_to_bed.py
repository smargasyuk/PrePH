import click
import polars as pl


def get_cut_params():
    cuts = [-35, -30, -25, -20]
    labels = [str(x) for x in range(len(cuts) + 1)]
    colors = ["255,0,255", "255,0,0", "255,100,0", "255,200,0", "0,100,0"]
    color_mapping = {k:v for k, v in zip(labels, colors)}
    return cuts, labels, color_mapping


@click.command()
@click.option("--input", required=True)
@click.option("--output", required=True)
@click.option("--header", default=None)
def main(input, output, header):
    df1 = pl.read_csv(input, separator='\t')\
        .with_columns(pl.col(x).cast(int) for x in ['panhandle_start',
        'panhandle_left_hand',
        'panhandle_right_hand',
        'panhandle_end',
        'al1_length',
        'al2_length'])
    
    cuts, labels, color_mapping = get_cut_params()

    df1bed = df1.select(
        pl.col("chr").alias("chrom"),
        pl.col("panhandle_start").alias("chromStart") - 1, 
        pl.col("panhandle_end").alias("chromEnd"),
        (pl.lit("id=") + pl.col("id").cast(str) + ",dG=" + pl.col("energy").cast(str)).alias("name"),
        pl.lit(1).alias("score"),
        pl.col("strand"),
        pl.col("panhandle_start").alias("thickStart") - 1, 
        pl.col("panhandle_end").alias("thickEnd"),   
        pl.col("energy").cut(cuts, labels=labels).replace_strict(color_mapping, return_dtype=pl.String).alias("itemRgb"),
        pl.lit(2).alias("blockCount"),
        (pl.col("al1_length").cast(str) + "," + pl.col("al2_length").cast(str)).alias("blockSizes"),
        (pl.lit("0,") + (pl.col("panhandle_right_hand") - pl.col("panhandle_start")).cast(str)).alias("blockStarts")
    )

    with open(output, "w") as f:
        if header:
            f.write(header + "\n")
        df1bed.write_csv(f, separator="\t", include_header=False)


if __name__ == "__main__":
    main()
