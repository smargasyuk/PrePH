import pandas as pd 
import click

def MakePretty(df, RNA_RNA_interaction = False):
    if RNA_RNA_interaction == True:
        print('I treat the results as RNA-RNA interaction')
    if df.shape[0] == 0:
        print('No structures found!')
    else:
        # Divide columns
        if  'gene' in list(df.columns.values):
            df[["start_gene", "end_gene"]] = df.gene.str.split('_', expand=True)[[1,2]].astype(int)

        if(RNA_RNA_interaction == False): 
            df[["strand"]] = df.interval1.str.split('_', expand=True)[[3]]
            df[["chr"]] = df.interval1.str.split('_', expand=True)[[0]]
        else:
            df[["strand_1"]] = df.interval1.str.split('_', expand=True)[[3]]
            df[["chr_1"]] = df.interval1.str.split('_', expand=True)[[0]]
            df[["strand_2"]] = df.interval2.str.split('_', expand=True)[[3]]
            df[["chr_2"]] = df.interval2.str.split('_', expand=True)[[0]]

        df[["interval1_start", "interval1_end"]] = df.interval1.str.split('_', expand=True)[[1,2]].astype(int)
        df[["interval2_start", "interval2_end"]] = df.interval2.str.split('_', expand=True)[[1,2]].astype(int)

        # Make absolute coordinates
        if(RNA_RNA_interaction == False): 
            # + strand
            df["panhandle_start"] = df.interval1_start + df.start_al1
            df["panhandle_left_hand"] = df.interval1_start + df.end_al1
            df["panhandle_right_hand"] = df.interval2_start + df.start_al2
            df["panhandle_end"] = df.interval2_start + df.end_al2
            # - strand (reverse complement)
            df.loc[df.strand == '-', 'panhandle_left_hand'] = df.loc[df.strand == '-'].interval1_end - df.loc[df.strand == '-'].start_al1
            df.loc[df.strand == '-', 'panhandle_start'] = df.loc[df.strand == '-'].interval1_end - df.loc[df.strand == '-'].end_al1
            df.loc[df.strand == '-', 'panhandle_end'] = df.loc[df.strand == '-'].interval2_end - df.loc[df.strand == '-'].start_al2
            df.loc[df.strand == '-', 'panhandle_right_hand'] = df.loc[df.strand == '-'].interval2_end - df.loc[df.strand == '-'].end_al2

            df = df.loc[~((df.panhandle_left_hand > df.panhandle_right_hand) & (df.panhandle_start < df.panhandle_end))]

            x = df.loc[df.panhandle_start > df.panhandle_end].copy()
            y = df.loc[~(df.panhandle_start > df.panhandle_end)].copy()
            if x.shape[0] != 0:
                x['al1'] = x[['alignment2']]
                x['alignment2'] = x[['alignment1']]
                x['alignment1'] = x[['al1']]
                x['lh'] = x[['panhandle_left_hand']]
                x['panhandle_left_hand'] = x[['panhandle_end']]
                x['panhandle_end'] = x[['lh']]
                x['st'] = x[['panhandle_start']]
                x['panhandle_start'] = x[['panhandle_right_hand']]
                x['panhandle_right_hand'] = x[['st']]
                x.drop(['al1', 'lh', 'st'], inplace = True, axis = 1)
                df = pd.concat([y, x])

        else:
            # + strand
            df.loc[df.strand_1 == '+',"panhandle_start"] = df.loc[df.strand_1 == '+'].interval1_start + df.loc[df.strand_1 == '+'].start_al1
            df.loc[df.strand_1 == '+',"panhandle_left_hand"] = df.loc[df.strand_1 == '+'].interval1_start + df.loc[df.strand_1 == '+'].end_al1
            df.loc[df.strand_2 == '+', "panhandle_right_hand"] = df.loc[df.strand_2 == '+'].interval2_start + df.loc[df.strand_2 == '+'].start_al2
            df.loc[df.strand_2 == '+', "panhandle_end"] = df.loc[df.strand_2 == '+'].interval2_start + df.loc[df.strand_2 == '+'].end_al2
            # - strand (reverse complement)
            df.loc[df.strand_1 == '-', 'panhandle_left_hand'] = df.loc[df.strand_1 == '-'].interval1_end - df.loc[df.strand_1 == '-'].start_al1
            df.loc[df.strand_1 == '-', 'panhandle_start'] = df.loc[df.strand_1 == '-'].interval1_end - df.loc[df.strand_1 == '-'].end_al1
            df.loc[df.strand_2 == '-', 'panhandle_end'] = df.loc[df.strand_2 == '-'].interval2_end - df.loc[df.strand_2 == '-'].start_al2
            df.loc[df.strand_2 == '-', 'panhandle_right_hand'] = df.loc[df.strand_2 == '-'].interval2_end - df.loc[df.strand_2 == '-'].end_al2



        # Calculate handle length
        df["al1_length"] = df.panhandle_left_hand - df.panhandle_start + 1
        df["al2_length"] = df.panhandle_end - df.panhandle_right_hand + 1

        # Remove broken phs (ph in one conserved interval that have end lefter than start)
        #if ((first_to_all == False) & ()):
            #df = df.loc[df.panhandle_start < df.panhandle_right_hand]
            #df = df.loc[df.panhandle_left_hand < df.panhandle_right_hand]
        df.drop(['interval1', 'interval2','start_al1','end_al1','start_al2','end_al2',
                 'interval1_start','interval2_start', 'interval1_end','interval2_end'], axis=1, inplace=True)
        if  'gene' in list(df.columns.values):
            df.drop(['gene'], axis=1, inplace=True)
        if RNA_RNA_interaction == False:
            df.sort_values(by=["chr", "panhandle_start", "panhandle_left_hand",
                                       "panhandle_right_hand", "panhandle_end"], inplace=True)
        else:
            df.sort_values(by=["chr_1", "panhandle_start", "panhandle_left_hand",
                                       "panhandle_right_hand", "panhandle_end"], inplace=True)
        df["id"] = range(1, df.shape[0] + 1)

        return df
    

@click.command()
@click.option("--input", required=True)
@click.option("--output", required=True)
def main(input, output):
    df1 = pd.read_csv(input, sep="\t")
    df2 = MakePretty(df1)
    df2.to_csv(output, sep="\t", index=False)

if __name__ == '__main__':
    main()